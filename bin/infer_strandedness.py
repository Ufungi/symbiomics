#!/usr/bin/env python3
"""Annotation-free strandedness inference from splice-junction motifs.

Why not Salmon or RSeQC: both need a transcriptome or a BED12 gene model, and
at this point in the pipeline there is no annotation -- producing one is the
reason the pipeline is running. Junction motifs need only the genome, so this
is the method that generalises.

How it works: for every gapped alignment, read the 2 bp at each end of the
intron out of the genome. GT..AG means the transcript runs on the plus strand,
CT..AC means minus (the same intron read from the other side). Compare that to
the strand the read itself aligned to; the fraction that agree is the library's
strandedness.

Bands, and why 'ambiguous' is separate from 'unstranded': a fraction of 0.72 is
not a weakly-stranded library, it is a signal that something is wrong (mixed
prep, gDNA carry-over, degraded RNA). Silently calling that 'forward' is worse
than refusing to call it.
"""
from __future__ import annotations

import argparse
import json
import sys

# Canonical and minor-spliceosome donor/acceptor pairs, and the same intron
# read from the opposite strand.
PLUS_MOTIFS = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
MINUS_MOTIFS = {("CT", "AC"), ("CT", "GC"), ("GT", "AT")}

BAM_CREF_SKIP = 3  # 'N' in the CIGAR


class Genome:
    """Random access to an indexed FASTA without loading it into memory."""

    def __init__(self, fasta: str):
        try:
            import pysam  # type: ignore
        except ImportError:
            sys.stderr.write("ERROR ~ [strandedness] pysam is required\n")
            raise
        self.fh = pysam.FastaFile(fasta)

    def dinuc(self, contig: str, start: int, end: int) -> str:
        # Uppercase before comparing: the genome is softmasked, and lowercase
        # repeats would otherwise never match a motif.
        try:
            return self.fh.fetch(contig, start, end).upper()
        except (KeyError, ValueError):
            return ""


def classify(fraction: float, informative: int, min_reads: int, threshold: float):
    if informative < min_reads:
        return "undetermined", (
            f"only {informative} informative spliced reads (need {min_reads}); "
            f"too few junctions to call"
        )
    if fraction >= threshold:
        return "forward", f"{fraction:.3f} of reads agree with the transcript strand"
    if fraction <= 1 - threshold:
        return "reverse", f"{fraction:.3f} of reads agree with the transcript strand"
    if 0.40 <= fraction <= 0.60:
        return "unstranded", f"{fraction:.3f} is consistent with an unstranded library"
    return "ambiguous", (
        f"{fraction:.3f} falls between the stranded and unstranded bands -- "
        f"possible mixed library prep, gDNA carry-over, or degraded RNA"
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--sample-id", default="sample")
    ap.add_argument("--out", default="strandedness.json")
    ap.add_argument("--max-reads", type=int, default=2_000_000)
    ap.add_argument("--min-reads", type=int, default=20_000)
    ap.add_argument("--threshold", type=float, default=0.85)
    ap.add_argument("--declared", default="auto",
                    help="value from the samplesheet; it wins, but disagreement is reported")
    ap.add_argument("--strict", action="store_true")
    args = ap.parse_args()

    import pysam  # type: ignore

    genome = Genome(args.fasta)
    agree = disagree = 0
    seen = spliced = 0
    discarded_conflict = 0
    discarded_motif = 0

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for record in bam.fetch(until_eof=True):
            if seen >= args.max_reads:
                break
            if record.is_unmapped or record.is_secondary or record.is_supplementary \
                    or record.is_duplicate:
                continue
            seen += 1
            # For pairs, read 1 carries the fragment's orientation.
            if record.is_paired and not record.is_read1:
                continue
            cigar = record.cigartuples
            if not cigar or not any(op == BAM_CREF_SKIP for op, _ in cigar):
                continue

            contig = record.reference_name
            pos = record.reference_start
            strands = set()
            for op, length in cigar:
                if op in (0, 2, 7, 8):          # M, D, =, X consume reference
                    pos += length
                elif op == BAM_CREF_SKIP:       # N: this is an intron
                    donor = genome.dinuc(contig, pos, pos + 2)
                    acceptor = genome.dinuc(contig, pos + length - 2, pos + length)
                    pair = (donor, acceptor)
                    if pair in PLUS_MOTIFS:
                        strands.add("+")
                    elif pair in MINUS_MOTIFS:
                        strands.add("-")
                    else:
                        discarded_motif += 1
                    pos += length

            if not strands:
                continue
            if len(strands) > 1:
                # A read whose junctions disagree tells us nothing.
                discarded_conflict += 1
                continue

            spliced += 1
            transcript_strand = strands.pop()
            read_strand = "-" if record.is_reverse else "+"
            if read_strand == transcript_strand:
                agree += 1
            else:
                disagree += 1

    informative = agree + disagree
    fraction = agree / informative if informative else 0.0
    call, reason = classify(fraction, informative, args.min_reads, args.threshold)

    declared = (args.declared or "auto").lower()
    effective = call
    conflict = False
    if declared not in ("auto", "", "-"):
        effective = declared
        if call in ("forward", "reverse", "unstranded") and call != declared:
            conflict = True
    elif call in ("undetermined", "ambiguous"):
        effective = "unstranded"

    result = {
        "sample_id": args.sample_id,
        "inferred": call,
        "declared": declared,
        "effective": effective,
        "conflict": conflict,
        "fraction_agreeing": round(fraction, 4),
        "informative_reads": informative,
        "spliced_reads": spliced,
        "reads_examined": seen,
        "discarded_conflicting_junctions": discarded_conflict,
        "discarded_noncanonical_motifs": discarded_motif,
        "reason": reason,
    }
    with open(args.out, "w") as fh:
        json.dump(result, fh, indent=2)

    sys.stderr.write(
        f"INFO  ~ [strandedness] {args.sample_id}: inferred={call} "
        f"f={fraction:.3f} n={informative} -> using {effective}\n"
    )
    if conflict:
        message = (
            f"[strandedness] {args.sample_id}: declared {declared!r} but inferred "
            f"{call!r} (f={fraction:.3f}, n={informative}). Using the declared "
            f"value. Re-check the library prep."
        )
        if args.strict:
            sys.stderr.write(f"ERROR ~ {message}\n")
            return 2
        sys.stderr.write(f"WARN  ~ {message}\n")
    if call in ("undetermined", "ambiguous"):
        message = f"[strandedness] {args.sample_id}: {reason}; falling back to {effective!r}"
        if args.strict:
            sys.stderr.write(f"ERROR ~ {message}\n")
            return 2
        sys.stderr.write(f"WARN  ~ {message}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
