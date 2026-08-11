#!/usr/bin/env python3
"""Parse an eukannot samplesheet, infer SE/PE layout, validate, and emit a
fully-resolved TSV plus a JSON channel description.

The resolved TSV is the reproducible record of what the pipeline actually ran
with, and can be fed straight back in as --input.

Layout resolution ladder (first match wins), per the design:
  R1 explicit `layout` column always wins, but is cross-checked -- a sheet that
     disagrees with itself fails loudly instead of being "helpfully" corrected.
  R2 fastq_2 present               -> PE
  R3 fastq_2 absent                -> filesystem mate inference on fastq_1
  R4 no candidate                  -> SE
  R5 layout=interleaved            -> never inferred, only declared
  R6 comma lists are technical replicates; PE lists must be equal length
  R7 a `bam` row takes its layout from the BAM header/flags
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import sys

NULLS = {"", "-", ".", "na", "n/a", "null", "none"}

# Accepted header spellings -> canonical name.
ALIASES = {
    "sample": "sample_id", "sample_id": "sample_id", "id": "sample_id", "name": "sample_id",
    "layout": "layout", "library_layout": "layout", "type": "layout",
    "fastq_1": "fastq_1", "fastq1": "fastq_1", "r1": "fastq_1", "read1": "fastq_1", "fq1": "fastq_1",
    "fastq_2": "fastq_2", "fastq2": "fastq_2", "r2": "fastq_2", "read2": "fastq_2", "fq2": "fastq_2",
    "strandedness": "strandedness", "strand": "strandedness",
    "condition": "condition", "group": "condition",
    "replicate": "replicate", "rep": "replicate",
    "read_type": "read_type",
    "bam": "bam", "alignment": "bam",
    "genome_id": "genome_id", "genome": "genome_id",
    "use_for": "use_for",
}

SAMPLE_RE = re.compile(r"^[A-Za-z][A-Za-z0-9._-]{0,63}$")
FQ_EXT_RE = re.compile(r"\.(fastq|fq)(\.(gz|bz2|zst))?$", re.IGNORECASE)

# (pattern, replacement) applied to the extension-stripped stem, in order.
MATE_RULES = [
    (re.compile(r"_R1_(\d+)$"), r"_R2_\1"),
    (re.compile(r"_R1$"), "_R2"),
    (re.compile(r"\.R1$"), ".R2"),
    (re.compile(r"_1$"), "_2"),
    (re.compile(r"\.1$"), ".2"),
    (re.compile(r"_fwd$"), "_rev"),
    (re.compile(r"_f$"), "_r"),
    (re.compile(r"-1$"), "-2"),
]

VALID_STRAND = {"auto", "unstranded", "forward", "reverse"}
VALID_READ_TYPE = {"short", "isoseq", "ont_cdna", "ont_drna"}
VALID_USE_FOR = {"braker", "stringtie", "count", "all"}


class SheetError(Exception):
    pass


def die(msg: str) -> "None":
    sys.stderr.write(f"ERROR ~ {msg}\n")
    sys.exit(1)


def warn(msg: str) -> None:
    sys.stderr.write(f"WARN  ~ {msg}\n")


def info(msg: str) -> None:
    sys.stderr.write(f"INFO  ~ {msg}\n")


def is_null(value: str) -> bool:
    return value.strip().lower() in NULLS


def read_rows(path: str):
    """Yield (line_number, [fields]) skipping comments and blank lines.

    The separator is taken from the extension, and the technical-replicate
    separator follows it: ',' for TSV, ';' for CSV.
    """
    sep = "," if path.lower().endswith(".csv") else "\t"
    rep_sep = ";" if sep == "," else ","
    with open(path, newline="") as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.rstrip("\r\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            yield lineno, [c.strip() for c in line.split(sep)], rep_sep


def parse_header(fields, lineno):
    canon, unknown = [], []
    for col in fields:
        key = col.strip().lower().lstrip("#").strip()
        if key in ALIASES:
            canon.append(ALIASES[key])
        else:
            canon.append(None)
            unknown.append(col)
    if "sample_id" not in canon:
        die(
            f"[samplesheet] line {lineno}: header is missing a sample id column.\n"
            f"        Found: {','.join(fields)}\n"
            f"        Accepted aliases: sample|sample_id|id|name ; "
            f"fastq_1|fastq1|R1|read1|fq1 ; fastq_2|fastq2|R2|read2|fq2"
        )
    if unknown:
        warn(f"[samplesheet] ignoring unrecognised column(s): {', '.join(unknown)}")
    return canon


def resolve_path(value: str, raw_dir: str | None) -> str:
    if os.path.isabs(value) or raw_dir is None:
        return value
    return os.path.join(raw_dir, value)


def first_read_ids(path: str, n: int = 4):
    """Return the first n read identifiers, mate suffixes stripped."""
    opener = gzip.open if path.endswith(".gz") else open
    ids = []
    try:
        with opener(path, "rt", errors="replace") as fh:  # type: ignore[operator]
            for i, line in enumerate(fh):
                if i % 4 != 0:
                    continue
                token = line[1:].split()[0] if len(line) > 1 else ""
                ids.append(re.sub(r"[/._][12]$", "", token))
                if len(ids) >= n:
                    break
    except OSError:
        return []
    return ids


def infer_mate(fastq1: str):
    """Return (mate_path, rule_index) or (None, reason)."""
    base = os.path.basename(fastq1)
    ext_match = FQ_EXT_RE.search(base)
    if not ext_match:
        return None, "filename does not look like a fastq"
    stem = base[: ext_match.start()]
    ext = base[ext_match.start():]
    directory = os.path.dirname(fastq1)

    for idx, (pattern, replacement) in enumerate(MATE_RULES):
        if not pattern.search(stem):
            continue
        # Guard: an ambiguous stem (e.g. sample_R1_lane_R1) is not inferable.
        if len(pattern.findall(stem)) > 1:
            return None, f"stem {stem!r} matches rule {idx} more than once"
        candidate_stem = pattern.sub(replacement, stem)
        if candidate_stem == stem:
            continue
        candidate = os.path.join(directory, candidate_stem + ext)
        if candidate == fastq1:
            continue
        if os.path.isfile(candidate) and os.path.getsize(candidate) > 0:
            return candidate, idx
    return None, "no mate candidate found"


def check_pair_ids(fastq1: str, fastq2: str, sample: str, lineno: int) -> None:
    ids1, ids2 = first_read_ids(fastq1), first_read_ids(fastq2)
    if not ids1 or not ids2:
        return  # unreadable or empty: alignment will complain more usefully
    if ids1 != ids2:
        die(
            f"[samplesheet] line {lineno}, sample {sample!r}: inferred mate "
            f"{fastq2} but read IDs disagree\n"
            f"        R1: {ids1[0]!r}   R2: {ids2[0]!r}\n"
            f"        Declare fastq_2 explicitly if this pairing is intended."
        )


def bam_is_paired(path: str) -> bool | None:
    try:
        import pysam  # type: ignore
    except ImportError:
        return None
    try:
        with pysam.AlignmentFile(path, "rb", check_sq=False) as handle:
            for i, record in enumerate(handle.fetch(until_eof=True)):
                if i >= 10000:
                    break
                if record.is_paired:
                    return True
        return False
    except Exception:  # noqa: BLE001 - any pysam failure means "unknown"
        return None


def bam_sort_order(path: str) -> str:
    try:
        import pysam  # type: ignore
    except ImportError:
        return "unknown"
    try:
        with pysam.AlignmentFile(path, "rb", check_sq=False) as handle:
            return handle.header.to_dict().get("HD", {}).get("SO", "unknown")
    except Exception:  # noqa: BLE001
        return "unknown"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True)
    ap.add_argument("--raw-dir", default=None)
    ap.add_argument("--out-tsv", default="samplesheet.resolved.tsv")
    ap.add_argument("--out-json", default="samplesheet.resolved.json")
    ap.add_argument("--default-strandedness", default="auto")
    ap.add_argument("--default-genome-id", default=None)
    ap.add_argument("--sample", default=None, help="keep only this sample_id")
    ap.add_argument("--known-genome-ids", default="", help="comma list for validation")
    ap.add_argument("--check-files", action="store_true", default=True)
    ap.add_argument("--no-check-files", dest="check_files", action="store_false")
    args = ap.parse_args()

    if not os.path.isfile(args.input):
        die(f"[samplesheet] file not found: {args.input}")

    known_genomes = {g for g in args.known_genome_ids.split(",") if g}
    header = None
    rows, seen = [], {}

    for lineno, fields, rep_sep in read_rows(args.input):
        if header is None:
            header = parse_header(fields, lineno)
            continue

        record = {}
        for name, value in zip(header, fields):
            if name:
                record[name] = value
        sample = record.get("sample_id", "").strip()
        if not sample:
            die(f"[samplesheet] line {lineno}: empty sample_id")
        if not SAMPLE_RE.match(sample):
            die(
                f"[samplesheet] line {lineno}: sample_id {sample!r} is invalid.\n"
                f"        Must match [A-Za-z][A-Za-z0-9._-]{{0,63}}"
            )
        if sample in seen:
            die(
                f"[samplesheet] line {lineno}: sample_id {sample!r} is duplicated "
                f"(first seen at line {seen[sample]}).\n"
                f"        Sample IDs must be unique; use the 'replicate' column "
                f"for technical replicates."
            )
        seen[sample] = lineno

        if args.sample and sample != args.sample:
            continue

        layout = record.get("layout", "").strip().lower()
        if layout in NULLS:
            layout = "auto"
        if layout not in {"auto", "se", "pe", "interleaved"}:
            die(
                f"[samplesheet] line {lineno}, sample {sample!r}: layout={layout!r} "
                f"is not one of SE, PE, interleaved, auto."
            )

        bam_raw = record.get("bam", "").strip()
        fq1_raw = record.get("fastq_1", "").strip()
        fq2_raw = record.get("fastq_2", "").strip()
        has_bam = not is_null(bam_raw)
        has_fq1 = not is_null(fq1_raw)
        has_fq2 = not is_null(fq2_raw)

        if has_bam and has_fq1:
            die(
                f"[samplesheet] line {lineno}, sample {sample!r}: both 'fastq_1' "
                f"and 'bam' are set. These are mutually exclusive; use 'bam' for "
                f"pre-aligned input."
            )
        if not has_bam and not has_fq1:
            die(
                f"[samplesheet] line {lineno}, sample {sample!r}: neither 'fastq_1' "
                f"nor 'bam' is set. One is required."
            )

        # ---------------------------------------------------------- R7: BAM
        if has_bam:
            bam = resolve_path(bam_raw, args.raw_dir)
            if args.check_files and not os.path.isfile(bam):
                die(
                    f"[samplesheet] line {lineno}, sample {sample!r}: bam not found.\n"
                    f"        Tried: {bam}"
                )
            if has_fq2:
                die(
                    f"[samplesheet] line {lineno}, sample {sample!r}: 'bam' is set "
                    f"so fastq columns must be '-'."
                )
            if args.check_files:
                order = bam_sort_order(bam)
                if order == "queryname":
                    die(
                        f"[samplesheet] line {lineno}, sample {sample!r}: bam is "
                        f"name-sorted (SO:queryname); the pipeline requires "
                        f"coordinate-sorted BAM.\n"
                        f"        fix: samtools sort -o sorted.bam {bam}"
                    )
                if not (os.path.isfile(bam + ".csi") or os.path.isfile(bam + ".bai")):
                    die(
                        f"[samplesheet] line {lineno}, sample {sample!r}: bam has no "
                        f"index.\n        fix: samtools index -c {bam}"
                    )
                paired = bam_is_paired(bam)
            else:
                paired = None
            resolved_layout = "PE" if paired else ("SE" if paired is False else "SE")
            files1, files2 = [], []
        else:
            bam = ""
            files1 = [resolve_path(p.strip(), args.raw_dir)
                      for p in fq1_raw.split(rep_sep) if p.strip()]
            files2 = [resolve_path(p.strip(), args.raw_dir)
                      for p in fq2_raw.split(rep_sep) if p.strip()] if has_fq2 else []

            if args.check_files:
                for path in files1 + files2:
                    if not os.path.isfile(path):
                        die(
                            f"[samplesheet] line {lineno}, sample {sample!r}: fastq "
                            f"not found.\n        Tried: {path}"
                            + (f"  (resolved against --raw_dir={args.raw_dir})"
                               if args.raw_dir else "")
                            + "\n        Give an absolute path, or fix --raw_dir."
                        )

            # ----------------------------------------------- R1: explicit wins
            if layout == "pe" and not has_fq2:
                die(
                    f"[samplesheet] line {lineno}, sample {sample!r}: layout=PE but "
                    f"fastq_2 is '-'.\n        Either give a mate file, or set "
                    f"layout to SE/auto."
                )
            if layout == "se" and has_fq2:
                die(
                    f"[samplesheet] line {lineno}, sample {sample!r}: layout=SE but "
                    f"fastq_2={fq2_raw!r} exists.\n        Refusing to silently "
                    f"discard a mate file. Set layout=PE or layout=auto."
                )
            if layout == "interleaved" and has_fq2:
                die(
                    f"[samplesheet] line {lineno}, sample {sample!r}: "
                    f"layout=interleaved but fastq_2 is also set."
                )

            if layout == "interleaved":
                resolved_layout = "interleaved"          # R5
            elif has_fq2:
                resolved_layout = "PE"                   # R2
            elif layout == "se":
                resolved_layout = "SE"
            else:
                # ------------------------------------- R3: filesystem inference
                mates, rules = [], []
                for path in files1:
                    mate, why = infer_mate(path)
                    if mate is None:
                        mates = []
                        if layout == "auto" and why.startswith("stem"):
                            warn(
                                f"[layout] sample {sample!r}: {why}; not inferring a "
                                f"mate. Declare fastq_2 explicitly if this is PE."
                            )
                        break
                    mates.append(mate)
                    rules.append(why)
                if mates and len(mates) == len(files1):
                    if args.check_files:
                        check_pair_ids(files1[0], mates[0], sample, lineno)
                    files2 = mates
                    resolved_layout = "PE"
                    info(
                        f"[layout] {sample}: inferred PE mate -> {mates[0]} "
                        f"(rule {rules[0]})"
                    )
                else:
                    resolved_layout = "SE"               # R4

            # -------------------------------- R6: technical replicate lengths
            if resolved_layout == "PE" and len(files1) != len(files2):
                die(
                    f"[samplesheet] line {lineno}, sample {sample!r}: fastq_1 has "
                    f"{len(files1)} file(s) but fastq_2 has {len(files2)}.\n"
                    f"        Technical-replicate lists must be the same length."
                )

        strandedness = record.get("strandedness", "").strip().lower()
        if strandedness in NULLS:
            strandedness = args.default_strandedness
        if strandedness not in VALID_STRAND:
            die(
                f"[samplesheet] line {lineno}, sample {sample!r}: strandedness="
                f"{strandedness!r} is not one of {sorted(VALID_STRAND)}."
            )

        read_type = record.get("read_type", "").strip().lower()
        if read_type in NULLS:
            read_type = "short"
        if read_type not in VALID_READ_TYPE:
            die(
                f"[samplesheet] line {lineno}, sample {sample!r}: read_type="
                f"{read_type!r} is not one of {sorted(VALID_READ_TYPE)}."
            )

        use_for = record.get("use_for", "").strip().lower()
        if use_for in NULLS:
            use_for = "all"
        use_set = {u.strip() for u in use_for.split(",") if u.strip()}
        bad = use_set - VALID_USE_FOR
        if bad:
            die(
                f"[samplesheet] line {lineno}, sample {sample!r}: use_for contains "
                f"{sorted(bad)}; allowed: {sorted(VALID_USE_FOR)}."
            )
        if "all" in use_set:
            use_set = {"braker", "stringtie", "count"}

        genome_id = record.get("genome_id", "").strip()
        if is_null(genome_id):
            genome_id = args.default_genome_id or ""
        if known_genomes and genome_id and genome_id not in known_genomes:
            die(
                f"[samplesheet] line {lineno}, sample {sample!r}: genome_id "
                f"{genome_id!r} not found in the genomes table.\n"
                f"        Known ids: {', '.join(sorted(known_genomes))}"
            )

        replicate = record.get("replicate", "").strip()
        if is_null(replicate):
            replicate = "1"

        rows.append({
            "sample_id": sample,
            "layout": resolved_layout,
            "fastq_1": files1,
            "fastq_2": files2,
            "bam": bam,
            "strandedness": strandedness,
            "condition": (record.get("condition", "").strip() or "NA"),
            "replicate": replicate,
            "read_type": read_type,
            "genome_id": genome_id,
            "use_for": sorted(use_set),
        })

    if header is None:
        die(f"[samplesheet] {args.input} has no header row")
    if not rows:
        target = f" matching --sample {args.sample}" if args.sample else ""
        die(f"[samplesheet] no usable rows found in {args.input}{target}")

    n_iso = sum(1 for r in rows if r["read_type"] != "short")
    if n_iso:
        info(f"[samplesheet] {n_iso} long-read row(s) detected")

    cols = ["sample_id", "layout", "fastq_1", "fastq_2", "strandedness",
            "condition", "replicate", "read_type", "bam", "genome_id", "use_for"]
    with open(args.out_tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join([
                row["sample_id"], row["layout"],
                ",".join(row["fastq_1"]) or "-",
                ",".join(row["fastq_2"]) or "-",
                row["strandedness"], row["condition"], row["replicate"],
                row["read_type"], row["bam"] or "-", row["genome_id"] or "-",
                ",".join(row["use_for"]),
            ]) + "\n")

    with open(args.out_json, "w") as fh:
        json.dump(rows, fh, indent=2)

    n_pe = sum(1 for r in rows if r["layout"] == "PE")
    n_se = sum(1 for r in rows if r["layout"] == "SE")
    n_bam = sum(1 for r in rows if r["bam"])
    info(f"[samplesheet] {len(rows)} sample(s): {n_pe} PE, {n_se} SE, {n_bam} pre-aligned")
    return 0


if __name__ == "__main__":
    sys.exit(main())
