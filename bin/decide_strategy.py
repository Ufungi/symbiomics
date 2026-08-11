#!/usr/bin/env python3
"""The policy engine.

Turns (genome stats, taxon, available evidence, resources, user budget) into an
explicit execution plan, written as strategy.yml.

Three properties matter more than the individual rules:

  1. Every decision carries a reason, the alternatives that were rejected, and
     the escape hatches. This is not a black box.
  2. --dry-run produces the plan without running anything, so a multi-week job
     can be reviewed before it is launched.
  3. Tracks carry time budgets. A track that blows its budget is dropped and
     the pipeline finishes with what is left; the omission is recorded here and
     surfaced in the report. Partial failure is not total failure.

Explicit user parameters always win over inference; when they do, the reason
records that.
"""
from __future__ import annotations

import argparse
import json
import sys

SIZE_ORDER = ["tiny", "small", "medium", "large", "huge"]

# Tools with no trained model / declared support for a clade. Recorded so the
# exclusion is documented rather than rediscovered. See docs/clade_notes.md.
CLADE_UNSUPPORTED = {
    "tiberius": {
        "gymnosperm": "Tiberius ships Angiosperms (Mesangiospermae); no gymnosperm model",
    },
    "egapx": {
        "gymnosperm": "NCBI EGAPx supports Magnoliopsida; gymnosperms are out of scope",
        "fungi": "NCBI EGAPx explicitly excludes fungi",
    },
}

DEFAULTS = {
    "masker_by_class": {
        "tiny":   {"any": "red"},
        "small":  {"any": "earlgrey"},
        "medium": {"plant": "edta", "any": "earlgrey"},
        "large":  {"any": "repeatmodeler_subsample"},
        "huge":   {"any": "red"},
    },
    "masker_blocked": {
        "large": ["earlgrey", "edta"],
        "huge":  ["earlgrey", "edta"],
    },
    "predictor_primary": {
        "tiny": "braker", "small": "braker", "medium": "braker",
        "large": "braker", "huge": "helixer",
    },
    "predictor_secondary": {"huge": "braker_sharded"},
    "predictor_blocked": {
        "large": ["funannotate2"],
        "huge":  ["funannotate2", "galba"],
    },
    "star_blocked_above": 2_000_000_000,
    "hisat2_large_index_above": 4_000_000_000,
    "hisat2_mm_index_above_gb": 8,
    "bai_max_contig": 536_870_912,
    "align_maxforks": {"tiny": 16, "small": 8, "medium": 6, "large": 4, "huge": 4},
    "braker_shard_threshold": 3_000_000_000,
    "braker_shard_mb": {"large": 20, "huge": 10},
    "braker_shard_overlap_kb": {"large": 200, "huge": 500},
    "disable_above": {
        "infernal": 2_000_000_000,
        "qualimap": 2_000_000_000,
        "tesorter": 2_000_000_000,
    },
    "track_budget": {
        "braker": "240.h", "braker_sharded": "240.h", "helixer": "96.h",
        "galba": "96.h", "funannotate2": "96.h", "repeats": "336.h",
    },
    "busco_lineage": {
        "fungi": "fungi_odb12", "plant": "embryophyta_odb12",
        "animal": "metazoa_odb12", "other": "eukaryota_odb12",
    },
    "odb_partition": {
        "fungi": "Fungi", "plant": "Viridiplantae",
        "animal": "Metazoa", "other": "Eukaryota",
    },
}


def decision(choice, reason, alternatives=None, escape=None, source="policy"):
    out = {"choice": choice, "reason": reason, "source": source}
    if alternatives:
        out["alternatives"] = alternatives
    if escape:
        out["escape"] = escape
    return out


def classify(total_bp: int, thresholds: dict) -> str:
    if total_bp <= thresholds["tiny"]:
        return "tiny"
    if total_bp <= thresholds["small"]:
        return "small"
    if total_bp <= thresholds["medium"]:
        return "medium"
    if total_bp <= thresholds["large"]:
        return "large"
    return "huge"


def pick_by_taxon(table: dict, taxon: str):
    if taxon in table:
        return table[taxon]
    return table.get("any")


def yaml_dump(obj, indent=0) -> str:
    """Minimal YAML writer -- avoids a PyYAML dependency in a tiny helper."""
    pad = "  " * indent
    lines = []
    if isinstance(obj, dict):
        if not obj:
            return pad + "{}\n"
        for key, value in obj.items():
            if isinstance(value, (dict, list)) and value:
                lines.append(f"{pad}{key}:")
                lines.append(yaml_dump(value, indent + 1).rstrip("\n"))
            else:
                lines.append(f"{pad}{key}: {scalar(value)}")
    elif isinstance(obj, list):
        if not obj:
            return pad + "[]\n"
        for item in obj:
            if isinstance(item, (dict, list)):
                nested = yaml_dump(item, indent + 1).rstrip("\n")
                lines.append(f"{pad}-")
                lines.append(nested)
            else:
                lines.append(f"{pad}- {scalar(item)}")
    else:
        return pad + scalar(obj) + "\n"
    return "\n".join(lines) + "\n"


def scalar(value) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    text = str(value)
    if text == "" or any(ch in text for ch in ":#{}[],&*?|<>=!%@`\"'") or text.strip() != text:
        return json.dumps(text)
    return text


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome-stats", required=True, help="genome_stats.json")
    ap.add_argument("--policy", default=None, help="policy overrides as JSON")
    ap.add_argument("--out", default="strategy.yml")
    ap.add_argument("--out-json", default="strategy.json")

    ap.add_argument("--genome-id", default="genome")
    ap.add_argument("--species", default=None)
    ap.add_argument("--taxon", default="auto")
    ap.add_argument("--clade", default=None,
                    help="finer clade hint, e.g. gymnosperm; gates clade-limited tools")
    ap.add_argument("--size-class", default="auto")
    ap.add_argument("--threshold-tiny", type=int, default=10_000_000)
    ap.add_argument("--threshold-small", type=int, default=500_000_000)
    ap.add_argument("--threshold-medium", type=int, default=2_000_000_000)
    ap.add_argument("--threshold-large", type=int, default=5_000_000_000)

    ap.add_argument("--masker", default="auto")
    ap.add_argument("--premasked", default="false")
    ap.add_argument("--repeat-lib", default=None)
    ap.add_argument("--aligner", default="auto")
    ap.add_argument("--hisat2-index", default=None)
    ap.add_argument("--hisat2-index-bytes", type=int, default=0)
    ap.add_argument("--bam-index", default="auto")

    ap.add_argument("--has-rna", default="false")
    ap.add_argument("--has-protein", default="false")
    ap.add_argument("--has-isoseq", default="false")
    ap.add_argument("--n-samples", type=int, default=0)

    ap.add_argument("--priority", default="balanced")
    ap.add_argument("--time-budget", default=None)
    ap.add_argument("--gpu", default="false")
    ap.add_argument("--max-cpus", type=int, default=0)
    ap.add_argument("--max-memory-gb", type=int, default=0)
    args = ap.parse_args()

    def flag(value: str) -> bool:
        return str(value).strip().lower() in {"true", "1", "yes", "on"}

    policy = dict(DEFAULTS)
    if args.policy:
        try:
            with open(args.policy) as fh:
                policy.update(json.load(fh))
        except (OSError, ValueError) as exc:
            sys.stderr.write(f"WARN  ~ [policy] ignoring --policy: {exc}\n")

    with open(args.genome_stats) as fh:
        stats = json.load(fh)

    total_bp = int(stats["total_bp"])
    max_contig = int(stats["max_contig_bp"])
    softmask_frac = float(stats.get("softmask_fraction", 0.0))

    notes: list[str] = []
    decisions: dict = {}

    # ------------------------------------------------------------- size class
    thresholds = {
        "tiny": args.threshold_tiny, "small": args.threshold_small,
        "medium": args.threshold_medium, "large": args.threshold_large,
    }
    if args.size_class in SIZE_ORDER:
        size_class = args.size_class
        size_source = "user"
    else:
        size_class = classify(total_bp, thresholds)
        size_source = "policy"
    decisions["size_class"] = decision(
        size_class,
        f"{total_bp/1e9:.3f} Gb against thresholds "
        f"tiny<={thresholds['tiny']/1e6:.0f}Mb small<={thresholds['small']/1e6:.0f}Mb "
        f"medium<={thresholds['medium']/1e9:.0f}Gb large<={thresholds['large']/1e9:.0f}Gb",
        source=size_source,
    )

    # ------------------------------------------------------------------ taxon
    taxon = args.taxon if args.taxon != "auto" else "other"
    taxon_source = "user" if args.taxon != "auto" else "default"
    if args.taxon == "auto":
        notes.append(
            "taxon could not be inferred (no --species_taxid); defaulting to 'other'. "
            "Set --taxon to enable clade-specific tracks."
        )
    decisions["taxon"] = decision(taxon, "from --taxon" if taxon_source == "user"
                                  else "no taxid given", source=taxon_source)

    clade = args.clade or ("gymnosperm" if taxon == "plant" and
                           (args.species or "").lower().startswith("pinus") else None)

    # ----------------------------------------------------------------- masker
    blocked_maskers = policy["masker_blocked"].get(size_class, [])
    premasked = flag(args.premasked)

    # The lowercase test catches the failure that actually bites -- a genome
    # declared premasked that carries no mask at all. It cannot prove the
    # opposite: some assemblers emit lowercase of their own, so a high fraction
    # is necessary but not sufficient evidence that RepeatMasker ever ran.
    if premasked:
        if softmask_frac < 0.02:
            sys.stderr.write(
                f"ERROR ~ [repeats] --premasked was given but {args.genome_id} is "
                f"{softmask_frac*100:.2f}% lowercase, i.e. not softmasked.\n"
                f"        Soft-masking is case-based; a genome with no lowercase "
                f"carries no mask.\n"
                f"        fix: drop --premasked, or supply a genuinely softmasked "
                f"genome (RepeatMasker -xsmall).\n"
            )
            return 2
        masker_choice = "none"
        masker_reason = f"--premasked and genome is {softmask_frac*100:.1f}% lowercase"
        masker_source = "user"
    elif args.repeat_lib:
        masker_choice = "repeatmasker_only"
        masker_reason = f"--repeat_lib given ({args.repeat_lib}); skipping de-novo modelling"
        masker_source = "user"
    elif args.masker != "auto":
        if args.masker in blocked_maskers:
            sys.stderr.write(
                f"ERROR ~ [repeats] masker={args.masker!r} cannot complete on a "
                f"{size_class} genome ({total_bp/1e9:.2f} Gb).\n"
                f"        EarlGrey's own guidance is 'weeks' at 25 Gb, and a "
                f"RepeatModeler run on a 21.7 Gb conifer genome burned 52 h and "
                f"emitted an empty library.\n"
                f"        alternatives: red (hours, no library), "
                f"repeatmodeler_subsample (model on a sampled fraction), "
                f"--repeat_lib <fasta>, --premasked\n"
            )
            return 2
        masker_choice = args.masker
        masker_reason = "explicitly requested"
        masker_source = "user"
    else:
        masker_choice = pick_by_taxon(policy["masker_by_class"][size_class], taxon)
        masker_reason = f"size_class={size_class}, taxon={taxon}"
        if size_class in ("large", "huge"):
            masker_reason += ("; de-novo modelling does not finish at this scale "
                              "(measured: 52 h to an empty library on a 21.7 Gb conifer)")
        masker_source = "policy"

    decisions["masker"] = decision(
        masker_choice, masker_reason,
        alternatives=[m for m in ["red", "repeatmodeler_subsample", "earlgrey", "edta"]
                      if m != masker_choice and m not in blocked_maskers],
        escape=["--repeat_lib <fasta>", "--premasked", "--masker <name>"],
        source=masker_source,
    )
    if blocked_maskers:
        decisions["masker"]["blocked_at_this_size"] = blocked_maskers

    # ---------------------------------------------------------------- aligner
    star_blocked = total_bp > policy["star_blocked_above"]
    if args.aligner == "star" and star_blocked:
        sys.stderr.write(
            f"ERROR ~ [align] aligner=star requested but the genome is "
            f"{total_bp/1e9:.2f} Gb.\n"
            f"        STAR's suffix-array index is not buildable at this scale on "
            f"any reasonable host.\n"
            f"        fix: --aligner hisat2 (uses --large-index automatically)\n"
        )
        return 2
    if args.aligner != "auto":
        aligner, aligner_reason, aligner_source = args.aligner, "explicitly requested", "user"
    elif star_blocked:
        aligner = "hisat2"
        aligner_reason = (f"genome {total_bp/1e9:.2f} Gb exceeds the STAR ceiling "
                          f"({policy['star_blocked_above']/1e9:.0f} Gb)")
        aligner_source = "policy"
    else:
        aligner = "hisat2"
        aligner_reason = "default short-read aligner"
        aligner_source = "policy"

    large_index = total_bp > policy["hisat2_large_index_above"]
    index_gb = args.hisat2_index_bytes / 1e9 if args.hisat2_index_bytes else 0
    use_mm = index_gb > policy["hisat2_mm_index_above_gb"] or large_index

    decisions["aligner"] = decision(
        aligner, aligner_reason,
        alternatives=[] if star_blocked else ["star"],
        source=aligner_source,
    )
    decisions["aligner"]["large_index"] = large_index
    decisions["aligner"]["mmap_index"] = use_mm
    if use_mm:
        decisions["aligner"]["mmap_reason"] = (
            "index is large; --mm shares one page-cache copy across concurrent "
            "aligners instead of one resident copy each"
        )
    decisions["aligner"]["reuse_index"] = args.hisat2_index or None
    decisions["aligner"]["maxforks"] = policy["align_maxforks"].get(size_class, 8)

    # -------------------------------------------------------------- BAM index
    needs_csi = max_contig > policy["bai_max_contig"]
    if args.bam_index != "auto":
        if args.bam_index == "bai" and needs_csi:
            sys.stderr.write(
                f"ERROR ~ [align] bam_index=bai requested but the longest contig is "
                f"{max_contig:,} bp, above the BAI limit of "
                f"{policy['bai_max_contig']:,}.\n        fix: --bam_index csi\n"
            )
            return 2
        bam_index, bam_reason, bam_source = args.bam_index, "explicitly requested", "user"
    elif needs_csi:
        bam_index = "csi"
        bam_reason = (f"longest contig {max_contig:,} bp exceeds the BAI ceiling "
                      f"{policy['bai_max_contig']:,} bp")
        bam_source = "policy"
    else:
        bam_index = "both"
        bam_reason = "all contigs are BAI-addressable; emitting both for tool compatibility"
        bam_source = "policy"
    decisions["bam_index"] = decision(bam_index, bam_reason, source=bam_source)

    # --------------------------------------------------------- evidence mode
    has_rna, has_protein, has_isoseq = flag(args.has_rna), flag(args.has_protein), flag(args.has_isoseq)
    if has_rna and has_protein and has_isoseq:
        evidence = "dual"
    elif has_rna and has_protein:
        evidence = "etp"
    elif has_rna:
        evidence = "et"
    elif has_protein:
        evidence = "ep"
    elif has_isoseq:
        evidence = "isoseq"
    else:
        evidence = "es"
        notes.append(
            "No RNA-seq and no protein evidence. BRAKER would fall back to "
            "GeneMark-ES self-training; expect substantially lower accuracy. "
            "Set --evidence_mode explicitly to silence this."
        )
    decisions["evidence_mode"] = decision(
        evidence,
        f"rna={has_rna} protein={has_protein} isoseq={has_isoseq} "
        f"({args.n_samples} sample(s))",
    )

    # -------------------------------------------------------------- predictors
    blocked_tracks = list(policy["predictor_blocked"].get(size_class, []))
    for tool, clades in CLADE_UNSUPPORTED.items():
        key = clade if clade in clades else (taxon if taxon in clades else None)
        if key:
            blocked_tracks.append(tool)
            notes.append(f"{tool} unavailable: {clades[key]}")

    primary = policy["predictor_primary"].get(size_class, "braker")
    secondary = policy["predictor_secondary"].get(size_class)
    if not flag(args.gpu) and primary == "helixer":
        primary = "braker_sharded" if size_class in ("large", "huge") else "braker"
        secondary = None
        notes.append(
            "Helixer would be the primary track at this genome size but no GPU is "
            "available; falling back to BRAKER. This is materially slower and the "
            "run may not finish -- consider a GPU host."
        )

    predictor_reason = f"size_class={size_class}"
    if size_class == "huge":
        predictor_reason += (
            "; BRAKER has no internal checkpointing and is not guaranteed to finish "
            "at this scale, so a GPU track that completes in days is run first and "
            "BRAKER follows in sharded, independently-resumable form"
        )
    decisions["predictors"] = decision(primary, predictor_reason)
    decisions["predictors"]["secondary"] = secondary
    decisions["predictors"]["blocked"] = sorted(set(blocked_tracks))

    braker_sharded = (total_bp > policy["braker_shard_threshold"]) or \
                     (secondary == "braker_sharded")
    decisions["braker"] = {
        "mode": "sharded" if braker_sharded else "native",
        "reason": (f"genome {total_bp/1e9:.2f} Gb exceeds the shard threshold "
                   f"{policy['braker_shard_threshold']/1e9:.0f} Gb; sharding turns one "
                   f"un-resumable multi-week process into independently resumable tasks"
                   if braker_sharded else "genome fits a single BRAKER run"),
        "fungus_model": taxon == "fungi",
    }
    if braker_sharded:
        decisions["braker"]["shard_mb"] = policy["braker_shard_mb"].get(size_class, 20)
        decisions["braker"]["shard_overlap_kb"] = policy["braker_shard_overlap_kb"].get(size_class, 200)

    # -------------------------------------------------------------- consensus
    tracks = [t for t in [primary, secondary] if t]
    if len(tracks) < 2:
        consensus, consensus_reason = "none", "only one predictor track"
    elif all(t.startswith("braker") or t == "galba" for t in tracks):
        consensus, consensus_reason = "tsebra", "all tracks are BRAKER-family"
    else:
        consensus, consensus_reason = "mikado", "tracks include a non-BRAKER isoform source"
    decisions["consensus"] = decision(consensus, consensus_reason)

    # ---------------------------------------------------------- disabled bits
    disabled = sorted({
        name for name, limit in policy["disable_above"].items() if total_bp > limit
    })
    if disabled:
        decisions["disabled_stages"] = decision(
            disabled, f"genome exceeds the size limit for these stages at {total_bp/1e9:.2f} Gb")

    # ---------------------------------------------------------------- budgets
    budget = {}
    for track in tracks + (["repeats"] if masker_choice != "none" else []):
        limit = args.time_budget or policy["track_budget"].get(track)
        if limit:
            budget[track] = {
                "max_time": limit,
                "on_exceed": "continue_without",
                "fallback": (f"drop {track}; finish with the remaining tracks and "
                             f"record the omission in the report"),
            }

    strategy = {
        "genome": {
            "id": args.genome_id,
            "species": args.species,
            "bp": total_bp,
            "contigs": stats["contigs"],
            "max_contig": max_contig,
            "n_fraction": stats.get("n_fraction"),
            "softmask_fraction": softmask_frac,
            "softmasked": softmask_frac >= 0.02,
            "size_class": size_class,
            "taxon": taxon,
            "clade": clade,
        },
        "reference_sets": {
            "busco_lineage": policy["busco_lineage"].get(taxon),
            "odb_partition": policy["odb_partition"].get(taxon),
        },
        "decisions": decisions,
        "budget": budget,
        "notes": notes,
        "priority": args.priority,
    }

    with open(args.out, "w") as fh:
        fh.write("# eukannot execution strategy -- generated, do not edit.\n")
        fh.write("# Every decision records why it was made and what it rejected.\n")
        fh.write(yaml_dump(strategy))
    with open(args.out_json, "w") as fh:
        json.dump(strategy, fh, indent=2)

    sys.stderr.write(
        f"INFO  ~ [strategy] {args.genome_id}: size_class={size_class} taxon={taxon} "
        f"masker={masker_choice} aligner={aligner} bam_index={bam_index} "
        f"evidence={evidence} primary={primary}"
        + (f" secondary={secondary}" if secondary else "") + "\n"
    )
    for note in notes:
        sys.stderr.write(f"WARN  ~ [strategy] {note}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
