#!/usr/bin/env python3
"""Unit tests for bin/decide_strategy.py -- the policy engine.

These are the tests that matter most: the engine is what makes the pipeline
decide rather than be told, so a silent regression here would produce a run that
looks fine and is wrong. Run directly:

    python3 tests/test_decide_strategy.py
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(os.path.dirname(HERE), "bin", "decide_strategy.py")

PASS, FAIL = 0, 0

# The two real reference genomes this pipeline was built against.
PINE = {  # Pinus densiflora: gymnosperm, 21.74 Gb, 2.008 Gb longest contig, unmasked
    "total_bp": 21737972491, "contigs": 2006, "max_contig_bp": 2007914973,
    "n_fraction": 0.0, "softmask_fraction": 0.0,
    "composition_scanned_bp": 200000000, "composition_is_sampled": True,
}
MATSUTAKE = {  # Tricholoma matsutake: fungus, 161 Mb, softmasked (71.3% per RepeatMasker)
    "total_bp": 161040721, "contigs": 13, "max_contig_bp": 19249005,
    "n_fraction": 2e-06, "softmask_fraction": 0.713173,
    "composition_scanned_bp": 161040721, "composition_is_sampled": False,
}


def check(name, condition, detail=""):
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  ok   {name}")
    else:
        FAIL += 1
        print(f"  FAIL {name}  {detail}")


def run(stats: dict, *args):
    tmp = tempfile.mkdtemp()
    stats_path = os.path.join(tmp, "stats.json")
    with open(stats_path, "w") as fh:
        json.dump(stats, fh)
    out_json = os.path.join(tmp, "s.json")
    proc = subprocess.run(
        [sys.executable, SCRIPT, "--genome-stats", stats_path,
         "--out", os.path.join(tmp, "s.yml"), "--out-json", out_json, *args],
        capture_output=True, text=True,
    )
    plan = None
    if os.path.exists(out_json) and os.path.getsize(out_json):
        with open(out_json) as fh:
            plan = json.load(fh)
    return proc, plan


def test_pine():
    print("Pinus densiflora (huge, gymnosperm, unmasked, GPU)")
    proc, plan = run(PINE, "--genome-id", "Pinde", "--species", "Pinus_densiflora",
                     "--taxon", "plant", "--has-rna", "true", "--gpu", "true")
    check("exit 0", proc.returncode == 0, proc.stderr)
    d = plan["decisions"]
    check("size_class=huge", d["size_class"]["choice"] == "huge")
    check("clade inferred as gymnosperm", plan["genome"]["clade"] == "gymnosperm")
    check("masker=red", d["masker"]["choice"] == "red")
    check("earlgrey/edta blocked", set(d["masker"]["blocked_at_this_size"]) == {"earlgrey", "edta"})
    check("masker reason cites the measured failure",
          "52 h" in d["masker"]["reason"], d["masker"]["reason"])
    check("aligner=hisat2", d["aligner"]["choice"] == "hisat2")
    check("large index", d["aligner"]["large_index"] is True)
    check("mmap index", d["aligner"]["mmap_index"] is True)
    check("bam_index=csi", d["bam_index"]["choice"] == "csi")
    check("csi reason shows the arithmetic",
          "2,007,914,973" in d["bam_index"]["reason"], d["bam_index"]["reason"])
    check("primary=helixer", d["predictors"]["choice"] == "helixer")
    check("secondary=braker_sharded", d["predictors"]["secondary"] == "braker_sharded")
    check("tiberius blocked", "tiberius" in d["predictors"]["blocked"])
    check("egapx blocked", "egapx" in d["predictors"]["blocked"])
    check("funannotate2 blocked", "funannotate2" in d["predictors"]["blocked"])
    check("braker sharded", d["braker"]["mode"] == "sharded")
    check("shard overlap 500 kb for conifer introns",
          d["braker"]["shard_overlap_kb"] == 500)
    check("busco lineage embryophyta", plan["reference_sets"]["busco_lineage"] == "embryophyta_odb12")
    check("odb partition Viridiplantae", plan["reference_sets"]["odb_partition"] == "Viridiplantae")
    check("every track has a budget with a fallback",
          all("fallback" in v for v in plan["budget"].values()), plan["budget"])


def test_matsutake():
    print("Tricholoma matsutake (small, fungus, premasked)")
    proc, plan = run(MATSUTAKE, "--genome-id", "Tmat", "--species", "Tricholoma_matsutake",
                     "--taxon", "fungi", "--premasked", "true",
                     "--has-rna", "true", "--has-protein", "true")
    check("exit 0", proc.returncode == 0, proc.stderr)
    d = plan["decisions"]
    check("size_class=small", d["size_class"]["choice"] == "small")
    check("masker=none (premasked accepted)", d["masker"]["choice"] == "none")
    check("bam_index=both", d["bam_index"]["choice"] == "both")
    check("evidence=etp", d["evidence_mode"]["choice"] == "etp")
    check("primary=braker", d["predictors"]["choice"] == "braker")
    check("braker native", d["braker"]["mode"] == "native")
    check("fungus model on", d["braker"]["fungus_model"] is True)
    check("consensus=none with one track", d["consensus"]["choice"] == "none")
    check("busco lineage fungi", plan["reference_sets"]["busco_lineage"] == "fungi_odb12")


def test_refusals():
    print("impossible requests are refused with the arithmetic")
    cases = [
        (["--premasked", "true"], "0.00% lowercase", "premasked on an unmasked genome"),
        (["--masker", "edta"], "cannot complete on a huge genome", "edta at 21.7 Gb"),
        (["--aligner", "star"], "not buildable at this scale", "star at 21.7 Gb"),
        (["--bam-index", "bai"], "above the BAI limit", "bai with a 2 Gb contig"),
    ]
    for args, needle, label in cases:
        proc, _ = run(PINE, "--genome-id", "Pinde", "--species", "Pinus_densiflora",
                      "--taxon", "plant", "--has-rna", "true", "--gpu", "true", *args)
        check(f"{label}: exit 2", proc.returncode == 2, proc.returncode)
        check(f"{label}: explains why", needle in proc.stderr, proc.stderr[:200])


def test_no_gpu_fallback():
    print("no GPU: Helixer is unavailable, fall back and say so")
    proc, plan = run(PINE, "--genome-id", "Pinde", "--species", "Pinus_densiflora",
                     "--taxon", "plant", "--has-rna", "true", "--gpu", "false")
    check("exit 0", proc.returncode == 0, proc.stderr)
    check("primary is no longer helixer",
          plan["decisions"]["predictors"]["choice"] != "helixer")
    check("the fallback is explained",
          any("no GPU is available" in n for n in plan["notes"]), plan["notes"])


def test_evidence_modes():
    print("evidence mode derivation")
    cases = [
        (("true", "true", "true"), "dual"),
        (("true", "true", "false"), "etp"),
        (("true", "false", "false"), "et"),
        (("false", "true", "false"), "ep"),
        (("false", "false", "true"), "isoseq"),
        (("false", "false", "false"), "es"),
    ]
    for (rna, prot, iso), expected in cases:
        _, plan = run(MATSUTAKE, "--taxon", "fungi", "--premasked", "true",
                      "--has-rna", rna, "--has-protein", prot, "--has-isoseq", iso)
        got = plan["decisions"]["evidence_mode"]["choice"]
        check(f"rna={rna[0]} prot={prot[0]} iso={iso[0]} -> {expected}", got == expected, got)
    _, plan = run(MATSUTAKE, "--taxon", "fungi", "--premasked", "true")
    check("no evidence warns about GeneMark-ES accuracy",
          any("lower accuracy" in n for n in plan["notes"]), plan["notes"])


def test_user_override_is_recorded():
    print("explicit choices are honoured and marked as user-sourced")
    _, plan = run(MATSUTAKE, "--taxon", "fungi", "--premasked", "true",
                  "--has-rna", "true", "--size-class", "medium")
    check("size_class overridden", plan["decisions"]["size_class"]["choice"] == "medium")
    check("source recorded as user", plan["decisions"]["size_class"]["source"] == "user")


def main() -> int:
    for test in [test_pine, test_matsutake, test_refusals, test_no_gpu_fallback,
                 test_evidence_modes, test_user_override_is_recorded]:
        test()
    print(f"\n{PASS} passed, {FAIL} failed")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
