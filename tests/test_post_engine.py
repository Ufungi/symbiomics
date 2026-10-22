"""Opt-in acceptance checks for the real R engine.

Run with SYMBIOMICS_TEST_R=1 and SYMBIOMICS_RSCRIPT=/path/to/Rscript.
These checks execute edgeR/tximport, not substitutes for their calculations.
"""
import csv
import json
import math
import os
from pathlib import Path
import random
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "post_analysis/post/R/engine.R"


def write_tsv(path, rows, fields=None):
    fields = fields or list(rows[0])
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with open(path, encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@unittest.skipUnless(os.environ.get("SYMBIOMICS_TEST_R") == "1", "opt in with SYMBIOMICS_TEST_R=1; requires the post R environment")
class PostEngineAcceptance(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rscript = os.environ.get("SYMBIOMICS_RSCRIPT") or shutil.which("Rscript")
        if not cls.rscript:
            raise unittest.SkipTest("Set SYMBIOMICS_RSCRIPT to the post environment's Rscript")
        check = subprocess.run([cls.rscript, "-e", "stopifnot(all(vapply(c('jsonlite','edgeR','limma','statmod','tximport','fgsea','sva'),requireNamespace,logical(1),quietly=TRUE)))"], capture_output=True, text=True)
        if check.returncode:
            raise unittest.SkipTest("R environment lacks required analysis packages: " + check.stderr[-600:])

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="symbiomics-engine-test-")
        self.addCleanup(self.tmp.cleanup)
        self.work = Path(self.tmp.name)
        self.samples = []
        for timepoint, time_h in [("early", "0"), ("late", "24")]:
            for treatment in ["control", "treated"]:
                for replicate in range(4):
                    sid = f"s{len(self.samples) + 1:02}"
                    self.samples.append(dict(sample_id=sid, biological_unit_id="bio_" + sid, subject_id="", treatment=treatment, timepoint=timepoint, time_h=time_h, batch=f"b{replicate % 2 + 1}"))
        for subject in range(3):
            for treatment in ["control", "treated"]:
                self.samples.append(dict(sample_id=f"p{subject}_{treatment}", biological_unit_id=f"subject_{subject}", subject_id=f"subject_{subject}", treatment=treatment, timepoint="", time_h="", batch="b1"))
        rng = random.Random(321)
        self.count_rows = []
        for gene in range(500):
            row = {"gene_id": f"g{gene}"}
            for sample in self.samples:
                value = max(1, round(rng.gauss(100, 24)))
                if gene < 25 and sample["treatment"] == "treated":
                    value *= 8
                if 25 <= gene < 50 and sample["timepoint"] == "late" and sample["treatment"] == "treated":
                    value *= 6
                if 50 <= gene < 75 and sample["timepoint"] == "late":
                    value *= 4
                row[sample["sample_id"]] = value
            self.count_rows.append(row)
        self.counts = self.work / "counts.tsv"
        write_tsv(self.counts, self.count_rows)
        self.basic = dict(analysis_id="case", assay_id="demo", analysis_mode="inferential", design_type="group", factor_cols="treatment", numeric_cols="", nuisance_factor_cols="batch", nuisance_numeric_cols="", interaction_terms="", subject_col="", reference_levels="treatment=control;batch=b1", numeric_units="", normalization="TMM", fdr_cut="0.05", lfc_cut="1", effect_test="qlf", enrichment="none", adjusted_view="none", seed="123")
        self.simple = self.cells("treated_vs_control", [("treatment=treated", 1), ("treatment=control", -1)])

    @staticmethod
    def cells(contrast_id, values, test_id=""):
        return [dict(analysis_id="case", contrast_id=contrast_id, basis="cell", term=term, weight=str(weight), test_id=test_id) for term, weight in values]

    def run_job(self, name, analysis=None, samples=None, contrasts=None, expected="ok", mode="run", **extra):
        directory = self.work / name
        directory.mkdir()
        selected = samples if samples is not None else self.samples[:8]
        write_tsv(directory / "samples.tsv", selected)
        analysis = dict(self.basic, **(analysis or {}))
        analysis["analysis_id"] = name
        contrasts = [dict(row, analysis_id=name) for row in (self.simple if contrasts is None else contrasts)]
        job = dict(mode=mode, counts_path=str(self.counts), samples_path=str(directory / "samples.tsv"), analysis=analysis, contrasts=contrasts, annotations=[], out_dir=str(directory / "out"))
        job.update(extra)
        job_path = directory / "job.json"
        job_path.write_text(json.dumps(job), encoding="utf-8")
        process = subprocess.run([self.rscript, str(ENGINE), str(job_path)], capture_output=True, text=True, timeout=120)
        (directory / "stdout.log").write_text(process.stdout)
        (directory / "stderr.log").write_text(process.stderr)
        status = json.loads((directory / "out/status.json").read_text())
        self.assertEqual(status["status"], expected, process.stderr + json.dumps(status, indent=2))
        self.assertEqual(process.returncode == 0, expected == "ok")
        return directory / "out", status

    def test_group_signal_and_input_permutation(self):
        out, _ = self.run_job("group")
        rows = {row["gene_id"]: row for row in read_tsv(out / "contrasts/treated_vs_control/DEG.tsv")}
        self.assertTrue(all(float(rows[f"g{i}"]["logFC"]) > 2 for i in range(25)))
        self.assertGreater(sum(float(rows[f"g{i}"]["FDR"]) < .05 for i in range(25)), 20)
        permuted, _ = self.run_job("permuted", samples=list(reversed(self.samples[:8])))
        other = {row["gene_id"]: row for row in read_tsv(permuted / "contrasts/treated_vs_control/DEG.tsv")}
        for gene, row in rows.items():
            self.assertAlmostEqual(float(row["logFC"]), float(other[gene]["logFC"]), places=8)
            self.assertAlmostEqual(float(row["FDR"]), float(other[gene]["FDR"]), places=8)

    def test_factorial_and_omnibus(self):
        analysis = dict(design_type="factorial", factor_cols="treatment|timepoint", interaction_terms="treatment:timepoint", reference_levels="treatment=control;timepoint=early;batch=b1")
        terms = self.cells("late_effect", [("treatment=treated;timepoint=late", 1), ("treatment=control;timepoint=late", -1)], "joint")
        terms += self.cells("interaction", [("treatment=treated;timepoint=late", 1), ("treatment=control;timepoint=late", -1), ("treatment=treated;timepoint=early", -1), ("treatment=control;timepoint=early", 1)], "joint")
        out, _ = self.run_job("factorial", analysis, self.samples[:16], terms)
        omnibus = read_tsv(out / "omnibus/joint/DEG.tsv")
        self.assertFalse(any(key.startswith("logFC") or key == "Regulation" for key in omnibus[0]))
        self.assertEqual(json.loads((out / "omnibus/joint/settings.json").read_text())["df"], 2)
        interaction = read_tsv(out / "contrasts/interaction/DEG.tsv")
        self.assertGreater(sum(float(row["logFC"]) > 1.5 for row in interaction if 25 <= int(row["gene_id"][1:]) < 50), 20)

    def test_paired_and_trend_coefficient_scaling(self):
        self.run_job("paired", dict(design_type="paired", nuisance_factor_cols="", subject_col="subject_id", reference_levels="treatment=control;subject_id=subject_0"), self.samples[16:])
        terms = [dict(analysis_id="case", contrast_id="per_hour", basis="coefficient", term="numeric:time_h", weight="1", test_id=""), dict(analysis_id="case", contrast_id="per_day", basis="coefficient", term="numeric:time_h", weight="24", test_id="")]
        out, _ = self.run_job("trend", dict(design_type="trend", numeric_cols="time_h", numeric_units="time_h=hour"), self.samples[:16], terms)
        hour = read_tsv(out / "contrasts/per_hour/DEG.tsv")
        day = read_tsv(out / "contrasts/per_day/DEG.tsv")
        for h, d in zip(hour, day):
            self.assertEqual(h["gene_id"], d["gene_id"])
            self.assertAlmostEqual(float(h["logFC"]) * 24, float(d["logFC"]), places=8)
            self.assertAlmostEqual(float(h["FDR"]), float(d["FDR"]), places=8)

    def test_guard_confounded_and_repeated_units(self):
        samples = [dict(row, batch=row["treatment"]) for row in self.samples[:8]]
        _, status = self.run_job("confounded", dict(reference_levels="treatment=control;batch=control"), samples, expected="failed", mode="validate")
        self.assertIn("Rank-deficient", status["error"])
        samples = [dict(row, biological_unit_id="same") for row in self.samples[:8]]
        _, status = self.run_job("pseudoreplication", samples=samples, expected="failed", mode="validate")
        self.assertIn("Repeated biological units", status["error"])

    def test_exploratory_without_experiment_assertion(self):
        samples = [{key: value for key, value in row.items() if key not in ("biological_unit_id", "subject_id")} for row in self.samples[:8]]
        out, status = self.run_job("exploratory", dict(analysis_mode="exploratory"), samples, [])
        self.assertEqual(status["stages"]["deg"]["status"], "skipped")
        self.assertTrue((out / "exploratory/pca_scores.tsv").exists())
        self.assertFalse(list(out.rglob("DEG.tsv")))

    def test_zero_retained_library_is_explicit(self):
        sparse = self.work / "sparse.tsv"
        names = [row["sample_id"] for row in self.samples[:3]]
        write_tsv(sparse, [dict(gene_id="common", **dict(zip(names, [100, 100, 0]))), dict(gene_id="private", **dict(zip(names, [0, 0, 1])))])
        out, status = self.run_job("sparse", dict(analysis_mode="exploratory"), self.samples[:3], [], expected="failed", counts_path=str(sparse))
        self.assertIn("zero retained counts", status["error"])
        rows = read_tsv(out / "normalized/sample_filter_status.tsv")
        self.assertEqual(rows[-1]["normalizable"], "FALSE")

    def test_treat_and_enrichment_and_partial_failure(self):
        self.run_job("treat", dict(effect_test="treat"))
        go = self.work / "go.tsv"
        pathway = self.work / "pathway.tsv"
        write_tsv(go, [dict(gene_id=f"g{i}", go_id=f"GO:{i // 25}") for i in range(500)])
        write_tsv(pathway, [dict(gene_id=f"g{i}", pathway_id=f"p{i // 25}", pathway_name=f"set{i // 25}") for i in range(500)])
        annotations = [dict(annotation_id="go", assay_id="demo", kind="gene_go", path=str(go)), dict(annotation_id="pathway", assay_id="demo", kind="gene_pathway", path=str(pathway))]
        analysis = dict(enrichment="go|pathway", go_method="fisher_ora", pathway_methods="camera|fgsea")
        out, _ = self.run_job("enrichment", analysis, annotations=annotations)
        dest = out / "contrasts/treated_vs_control"
        for name in ["go_go_ora.tsv", "pathway_camera.tsv", "pathway_fgsea.tsv"]:
            self.assertTrue((dest / name).exists())
            self.assertGreater(len(read_tsv(dest / name)), 0)
        annotations[0]["path"] = str(self.work / "absent.tsv")
        out, status = self.run_job("partial", analysis, annotations=annotations, expected="partial_failure")
        self.assertEqual(status["stages"]["contrast:treated_vs_control"]["status"], "ok")
        self.assertTrue((out / "contrasts/treated_vs_control/DEG.tsv").exists())

    def test_deferred_annotation_schema_failure_preserves_inference(self):
        out, status = self.run_job("annotation_schema_failure", dict(enrichment="go", go_method="fisher_ora"), annotations=[], annotation_error="annotations.tsv: unknown columns: wrong_column", expected="partial_failure")
        self.assertEqual(status["stages"]["contrast:treated_vs_control"]["status"], "ok")
        enrichment = status["stages"]["enrichment:treated_vs_control"]
        self.assertEqual(enrichment["status"], "failed")
        self.assertIn("unknown columns", enrichment["reason"])
        self.assertTrue((out / "contrasts/treated_vs_control/DEG.tsv").exists())
        self.assertTrue((out / "exploratory/pca_scores.tsv").exists())

    def test_combat_preserves_raw_statistical_input(self):
        out, status = self.run_job("combat", dict(adjusted_view="combat_seq", adjusted_batch_col="batch"))
        self.assertEqual(status["stages"]["adjusted_view"]["status"], "ok")
        self.assertTrue((out / "exploratory_adjusted/logCPM_visualization_only.tsv").exists())
        baseline, _ = self.run_job("baseline")
        left = read_tsv(out / "contrasts/treated_vs_control/DEG.tsv")
        right = read_tsv(baseline / "contrasts/treated_vs_control/DEG.tsv")
        for a, b in zip(left, right):
            self.assertEqual(a["gene_id"], b["gene_id"])
            self.assertEqual(a["logFC"], b["logFC"])
            self.assertEqual(a["PValue"], b["PValue"])

    def test_salmon_import_retains_length_offsets(self):
        txmap = self.work / "tx2gene.tsv"
        write_tsv(txmap, [dict(tx_id=f"t{i}", gene_id=f"g{i}") for i in range(500)])
        salmon = []
        for sample in self.samples[:8]:
            rows = []
            for i, gene in enumerate(self.count_rows):
                effective = (300 + 3 * i) * (1.2 if sample["treatment"] == "treated" else 1)
                rows.append(dict(Name=f"t{i}", Length=effective + 100, EffectiveLength=effective, TPM=gene[sample["sample_id"]] / effective, NumReads=gene[sample["sample_id"]]))
            total = sum(row["TPM"] for row in rows)
            for row in rows:
                row["TPM"] *= 1e6 / total
            quant = self.work / f"{sample['sample_id']}_quant.sf"
            write_tsv(quant, rows)
            salmon.append(dict(sample_id=sample["sample_id"], quant_path=str(quant), gene_map_path=str(txmap)))
        out, _ = self.run_job("salmon", salmon=salmon)
        settings = json.loads((out / "normalized/settings.json").read_text())
        self.assertTrue(settings["length_offset"])
        self.assertEqual(settings["source"], "salmon_tximport_length_offset")
        offsets = read_tsv(out / "normalized/log_offsets.tsv")
        self.assertEqual(len(offsets), 500)
        self.assertTrue(all(math.isfinite(float(value)) for row in offsets for key, value in row.items() if key != "gene_id"))

    def test_quoted_tsv_gene_and_metadata_roundtrip(self):
        self.count_rows[0]["gene_id"] = 'gene"one'
        write_tsv(self.counts, self.count_rows)
        samples = [dict(row, treatment='control "arm"' if row["treatment"] == "control" else 'treated "arm"') for row in self.samples[:8]]
        contrasts = self.cells("quoted", [('treatment=treated "arm"', 1), ('treatment=control "arm"', -1)])
        go = self.work / "quoted_go.tsv"
        write_tsv(go, [dict(gene_id=row["gene_id"], go_id=f"GO:{i // 25}") for i, row in enumerate(self.count_rows)])
        annotations = [dict(annotation_id="quoted_go", assay_id="demo", kind="gene_go", path=str(go))]
        out, _ = self.run_job("quoted", dict(reference_levels='treatment=control "arm";batch=b1', enrichment="go", go_method="fisher_ora"), samples, contrasts, annotations=annotations)
        genes = read_tsv(out / "contrasts/quoted/DEG.tsv")
        self.assertIn('gene"one', [row["gene_id"] for row in genes])
        levels = read_tsv(out / "design/factor_levels.tsv")
        self.assertIn('control "arm"', [row["label"] for row in levels])
        background = read_tsv(out / "contrasts/quoted/quoted_go_background.tsv")
        self.assertIn('gene"one', [row["gene_id"] for row in background])
        runtime = json.loads((out / "runtime.json").read_text())
        self.assertTrue(runtime["packages"]["edgeR"])

    def test_encoded_factor_levels_are_decoded_once(self):
        samples = [dict(row, treatment="control;100%" if row["treatment"] == "control" else "treated+value") for row in self.samples[:8]]
        contrasts = self.cells("encoded", [("treatment=treated%2Bvalue", 1), ("treatment=control%3B100%25", -1)])
        out, _ = self.run_job("encoded", dict(reference_levels="treatment=control%3B100%25;batch=b1"), samples, contrasts)
        levels = read_tsv(out / "design/factor_levels.tsv")
        self.assertIn("control;100%", [row["label"] for row in levels])
        self.assertIn("treated+value", [row["label"] for row in levels])


if __name__ == "__main__":
    unittest.main()
