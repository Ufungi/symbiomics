"""Meaningful contracts for the independent stdlib count/QC import layer."""
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "post_analysis" / "post"))
from importers import import_inputs, read_table, write_tsv
from qc import emit_qc


class ImportContracts(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.base = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def put(self, name, text):
        path = self.base / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
        return name

    def inp(self, path="counts.tsv", **values):
        return dict(input_id="in1", source_project="project1", assay_id="a1", format="htseq_matrix",
                    counts_path=path, stats_path="", stats_format="none", gene_universe="full",
                    reference_id="ref1", annotation_id="ann1", count_unit="fragment",
                    assignment_mode="exclusive", provenance_level="declared", **values)

    def sample(self, sample_id="s1", **values):
        row = dict(sample_id=sample_id, input_id="in1", source_sample_id=sample_id,
                   biological_unit_id="bio_" + sample_id, library_id="lib_" + sample_id,
                   study_id="study1", sample_sets="main", include_analysis="true", treatment="A")
        row.update(values)
        return row

    def run_import(self, inputs, samples):
        return import_inputs(inputs, samples, self.base, self.base / "out")

    def test_reorders_genes_and_preserves_excluded_zero_qc(self):
        self.put("counts.tsv", "gene_id\ts1\ts2\ng2\t3\t0\ng1\t5\t0\n")
        result = self.run_import([self.inp()], [self.sample("s2", include_analysis="false", exclusion_reason="zero"), self.sample()])
        self.assertEqual({m["sample_id"]: m["gene_count_sum"] for m in result["metrics"]}, {"s1": 8, "s2": 0})
        header, rows = read_table(result["assays"]["a1"]["counts_path"])
        self.assertEqual(header, ["gene_id", "s2", "s1"])
        self.assertEqual(rows[1]["s1"], "5")

    def test_missing_column_mapping_and_duplicate_genes_fail(self):
        self.put("counts.tsv", "gene_id\ts1\ts2\ng1\t1\t2\n")
        with self.assertRaisesRegex(ValueError, "1:1"):
            self.run_import([self.inp()], [self.sample()])
        self.put("counts.tsv", "gene_id\ts1\ng1\t1\ng1\t2\n")
        with self.assertRaisesRegex(ValueError, "duplicate gene"):
            self.run_import([self.inp()], [self.sample()])

    def test_invalid_and_fractional_integer_inputs_fail(self):
        for value in ["NaN", "-1", "1.1", "inf"]:
            self.put("counts.tsv", "gene_id\ts1\ng1\t" + value + "\n")
            with self.assertRaises(ValueError):
                self.run_import([self.inp()], [self.sample()])

    def test_htseq_raw_rate_has_exact_denominator(self):
        self.put("raw.tsv", "g1\t11\ng2\t0\n__no_feature\t4752\n__ambiguous\t0\n__alignment_not_unique\t86124\n__not_aligned\t35041973\n__too_low_aQual\t108062\n")
        inp = self.inp("raw.tsv"); inp.update(format="htseq_raw", stats_format="htseq_inline")
        result = self.run_import([inp], [self.sample()])
        metric = result["metrics"][0]
        self.assertEqual(metric["assignment_denominator"], 35240922)
        self.assertGreater(metric["assignment_rate_pct"], 0)
        self.assertLess(metric["assignment_rate_pct"], .0001)
        self.assertFalse(metric["zero_assigned"])
        paths = emit_qc(result, self.base / "out")
        _, rows = read_table(paths["sample_metrics"])
        self.assertEqual(rows[0]["assignment_denominator"], "35240922")

    def test_missing_counter_is_not_zero(self):
        self.put("raw.tsv", "g1\t8\n__no_feature\t2\n")
        inp = self.inp("raw.tsv"); inp.update(format="htseq_raw", stats_format="htseq_inline")
        metric = self.run_import([inp], [self.sample()])["metrics"][0]
        self.assertIsNone(metric["assignment_rate_pct"])
        self.assertEqual(metric["assigned_count"], 8)

    def test_native_counter_rate_retains_unknown_physical_unit(self):
        self.put("raw.tsv", "g1\t8\n__no_feature\t2\n__ambiguous\t0\n__alignment_not_unique\t0\n__not_aligned\t0\n__too_low_aQual\t0\n")
        inp = self.inp("raw.tsv"); inp.update(format="htseq_raw", count_unit="unknown")
        metric = self.run_import([inp], [self.sample()])["metrics"][0]
        self.assertEqual(metric["assignment_counter_rate_pct"], 80)
        self.assertEqual(metric["assignment_counter_status"], "ok")
        self.assertIsNone(metric["assignment_rate_pct"])
        self.assertIn("unverified", metric["assignment_counter_reason"])
        self.put("raw.tsv", "g1\t8\n__no_feature\t2\n")
        metric = self.run_import([inp], [self.sample()])["metrics"][0]
        self.assertIsNone(metric["assignment_counter_rate_pct"])
        self.assertIn("Incomplete HTSeq", metric["assignment_reason"])
        self.assertIn("Unknown physical", metric["assignment_reason"])

    def test_prefiltered_htseq_uses_reported_assigned_not_supplied_sum(self):
        self.put("counts.tsv", "gene_id\ts1\ng1\t4\n")
        self.put("stats.tsv", "sample\tassigned\tno_feature\tambiguous\talignment_not_unique\tnot_aligned\ttoo_low_aQual\ns1\t8\t2\t0\t0\t0\t0\n")
        inp = self.inp(); inp.update(gene_universe="prefiltered",stats_format="htseq_legacy",stats_path="stats.tsv")
        metric = self.run_import([inp], [self.sample()])["metrics"][0]
        self.assertEqual(metric["gene_count_sum"], 4)
        self.assertEqual(metric["assigned_count"], 8)
        self.assertEqual(metric["assignment_counter_rate_pct"], 80)
        self.assertIsNone(metric["genes_detected_pct"])
        self.put("stats.tsv", "sample\tno_feature\tambiguous\talignment_not_unique\tnot_aligned\ttoo_low_aQual\ns1\t2\t0\t0\t0\t0\n")
        metric = self.run_import([inp], [self.sample()])["metrics"][0]
        self.assertIsNone(metric["assignment_counter_rate_pct"])

    def test_prefiltered_featurecounts_uses_full_summary(self):
        self.put("counts.tsv", "gene_id\ts1\ng1\t4\n")
        self.put("summary.tsv", "Status\ts1\nAssigned\t8\nUnassigned_Unmapped\t2\n")
        inp = self.inp(); inp.update(format="featurecounts_matrix",gene_universe="prefiltered",stats_format="featurecounts_summary",stats_path="summary.tsv")
        metric = self.run_import([inp], [self.sample()])["metrics"][0]
        self.assertEqual(metric["assignment_counter_rate_pct"], 80)

    def test_nonexclusive_count_does_not_create_rate(self):
        self.put("raw.tsv", "g1\t2\ng2\t2\n__no_feature\t0\n__ambiguous\t2\n__alignment_not_unique\t0\n__not_aligned\t0\n__too_low_aQual\t0\n")
        inp = self.inp("raw.tsv"); inp.update(format="htseq_raw", assignment_mode="multiassign")
        metric = self.run_import([inp], [self.sample()])["metrics"][0]
        self.assertEqual(metric["gene_count_sum"], 4)
        self.assertIsNone(metric["assignment_rate_pct"])
        self.assertFalse(metric["composition_valid"])

    def test_featurecounts_count_and_stats_have_distinct_ids(self):
        self.put("counts.tsv", "gene_id\ts1.count\ng1\t8\n")
        self.put("summary.tsv", "Status\t/path/s1.count.bam\nAssigned\t8\nUnassigned_Unmapped\t2\n")
        inp = self.inp(); inp.update(format="featurecounts_matrix", stats_path="summary.tsv", stats_format="featurecounts_summary")
        result = self.run_import([inp], [self.sample(source_sample_id="s1.count", stats_sample_id="/path/s1.count.bam")])
        self.assertEqual(result["metrics"][0]["assignment_rate_pct"], 80)
        with self.assertRaisesRegex(ValueError, "sample field"):
            self.run_import([inp], [self.sample(source_sample_id="s1.count")])

    def test_raw_featurecounts_parser_preserves_exact_column_identity(self):
        self.put("raw.txt", "# Program:featureCounts\nGeneid\tChr\tStart\tEnd\tStrand\tLength\t/path/s1.count.bam\ng1\t1\t1\t20\t+\t20\t8\n")
        inp = self.inp("raw.txt"); inp["format"] = "featurecounts_raw"
        result = self.run_import([inp], [self.sample(source_sample_id="/path/s1.count.bam")])
        self.assertEqual(result["metrics"][0]["gene_count_sum"], 8)
        self.assertIsNone(result["metrics"][0]["assignment_rate_pct"])

    def test_csv_matrix_accepts_quoted_exact_gene_ids(self):
        self.put("counts.csv", 'gene_id,s1\n"gene,one",3\n')
        result = self.run_import([self.inp("counts.csv")], [self.sample()])
        _, rows = read_table(result["assays"]["a1"]["counts_path"])
        self.assertEqual(rows[0]["gene_id"], "gene,one")

    def test_legacy_source_id_and_large_integer_precision(self):
        self.put("counts.tsv", "gene_id\ts1\ng1\t9007199254740993\n")
        self.put("stats.tsv", "sample_id\tsample_id_original\tassigned\tno_feature\tambiguous\talignment_not_unique\tnot_aligned\ttoo_low_aQual\nlabel\toriginal\t9007199254740993\t1\t0\t0\t0\t0\n")
        inp = self.inp(); inp.update(stats_path="stats.tsv", stats_format="htseq_legacy")
        result = self.run_import([inp], [self.sample(stats_sample_id="original")])
        metric = result["metrics"][0]
        self.assertEqual(metric["gene_count_sum"], 9007199254740993)
        self.assertEqual(metric["assignment_denominator"], 9007199254740994)
        self.assertEqual(result["assays"]["a1"]["counts_path"], "")
        self.assertIn("double-precision", result["errors"][0])
        paths = emit_qc(result, self.base / "out")
        self.assertIn("9007199254740993", Path(paths["sample_metrics"]).read_text())

    def test_legacy_display_label_swap_uses_original_id_only(self):
        self.put("counts.tsv", "gene_id\ts1\ts2\ng1\t1\t2\n")
        self.put("stats.tsv", "sample_id\tsample_id_original\tassigned\tno_feature\tambiguous\talignment_not_unique\tnot_aligned\ttoo_low_aQual\ns2\ts1\t1\t0\t0\t0\t0\t0\ns1\ts2\t2\t0\t0\t0\t0\t0\n")
        inp = self.inp(); inp.update(stats_path="stats.tsv", stats_format="htseq_legacy")
        result = self.run_import([inp], [self.sample(), self.sample("s2")])
        self.assertEqual([m["assigned_count"] for m in result["metrics"]], [1, 2])

    def test_reference_conflict_preserves_raw_qc(self):
        self.put("a.tsv", "gene_id\ts1\ng1\t1\n")
        self.put("b.tsv", "gene_id\ts2\ng1\t2\n")
        a = self.inp("a.tsv"); b = self.inp("b.tsv"); b.update(input_id="in2", reference_id="different")
        result = self.run_import([a, b], [self.sample(), self.sample("s2", input_id="in2")])
        self.assertEqual(len(result["metrics"]), 2)
        self.assertTrue(result["assays"]["a1"]["errors"])
        self.assertEqual(result["assays"]["a1"]["counts_path"], "")

    def test_technical_merge_recalculates_rate_and_detection(self):
        self.put("counts.tsv", "gene_id\ts1\ts2\ng1\t8\t0\ng2\t0\t2\n")
        self.put("stats.tsv", "sample\tno_feature\tambiguous\talignment_not_unique\tnot_aligned\ttoo_low_aQual\ns1\t2\t0\t0\t0\t0\ns2\t98\t0\t0\t0\t0\n")
        inp = self.inp(); inp.update(stats_path="stats.tsv", stats_format="htseq_legacy")
        members = [self.sample(sid, biological_unit_id="bio", library_id="lib", technical_merge_id="merged",
                               technical_partition_id=sid, technical_disjoint="true") for sid in ["s1", "s2"]]
        result = self.run_import([inp], members)
        metric = next(m for m in result["metrics"] if m["sample_id"] == "merged")
        self.assertAlmostEqual(metric["assignment_rate_pct"], 100 * 10 / 110)
        self.assertEqual(metric["genes_detected"], 2)
        self.assertEqual(len(result["metrics"]), 3)
        self.assertEqual(result["assays"]["a1"]["samples"][0]["sample_id"], "merged")
        inp["count_unit"] = "unknown"
        unknown = self.run_import([inp], members)
        aggregate = next(m for m in unknown["metrics"] if m["sample_id"] == "merged")
        self.assertIsNone(aggregate["assignment_rate_pct"])
        self.assertIsNone(aggregate["assignment_counter_rate_pct"])
        members[1]["treatment"] = "B"
        with self.assertRaisesRegex(ValueError, "disagree on treatment"):
            self.run_import([inp], members)

    def test_explicit_metadata_conflict_rejects_false_exclusive_claim(self):
        self.put("counts.tsv", "gene_id\ts1\ng1\t1\n")
        self.put("run.json", json.dumps({"nonunique": "all"}))
        inp = self.inp(); inp["run_metadata_path"] = "run.json"
        with self.assertRaisesRegex(ValueError, "exclusive assignment conflicts"):
            self.run_import([inp], [self.sample()])

    def test_prefiltered_matrix_does_not_claim_full_detection(self):
        self.put("counts.tsv", "gene_id\ts1\ng1\t1\n")
        inp = self.inp(); inp["gene_universe"] = "prefiltered"
        metric = self.run_import([inp], [self.sample()])["metrics"][0]
        self.assertIsNone(metric["genes_detected_pct"])
        self.assertEqual(metric["supplied_gene_detection_pct"], 100)

    def test_salmon_mapping_and_required_gene_map(self):
        self.put("quant.sf", "Name\tLength\tEffectiveLength\tTPM\tNumReads\nt1\t100\t80\t500000\t2.25\nt2\t100\t80\t500000\t1.75\n")
        self.put("map.tsv", "transcript_id\tgene_id\nt1\tg1\nt2\tg1\n")
        self.put("meta.json", json.dumps({"num_processed": 100, "num_mapped": 20}))
        inp = self.inp("quant.sf"); inp.update(format="salmon_quant", stats_path="meta.json", stats_format="salmon_meta", gene_map_path="map.tsv", assignment_mode="fractional", count_unit="estimated_fragment")
        result = self.run_import([inp], [self.sample()])
        self.assertEqual(result["metrics"][0]["mapping_rate_pct"], 20)
        self.assertEqual(result["metrics"][0]["gene_count_sum"], 4)
        self.assertIsNone(result["metrics"][0]["assignment_rate_pct"])
        self.put("map.tsv", "transcript_id\tgene_id\nt1\tg1\n")
        with self.assertRaisesRegex(ValueError, "unmapped transcript"):
            self.run_import([inp], [self.sample()])


if __name__ == "__main__":
    unittest.main()
