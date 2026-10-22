"""Report contracts: availability, units, safe markup and metadata-driven projections."""
import csv
import json
import re
import sys
import tempfile
import unittest
from pathlib import Path
from xml.etree import ElementTree

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "post_analysis" / "post"))
from report import render_report


class ReportContracts(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="post-report-test-")
        self.run = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def table(self, relative, rows, fields=None):
        path = self.run / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        return path

    def read_table(self, path):
        with Path(path).open(encoding="utf-8", newline="") as stream:
            return list(csv.DictReader(stream, delimiter="\t"))

    def metric(self, sample_id, value, **extra):
        row = dict(sample_id=sample_id, assay_id="fungus", study_id="study", count_unit="fragment",
                   quant_tool="htseq", assignment_scope="counter_input", gene_universe="full", qc_level="library",
                   assignment_rate_pct=value, assigned_count=1 if value else value, gene_count_sum=1 if value else value,
                   counting_processed_count=1000000000, genes_detected=1 if value else value,
                   quantified_gene_count=1000, genes_detected_pct=.1 if value else value,
                   biological_unit_id="bio_"+sample_id, tissue="root")
        row.update(extra)
        return row

    def test_zero_missing_and_tiny_positive_are_distinct_and_escaped(self):
        metrics = [self.metric("positive<script>", 1e-7), self.metric("zero", 0), self.metric("missing", None)]
        metrics.append(self.metric("aggregate", 1e-7, qc_level="technical_aggregate"))
        path = render_report(self.run, metrics, [], [dict(stage="<script>", status="failed", reason='bad <counter> & "quote"', path="javascript:alert(1)")], {})
        document = path.read_text()
        self.assertNotIn("<script>", document)
        self.assertIn("&lt;script&gt;", document)
        self.assertNotIn('href="javascript:', document)
        self.assertNotRegex(document, r'(?:src|href)="https?://')
        self.assertIn("1.0000e-07", document)
        self.assertIn('<b>3</b><span>Input observations', document)
        groups = self.read_table(self.run/"00_quant_qc/group_summary.tsv")
        raw = [g for g in groups if g["metric"]=="assignment_rate_pct" and "QC=library" in g["denominator_scope"]][0]
        self.assertEqual((raw["n_total"], raw["n_valid"], raw["n_missing"], raw["n_zero"]), ("3","2","1","1"))
        aggregates = [g for g in groups if g["metric"]=="assignment_rate_pct" and "QC=technical_aggregate" in g["denominator_scope"]]
        self.assertEqual(aggregates[0]["n_total"], "1")
        sidecar = [r for p in (self.run/"00_quant_qc/plots").glob("*.tsv") for r in self.read_table(p) if "assignment_rate_pct" in r]
        self.assertTrue(any(r["assignment_rate_pct"]=="1e-07" for r in sidecar))
        self.assertTrue(any(r["assignment_rate_pct"]=="NA" for r in sidecar))
        for svg in (self.run/"00_quant_qc/plots").glob("*.svg"):
            ElementTree.parse(svg)

    def test_native_counter_fraction_keeps_unknown_physical_unit_explicit(self):
        row = self.metric("native", None, count_unit="unknown", assignment_mode="exclusive",
                          assignment_counter_rate_pct=1.23e-5, assignment_counter_status="ok",
                          assignment_counter_reason="Declared exclusive counters; physical unit unverified",
                          assigned_count=123, gene_count_sum=123, counting_processed_count=1000000000)
        doc = render_report(self.run, [row], [], [], {}, [dict(plot_id="native_fraction",scope="quant_qc",metric="assignment_counter_rate_pct",x="study_id",scale="log10")]).read_text()
        self.assertIn("Native counter assignment (%)", doc)
        self.assertIn("1.2300e-05", doc)
        self.assertIn("physical unit unverified", doc)
        self.assertIn("unknown physical units are not pooled", doc)
        summary = self.read_table(self.run/"00_quant_qc/group_summary.tsv")
        physical = [r for r in summary if r["metric"]=="assignment_rate_pct"][0]
        native = [r for r in summary if r["metric"]=="assignment_counter_rate_pct"][0]
        self.assertEqual(physical["n_missing"], "1")
        self.assertEqual(native["n_valid"], "1")
        self.assertEqual(native["count_unit"], "unknown")
        self.assertEqual(self.read_table(self.run/"plot_status.tsv")[0]["status"], "ok")
        # The first table is the pooled summary; all three unknown-unit totals are NA.
        pooled = doc.split("<table>",1)[1].split("</table>",1)[0]
        self.assertIn("<td>NA</td><td>NA</td><td>NA</td>", pooled)

    def test_mapping_units_and_incomplete_gene_universe(self):
        row = self.metric("salmon", None, quant_tool="salmon", count_unit="estimated_fragment",
                          mapping_unit="fragment", mapping_scope="salmon_input", mapping_rate_pct=4,
                          gene_universe="unknown", supplied_gene_detection_pct=3)
        page = render_report(self.run, [row], [], [], {}).read_text()
        self.assertIn("Detected out of 1,000 supplied genes (%)", page)
        summaries = self.read_table(self.run/"00_quant_qc/group_summary.tsv")
        mapping = [r for r in summaries if r["metric"]=="mapping_rate_pct"]
        self.assertEqual(mapping[0]["count_unit"], "fragment")
        self.assertIn("salmon_input", mapping[0]["denominator_scope"])
        assigned = [r for r in summaries if r["metric"]=="assignment_rate_pct"]
        self.assertEqual(assigned[0]["count_unit"], "estimated_fragment")

    def setup_coordinates(self):
        self.table("analyses/main/input/samples.tsv", [
            dict(sample_id="b", tissue="stem", condition="treated", biological_unit_id="b"),
            dict(sample_id="a", tissue="root", condition="control<script>", biological_unit_id="a"),
            dict(sample_id="c", tissue="stem", condition="control<script>", biological_unit_id="c"),
        ])
        self.table("analyses/main/exploratory/pca_scores.tsv", [
            dict(sample_id="c", PC1="1.25", PC2="-0.5", group="engine_default"),
            dict(sample_id="a", PC1="-2.5", PC2="0.125", group="engine_default"),
            dict(sample_id="b", PC1="0.75", PC2="0.375", group="engine_default"),
        ])
        self.table("analyses/main/exploratory/pca_variance.tsv", [
            dict(component=1, variance_fraction=.75), dict(component=2, variance_fraction=.25)])
        self.table("analyses/main/exploratory/mds_scores.tsv", [
            dict(sample_id="c", MDS1="-3", MDS2="2"), dict(sample_id="a", MDS1="4", MDS2="-5"), dict(sample_id="b", MDS1="1", MDS2="3")])
        pdf=self.run/"analyses/main/exploratory/PCA.pdf"
        pdf.write_bytes(b"%PDF-1.4")
        return [dict(analysis_id="main", analysis_mode="exploratory")]

    def test_pca_joins_by_id_preserves_scores_and_facets_colors(self):
        analyses=self.setup_coordinates()
        request=dict(plot_id="pca_by_tissue", scope="analysis", analysis_id="main", metric="pca", x="ignored_metadata_axis", color="condition", facet="tissue", order="treated|control%3Cscript%3E", scale="linear")
        doc=render_report(self.run, [], analyses, [], {}, [request]).read_text()
        dest=self.run/"analyses/main/exploratory/configured_plots"
        self.assertEqual(len(list(dest.glob("*.svg"))), 2)
        settings=json.loads(next(dest.glob("*_settings.json")).read_text())
        self.assertEqual(settings["axes"], ["PC1","PC2"])
        self.assertEqual(settings["color_order"], ["treated", "control<script>"])
        self.assertEqual(settings["x_setting_ignored"], "ignored_metadata_axis")
        self.assertEqual(settings["sample_count"], 3)
        points={r["sample_id"]: r for p in dest.glob("*.tsv") for r in self.read_table(p)}
        self.assertEqual(points["c"]["condition"], "control<script>")
        self.assertEqual(points["a"]["tissue"], "root")
        self.assertEqual((points["a"]["PC1"],points["a"]["PC2"]), ("-2.5","0.125"))
        colors={}
        for path in dest.glob("*.svg"):
            text=path.read_text()
            self.assertIn("75.0% variance", text)
            self.assertNotIn("<script>", text)
            tree=ElementTree.parse(path)
            for circle in tree.findall('.//{http://www.w3.org/2000/svg}circle'):
                title=circle.find('{http://www.w3.org/2000/svg}title')
                if title is not None:
                    colors[title.text.split(" | ")[0]]=circle.attrib["fill"]
        self.assertEqual(colors["a"],colors["c"])
        self.assertNotEqual(colors["a"],colors["b"])
        self.assertIn("analyses/main/exploratory/PCA.pdf",doc)
        self.assertIn("configured_plots",doc)
        self.assertIn("x setting is ignored",doc)
        self.assertEqual(self.read_table(self.run/"analysis_plot_status.tsv")[0]["status"],"ok")

    def test_mds_negative_coordinates_are_not_logged_or_refit(self):
        analyses=self.setup_coordinates()
        render_report(self.run, [], analyses, [], {}, [dict(plot_id="mds",scope="analysis",analysis_id="main",metric="mds",color="tissue")])
        dest=self.run/"analyses/main/exploratory/configured_plots"
        data=self.read_table(next(dest.glob("*.tsv")))
        self.assertEqual({r["sample_id"]: r["MDS2"] for r in data}, {"a":"-5","b":"3","c":"2"})
        settings=json.loads(next(dest.glob("*_settings.json")).read_text())
        self.assertEqual(settings["scale"],"linear")
        self.assertEqual(settings["axes"],["MDS1","MDS2"])

    def test_requested_plot_errors_are_visible_and_do_not_erase_qc(self):
        analyses=self.setup_coordinates()
        requests=[dict(plot_id="wrong_metric",scope="analysis",analysis_id="main",metric="volcano"),
                  dict(plot_id="wrong_scale",scope="analysis",analysis_id="main",metric="pca",scale="log10"),
                  dict(plot_id="wrong_column",scope="analysis",analysis_id="main",metric="pca",color="made_up")]
        doc=render_report(self.run,[self.metric("raw",0)],analyses,[],{},requests).read_text()
        statuses=self.read_table(self.run/"analysis_plot_status.tsv")
        self.assertTrue(all(r["status"]=="failed" for r in statuses))
        self.assertIn("unsupported analysis metric",doc)
        self.assertIn("missing metadata column",doc)
        self.assertIn("Sample metrics and availability",doc)
        self.assertFalse((self.run/"analyses/main/exploratory/configured_plots").exists())

    def test_unknown_coordinate_ids_fail_instead_of_positional_join(self):
        analyses=self.setup_coordinates()
        self.table("analyses/main/exploratory/pca_scores.tsv",[dict(sample_id="unknown",PC1=1,PC2=2)])
        render_report(self.run,[],analyses,[],{},[dict(plot_id="pca",scope="analysis",analysis_id="main",metric="pca")])
        status=self.read_table(self.run/"analysis_plot_status.tsv")[0]
        self.assertEqual(status["status"],"failed")
        self.assertIn("missing from recorded metadata",status["reason"])

    def test_unified_status_records_every_qc_and_analysis_request(self):
        analyses = self.setup_coordinates()
        requests = [
            dict(plot_id="qc_ok", scope="quant_qc", metric="gene_count_sum", x="tissue"),
            dict(plot_id="pca_ok", scope="analysis", analysis_id="main", metric="pca", color="tissue"),
            dict(plot_id="bad_metric", scope="quant_qc", metric="fiction"),
            dict(plot_id="bad_scale", scope="quant_qc", metric="gene_count_sum", scale="sqrt"),
            dict(plot_id="bad_column", scope="quant_qc", metric="gene_count_sum", facet="does_not_exist"),
            dict(plot_id="bad_scope", scope="unrecognized", metric="gene_count_sum"),
            dict(plot_id="duplicate", scope="quant_qc", metric="gene_count_sum"),
            dict(plot_id="duplicate", scope="analysis", analysis_id="main", metric="mds"),
            dict(plot_id="unavailable", scope="quant_qc", metric="mapping_rate_pct"),
        ]
        doc = render_report(self.run, [self.metric("a", 0)], analyses, [], {}, requests).read_text()
        outcomes = self.read_table(self.run/"plot_status.tsv")
        self.assertEqual(len(outcomes), len(requests))
        self.assertEqual([o["request_row"] for o in outcomes], [str(i) for i in range(1, len(requests)+1)])
        self.assertEqual([o["status"] for o in outcomes], ["ok", "ok", "failed", "failed", "failed", "failed", "failed", "failed", "blocked"])
        self.assertTrue(all("duplicate plot_id" in o["reason"] for o in outcomes if o["plot_id"]=="duplicate"))
        self.assertIn("No available mapping_rate_pct", outcomes[-1]["reason"])
        backwards = self.read_table(self.run/"analysis_plot_status.tsv")
        self.assertEqual([(o["plot_id"],o["status"]) for o in backwards], [("pca_ok","ok"),("duplicate","failed")])
        self.assertIn("All plot request outcomes TSV", doc)
        self.assertIn("unsupported QC scale", doc)
        self.assertTrue(list((self.run/"00_quant_qc/plots").glob("*.svg")))

    def test_missing_plot_id_and_scope_fail_explicitly(self):
        requests = [dict(scope="quant_qc", metric="gene_count_sum"),
                    dict(plot_id="no_scope", metric="gene_count_sum")]
        render_report(self.run, [self.metric("a", 0)], [], [], {}, requests)
        outcomes = self.read_table(self.run/"plot_status.tsv")
        self.assertEqual([o["status"] for o in outcomes], ["failed", "failed"])
        self.assertIn("plot_id is required", outcomes[0]["reason"])
        self.assertIn("unsupported plot scope", outcomes[1]["reason"])

    def test_initial_report_records_pending_coordinates(self):
        render_report(self.run,[],[dict(analysis_id="main")],[],{},[dict(plot_id="pca",scope="analysis",analysis_id="main",metric="pca")])
        status=self.read_table(self.run/"analysis_plot_status.tsv")[0]
        self.assertEqual(status["status"],"blocked")
        self.assertIn("Coordinates unavailable",status["reason"])


if __name__ == "__main__":
    unittest.main()
