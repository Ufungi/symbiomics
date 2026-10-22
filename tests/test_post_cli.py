"""Fast end-to-end CLI lifecycle checks; no R or Nextflow is required."""
import csv,json,os,subprocess,sys,tempfile,unittest
from pathlib import Path
REPO=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO/"post_analysis"))
from post.cli import write_tsv
from post.schema import validate_import_tables, read_table, path_value, validate_analysis, analysis_defaults
CLI=REPO/"bin"/"symbiomics_post.py"

class Composition(unittest.TestCase):
    def test_composition_counts_and_groups(self):
        with tempfile.TemporaryDirectory() as t:
            t=Path(t)
            def w(name,rows): (t/name).write_text("".join(f"{a}\t{b}\n" for a,b in rows))
            w("a_pine.tsv",[("g1",60),("g2",20),("__no_feature",10),("__not_aligned",10)])
            w("a_fun.tsv",[("f1",5),("__no_feature",5),("__not_aligned",90)])
            w("b_pine.tsv",[("g1",80),("__not_aligned",20)])
            w("b_fun.tsv",[("f1",0),("__not_aligned",100)])
            (t/"libs.tsv").write_text("sample_id\torganism\tfile\n"+"".join(f"{s}\t{o}\t{f}\n" for s,o,f in
                [("a","pine","a_pine.tsv"),("a","fun","a_fun.tsv"),("b","pine","b_pine.tsv"),("b","fun","b_fun.tsv")]))
            p=subprocess.run([sys.executable,str(CLI),"post","composition","--libraries",str(t/"libs.tsv"),"--out",str(t/"o"),
                              "--primary","pine","--signal-organism","fun","--signal-ppm","50000","--min-aligned-pairs","85"],capture_output=True,text=True,timeout=40)
            self.assertEqual(p.returncode,0,p.stderr)
            rows={r["sample_id"]:r for r in csv.DictReader((t/"o"/"composition.tsv").open(),delimiter="\t")}
            self.assertEqual(rows["a"]["input_pairs"],"100");self.assertEqual(rows["a"]["pine_aligned_pairs"],"90")
            self.assertEqual(rows["a"]["fun_aligned_pct"],"10.0");self.assertEqual(rows["a"]["signal_group"],"high")
            self.assertEqual(rows["a"]["qc_pass"],"True");self.assertEqual(rows["b"]["qc_pass"],"False")
            self.assertEqual(rows["b"]["signal_group"],"low");self.assertTrue((t/"o"/"composition.svg").exists())

class PostCLI(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix="symbiomics-post-cli-")
        self.project=Path(self.tmp.name)/"project"
        p=self.call("post","init","--project-dir",str(self.project))
        self.assertEqual(p.returncode,0,p.stderr)
        with (self.project/"counts.tsv").open("w",newline="",encoding="utf-8") as h:
            w=csv.writer(h,delimiter="\t");w.writerow(["gene_id"]+[f"sample{i}" for i in range(1,7)])
            for g in range(30):w.writerow([f"g{g}"]+[10+g+i for i in range(6)])
    def tearDown(self):self.tmp.cleanup()
    def call(self,*args):
        return subprocess.run([sys.executable,str(CLI),*args],capture_output=True,text=True,timeout=40)
    def run_qc(self):
        return self.call("post","run","--project-dir",str(self.project),"--qc-only")
    def latest(self):
        ref=json.loads((self.project/"output"/"latest_run.json").read_text())
        return Path(ref["path"])
    def test_scaffold_refuses_overwrite(self):
        p=self.call("post","init","--project-dir",str(self.project))
        self.assertNotEqual(p.returncode,0)
        self.assertIn("already exists",p.stderr)
    def test_qc_without_analysis_and_immutable_runs(self):
        (self.project/"analyses.tsv").unlink();(self.project/"contrasts.tsv").unlink()
        p=self.run_qc();self.assertEqual(p.returncode,0,p.stderr)
        first=self.latest()
        self.assertTrue((first/"00_quant_qc"/"sample_metrics.tsv").exists())
        self.assertTrue((first/"report.html").exists())
        self.assertFalse((first/"00_quant_qc"/"00_quant_qc").exists())
        self.assertFalse((first/"imported"/"imported").exists())
        p=self.run_qc();self.assertEqual(p.returncode,0,p.stderr)
        self.assertNotEqual(first,self.latest());self.assertTrue(first.exists())
    def test_missing_r_preserves_qc_and_last_success(self):
        self.assertEqual(self.run_qc().returncode,0)
        old=json.loads((self.project/"output"/"latest_success.json").read_text())
        with (self.project/"config.yaml").open("a") as h:h.write("rscript: /nonexistent/post/Rscript\n")
        p=self.call("post","run","--project-dir",str(self.project))
        self.assertNotEqual(p.returncode,0)
        run=self.latest();manifest=json.loads((run/"manifest.json").read_text())
        self.assertEqual(manifest["status"],"partial_failure")
        self.assertTrue((run/"00_quant_qc"/"sample_metrics.tsv").exists())
        self.assertEqual(old,json.loads((self.project/"output"/"latest_success.json").read_text()))
    def test_strict_provenance_affects_qc_only(self):
        p=self.project/"config.yaml";p.write_text(p.read_text().replace("strict_provenance: false","strict_provenance: true"))
        run=self.run_qc();self.assertNotEqual(run.returncode,0)
        self.assertTrue((self.latest()/"00_quant_qc"/"sample_metrics.tsv").exists())
        status=json.loads((self.latest()/"manifest.json").read_text())
        self.assertTrue(any(s["stage"]=="provenance" and s["status"]=="failed" for s in status["stages"]))
    def test_output_lock_is_not_removed_by_nonowner(self):
        out=self.project/"output";out.mkdir()
        lock=out/".post.lock";lock.write_text('{"pid":123,"owner":"different"}')
        run=self.run_qc();self.assertNotEqual(run.returncode,0);self.assertTrue(lock.exists())
    def test_merge_preserves_qc_and_refuses_reuse(self):
        out=Path(self.tmp.name)/"merged"
        args=("merge","--inputs",str(self.project/"inputs.tsv"),"--samples",str(self.project/"samples.tsv"),"--out",str(out))
        p=self.call(*args);self.assertEqual(p.returncode,0,p.stderr)
        self.assertTrue((out/"merge_report.json").exists())
        self.assertTrue((out/"00_quant_qc"/"sample_metrics.tsv").exists())
        self.assertNotEqual(self.call(*args).returncode,0)
    def test_bad_override_and_hash_and_status_rejected(self):
        inputs=read_table(self.project/"inputs.tsv","inputs");samples=read_table(self.project/"samples.tsv","samples")
        samples[0]["count_unit_override"]="banana"
        with self.assertRaisesRegex(ValueError,"count_unit_override"):validate_import_tables(inputs,samples)
        samples[0]["count_unit_override"]="";inputs[0]["reference_sha256"]="not-a-hash"
        with self.assertRaisesRegex(ValueError,"SHA256"):validate_import_tables(inputs,samples)
        inputs[0]["reference_sha256"]="";samples[0]["biological_unit_status"]="Unverified"
        with self.assertRaisesRegex(ValueError,"biological_unit_status"):validate_import_tables(inputs,samples)
    def test_duplicate_config_is_rejected(self):
        with (self.project/"config.yaml").open("a") as h:h.write("threads: 2"+chr(10))
        result=self.run_qc();self.assertNotEqual(result.returncode,0)
        self.assertIn("duplicate config key",result.stderr)
    def test_validate_can_preview_coefficients_before_contrasts(self):
        analysis=analysis_defaults(read_table(self.project/"analyses.tsv","analyses")[0])
        samples=read_table(self.project/"samples.tsv","samples")
        validate_analysis(analysis,[],samples,allow_empty_contrasts=True)
        with self.assertRaisesRegex(ValueError,"needs declared contrasts"):
            validate_analysis(analysis,[],samples)
    def test_invalid_requested_plot_is_not_a_success(self):
        p=self.project/"plots.tsv"
        text=p.read_text().replace("assignment_rate_pct","unsupported_metric")
        p.write_text(text)
        result=self.run_qc();self.assertNotEqual(result.returncode,0)
        self.assertTrue((self.latest()/"00_quant_qc"/"sample_metrics.tsv").exists())
        self.assertFalse((self.project/"output"/"latest_success.json").exists())
    @unittest.skipUnless(os.environ.get("SYMBIOMICS_TEST_R")=="1","opt in for R CLI validation")
    def test_real_validate_before_declaring_contrasts(self):
        contrasts=self.project/"contrasts.tsv"
        contrasts.write_text(contrasts.read_text().splitlines()[0]+chr(10))
        result=self.call("post","validate","--project-dir",str(self.project))
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        run=self.latest()
        self.assertTrue(list(run.glob("analyses/*/design/coefficients.tsv")))
        self.assertFalse(list(run.glob("analyses/*/contrasts/*/DEG.tsv")))
    def test_windows_paths_resolve_under_wsl(self):
        if os.name!="nt":
            self.assertEqual(str(path_value("C:"+chr(92)+"AI"+chr(92)+"a.tsv",self.project)),"/mnt/c/AI/a.tsv")
    def test_cross_assay_aggregate_id_rejected(self):
        inputs=read_table(self.project/"inputs.tsv","inputs");samples=read_table(self.project/"samples.tsv","samples")
        inputs.append({**inputs[0],"input_id":"input2","assay_id":"assay2"})
        for n in (0,1):
            samples[n].update(technical_merge_id="aggregate",technical_partition_id=str(n),technical_disjoint="true")
        samples[1]["input_id"]="input2"
        with self.assertRaisesRegex(ValueError,"across assays"):validate_import_tables(inputs,samples)

if __name__=="__main__":unittest.main()
