#!/usr/bin/env python3
"""Build a read-only-source Piden example: complete QC and exploratory PCA, no DEG."""
from __future__ import annotations
import argparse,csv,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"post_analysis"))
from post.cli import write_tsv
from post.schema import INPUT_COLUMNS,SAMPLE_COLUMNS,ANALYSIS_COLUMNS,CONTRAST_COLUMNS,ANNOTATION_COLUMNS,PLOT_COLUMNS,path_value

def rows(path):
    raw=path.read_bytes()
    try:text=raw.decode("utf-8-sig")
    except UnicodeDecodeError:text=raw.decode("cp949")
    reader=csv.DictReader(text.splitlines(),delimiter="\t")
    return [{k:v for k,v in row.items() if k} for row in reader]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",default="/mnt/c/AI/Research/Piden_transcriptome/htseq-stat")
    p.add_argument("--project-dir",default="post/piden_htseq_example")
    p.add_argument("--min-assigned",type=int,default=10000)
    p.add_argument("--min-detected",type=int,default=1000)
    args=p.parse_args()
    source=path_value(args.source,Path.cwd());dest=path_value(args.project_dir,Path.cwd())
    if dest.exists():raise SystemExit(f"Refusing to replace existing example project: {dest}")
    if min(args.min_assigned,args.min_detected)<1:raise SystemExit("positive exploratory thresholds required")
    groups=rows(source/"sample_groups.tsv")
    metadata={r["sample_id_original"]:r for r in groups}
    if len(metadata)!=len(groups):raise SystemExit("duplicate source metadata ID")
    with (source/"piden2trima_all_htseq.csv").open(encoding="utf-8-sig",newline="") as h:
        r=csv.reader(h);header=next(r);names=header[1:]
        sums={n:0 for n in names};detected={n:0 for n in names};ngenes=0
        for row in r:
            ngenes+=1
            for name,text in zip(names,row[1:]):
                value=int(text);sums[name]+=value;detected[name]+=int(value>0)
    if len(names)!=len(set(names)):raise SystemExit("duplicate count columns")
    # Explicit exploratory eligibility rule; no biological effect is used.
    known_positive=[n for n in names if n in metadata and sums[n]>0]
    retained_sum={n:0 for n in known_positive}
    index={n:names.index(n)+1 for n in known_positive}
    with (source/"piden2trima_all_htseq.csv").open(encoding="utf-8-sig",newline="") as h:
        rr=csv.reader(h);next(rr)
        for values in rr:
            vals={n:int(values[index[n]]) for n in known_positive}
            if sum(vals[n]*1000000/sums[n]>=1 for n in known_positive)>=min(2,len(known_positive)):
                for n in known_positive:retained_sum[n]+=vals[n]
    dest.mkdir(parents=True)
    inp=dict(input_id="piden_counts",source_project="piden_combined",assay_id="matsutake_genes",
             format="htseq_matrix",counts_path=str(source/"piden2trima_all_htseq.csv"),stats_path="",
             stats_format="htseq_legacy",gene_universe="unknown",run_metadata_path="",
             reference_id="Tricholoma_matsutake_declared",annotation_id="Tricholoma_matsutake_16046_declared",
             count_unit="unknown",assignment_mode="exclusive",provenance_level="declared")
    write_tsv(dest/"inputs.tsv",[inp],INPUT_COLUMNS)
    sample_rows=[]
    for name in names:
        meta=metadata.get(name);sets=[];reason=""
        if meta and sums[name]>0 and retained_sum[name]>0:sets.append("positive_known")
        if not meta:reason="No row in supplied sample_groups.tsv; biological metadata unverified"
        elif sums[name]==0:reason="Raw gene count sum is zero; retained in quantification QC"
        elif retained_sum[name]==0:reason="No retained counts at CPM >=1 in at least 2 positive metadata-known samples; QC retained"
        if meta and meta.get("tissue")=="root" and meta.get("source")=="2023_shiro_colonized" and sums[name]>=args.min_assigned and detected[name]>=args.min_detected:
            sets.append("root_signal")
        sid=meta["sample_id"] if meta else name
        row=dict(sample_id=sid,input_id="piden_counts",source_sample_id=name,
                 stats_sample_id=name,stats_source_path=str(source/"sample_metrics.tsv") if meta else "",
                 biological_unit_id=sid,library_id=sid,study_id=meta["source"] if meta else "metadata_missing",
                 sample_sets="|".join(sets),include_analysis="true" if sets else "false",
                 exclusion_reason=reason if reason else ("" if sets else "Outside exploratory sample sets"),
                 biological_unit_status="unverified",metadata_status="supplied" if meta else "missing",
                 group=meta.get("group","unknown") if meta else "unknown",
                 tissue=meta.get("tissue","unknown") if meta else "unknown",
                 group_sub=meta.get("group_sub","") if meta else "",
                 target_signal_status=("high_target_signal" if meta and meta.get("source")=="2023_shiro_colonized" and meta.get("tissue")=="root" and sums[name]>=args.min_assigned and detected[name]>=args.min_detected else ("low_or_unconfirmed_target_signal" if meta and meta.get("source")=="2023_shiro_colonized" and meta.get("tissue")=="root" else "outside_2023_root_evaluation")),
                 source_project_original=meta.get("project","") if meta else "",
                 raw_count_sum=sums[name],raw_genes_detected=detected[name])
        sample_rows.append(row)
    sample_fields=SAMPLE_COLUMNS+["biological_unit_status","metadata_status","group","tissue","group_sub",
                                 "source_project_original","raw_count_sum","raw_genes_detected","target_signal_status"]
    write_tsv(dest/"samples.tsv",sample_rows,sample_fields)
    analyses=[]
    for aid,sset in [("positive_overview","positive_known"),("root_signal_pca","root_signal")]:
        analyses.append(dict(analysis_id=aid,sample_set=sset,assay_id="matsutake_genes",analysis_mode="exploratory",
                             design_type="group",normalization="TMM",fdr_cut="0.05",lfc_cut="1",
                             effect_test="qlf",enrichment="none",adjusted_view="none",seed="123"))
    write_tsv(dest/"analyses.tsv",analyses,ANALYSIS_COLUMNS)
    write_tsv(dest/"contrasts.tsv",[],CONTRAST_COLUMNS)
    write_tsv(dest/"annotations.tsv",[],ANNOTATION_COLUMNS)
    plots=[]
    for metric in ("assignment_counter_rate_pct","supplied_gene_detection_pct","gene_count_sum"):
        plots.append(dict(plot_id=metric,scope="quant_qc",metric=metric,x="group",color="group",
                          facet="tissue",order="control|treated|unknown",scale="linear"))
    for aid,color in [("positive_overview","study_id"),("root_signal_pca","group_sub")]:
        for metric in ("pca","mds"):
            plots.append(dict(plot_id=aid+"_"+metric,scope="analysis",analysis_id=aid,
                              metric=metric,color=color,facet="",order="",scale="linear"))
    plots.append(dict(plot_id="positive_overview_signal_pca",scope="analysis",analysis_id="positive_overview",
                      metric="pca",color="target_signal_status",facet="study_id",order="high_target_signal|low_or_unconfirmed_target_signal|outside_2023_root_evaluation",scale="linear"))
    write_tsv(dest/"plots.tsv",plots,PLOT_COLUMNS)
    (dest/"config.yaml").write_text("schema_version: 1\nout_dir: output\nthreads: 1\nstrict_provenance: false\ntimeout_seconds: 3600\n",encoding="utf-8")
    summary=dict(input_samples=len(names),metadata_samples=len(groups),genes=ngenes,
                 zero_count_samples=sum(v==0 for v in sums.values()),
                 positive_known=sum("positive_known" in r["sample_sets"].split("|") for r in sample_rows),
                 root_signal=sum("root_signal" in r["sample_sets"].split("|") for r in sample_rows),
                 root_signal_rule=dict(tissue="root",study_id="2023_shiro_colonized",min_assigned=args.min_assigned,min_detected=args.min_detected),
                 positive_filter_rule="CPM >=1 in at least 2 metadata-known positive libraries; retain samples with nonzero filtered counts",
                 positive_filter_exclusions=[n for n in known_positive if retained_sum[n]==0],
                 biological_independence="unverified; no DEG",count_unit="unknown in supplied stats; do not call denominator raw sequenced reads",
                 missing_metadata=sorted(set(names)-set(metadata)))
    (dest/"example_selection.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    (dest/"README.md").write_text("# Piden HTSeq exploratory example\n\n"
        "All input columns are retained in raw QC. Original external files are read-only.\n"
        "Metadata come from sample_groups.tsv with its explicit label corrections; no names are parsed into biology.\n"
        "Quantification unit and full run options are not verified from the supplied summary. Native-counter assignment ratios use declared exclusive HTSeq semantics. Physical read/fragment assignment and sequencing mapping rates stay unavailable until their units/denominators are verified.\n"
        "Sample IDs stand in only as observation identifiers; biological independence remains unverified. No DEG is requested.\n"
        "positive_overview uses metadata-known libraries with raw sum > 0 and nonzero retained counts after CPM>=1 in at least two positive samples. Exclusions are in samples.tsv and example_selection.json.\n"
        f"root_signal_pca uses 2023 root samples with assigned >= {args.min_assigned} and detected genes >= {args.min_detected}; these are explicit exploratory thresholds, not infection or quality classifications.\n"
        "Gene detection percentages use the supplied matrix rows; reference annotation completeness is unverified.\n"
        "Target fungal read proportion and relative gene expression answer different questions.\n",encoding="utf-8")
    print(json.dumps(summary,ensure_ascii=False,indent=2))
    return 0
if __name__=="__main__":raise SystemExit(main())
