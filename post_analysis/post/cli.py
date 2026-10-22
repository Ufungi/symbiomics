"""CLI lifecycle for table-driven RNA-seq post-analysis."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
import shutil
import socket
import subprocess
import sys
import uuid
from datetime import datetime, timezone
from pathlib import Path

from .schema import (COLUMNS, read_table, safe_id, path_value, validate_import_tables,
                     analysis_defaults, validate_analysis)
from .runtime import resolve_rscript

REPO=Path(__file__).resolve().parents[2]
CONFIG_KEYS={"schema_version","out_dir","threads","strict_provenance","rscript","timeout_seconds"}

def dump_json(path, obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+".tmp")
    tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=str)+"\n",encoding="utf-8")
    os.replace(tmp,path)

def write_tsv(path, rows, fields=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fields=fields or list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w",encoding="utf-8",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",extrasaction="ignore",lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in row.items()})

def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):h.update(block)
    return h.hexdigest()

def config_load(project):
    p=project/"config.yaml"
    defaults={"schema_version":1,"out_dir":"output","threads":1,"strict_provenance":False,"timeout_seconds":3600}
    if p.exists():
        try:import yaml
        except ImportError:raise ValueError("config.yaml needs PyYAML; install python3-yaml or use the post environment")
        class UniqueLoader(yaml.SafeLoader):
            pass
        def unique_mapping(loader, node, deep=False):
            result={}
            for key_node, value_node in node.value:
                key=loader.construct_object(key_node, deep=deep)
                if key in result:raise ValueError(f"duplicate config key: {key}")
                result[key]=loader.construct_object(value_node, deep=deep)
            return result
        UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,unique_mapping)
        try:data=yaml.load(p.read_text(encoding="utf-8"),Loader=UniqueLoader) or {}
        except yaml.YAMLError as exc:raise ValueError(f"invalid config.yaml: {exc}") from exc
        if not isinstance(data,dict):raise ValueError("config.yaml must be a mapping")
        extra=set(data)-CONFIG_KEYS
        if extra:raise ValueError(f"unknown config keys: {sorted(extra)}")
        defaults.update(data)
    if isinstance(defaults["schema_version"],bool) or defaults["schema_version"]!=1:raise ValueError("unsupported schema_version")
    for k in ("threads","timeout_seconds"):
        if isinstance(defaults[k],bool) or not isinstance(defaults[k],int) or defaults[k]<1:
            raise ValueError(f"config {k} must be a positive integer")
    if not isinstance(defaults["strict_provenance"],bool):raise ValueError("strict_provenance must be a boolean")
    if not isinstance(defaults["out_dir"],str) or not defaults["out_dir"].strip():
        raise ValueError("config out_dir must be a nonempty path")
    if "rscript" in defaults and (not isinstance(defaults["rscript"],str) or not defaults["rscript"].strip()):
        raise ValueError("config rscript must be a nonempty path")
    return defaults

def project_path(args):
    if args.project_dir:
        return Path(args.project_dir).expanduser().resolve()
    if not args.project:raise ValueError("-p NAME or --project-dir PATH is required")
    return REPO/"post"/safe_id(args.project,"project")

def scaffold(args):
    project=project_path(args)
    if project.exists():raise ValueError(f"project already exists: {project}")
    project.mkdir(parents=True)
    template=REPO/"assets"/"post"
    for src in sorted(template.iterdir()):
        if src.is_file():shutil.copyfile(src,project/src.name)
    print(f"Created {project}\nEdit the TSV tables, then run: symbiomics post validate --project-dir '{project}'")
    return 0

def load_deferred(project,name):
    try:return read_table(project/(name+".tsv"),name,optional=True),None
    except ValueError as e:return [],str(e)

def check_strict_provenance(inputs):
    bad=[r["input_id"] for r in inputs if not r.get("reference_sha256") or not r.get("annotation_sha256")]
    if bad:raise ValueError(f"strict provenance needs reference/annotation hashes: {bad}")

def selected_samples(assay,analysis):
    return [s for s in assay["samples"] if str(s.get("include_analysis","true")).lower()=="true"
            and analysis["sample_set"] in str(s.get("sample_sets","")).split("|")]

def export_subset(counts_path,samples,out_path):
    ids=[s["sample_id"] for s in samples]
    with Path(counts_path).open(encoding="utf-8-sig",newline="") as src:
        r=csv.reader(src,delimiter="\t");header=next(r)
        try:indices=[header.index(s) for s in ids]
        except ValueError as e:raise ValueError("resolved sample missing from imported matrix") from e
        with Path(out_path).open("w",encoding="utf-8",newline="") as dst:
            w=csv.writer(dst,delimiter="\t",lineterminator="\n");w.writerow(["gene_id"]+ids)
            for row in r:w.writerow([row[0]]+[row[i] for i in indices])

def render(run, imported, analyses, statuses, manifest,plots):
    from .report import render_report
    return render_report(run,imported.get("metrics",[]),analyses,statuses,manifest,plots)

def resolve_annotations(rows,project):
    out=[]
    for row in rows:
        safe_id(row.get("annotation_id",""),"annotation_id")
        if row.get("kind") not in {"gene_go","gene_pathway","gene_labels"}:
            raise ValueError(f"unsupported annotation kind: {row.get('kind')}")
        r=dict(row)
        for k in ("path","id_map_path","universe_path"):
            if r.get(k):r[k]=str(path_value(r[k],project))
        out.append(r)
    return out

def execute_project(args):
    from .importers import import_inputs
    from .qc import emit_qc
    project=project_path(args)
    cfg=config_load(project)
    output=path_value(cfg["out_dir"],project)
    output.mkdir(parents=True,exist_ok=True)
    lock=output/".post.lock"
    try:fd=os.open(str(lock),os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o600)
    except FileExistsError:raise ValueError(f"output is locked: {lock}; inspect its owner before removing a stale lock")
    with os.fdopen(fd,"w",encoding="utf-8") as h:
        json.dump({"pid":os.getpid(),"host":socket.gethostname(),"started":datetime.now(timezone.utc).isoformat()},h)
    run_id=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"_"+uuid.uuid4().hex[:8]
    run=output/"runs"/run_id;run.mkdir(parents=True)
    statuses=[];analyses=[];plots=[];imported={"metrics":[],"provenance":[],"assays":{}}
    manifest={"schema_version":1,"run_id":run_id,"project":str(project),"mode":args.command,
              "qc_only":bool(args.qc_only),"started_at":datetime.now(timezone.utc).isoformat(),
              "config":cfg,"files":{},"code_hashes":{},"status":"running","analyses":[]}
    for p in sorted((REPO/"post_analysis"/"post").rglob("*")):
        if p.is_file() and p.suffix in {".py",".R"}:
            manifest["code_hashes"][str(p.relative_to(REPO))]=digest(p)
    error=None
    try:
        resolved=run/"resolved";resolved.mkdir()
        for name in COLUMNS:
            src=project/(name+".tsv")
            if src.exists():
                shutil.copyfile(src,resolved/src.name);manifest["files"][src.name]=digest(src)
        if (project/"config.yaml").exists():
            shutil.copyfile(project/"config.yaml",resolved/"config.yaml")
            manifest["files"]["config.yaml"]=digest(project/"config.yaml")
        inputs=read_table(project/"inputs.tsv","inputs")
        samples=read_table(project/"samples.tsv","samples")
        validate_import_tables(inputs,samples)
        for row in inputs:
            for k in ("counts_path","stats_path","run_metadata_path","gene_map_path"):
                if row.get(k):row[k]=str(path_value(row[k],project))
        for row in samples:
            if row.get("stats_source_path"):row["stats_source_path"]=str(path_value(row["stats_source_path"],project))
        imported=import_inputs(inputs,samples,project,run)
        manifest["provenance"]=imported.get("provenance",[])
        manifest["warnings"]=imported.get("warnings",[])
        manifest["import_errors"]=imported.get("errors",[])
        emit_qc(imported,run)
        analyses,analysis_table_error=load_deferred(project,"analyses")
        analyses=[analysis_defaults(a) for a in analyses]
        plots,plot_error=load_deferred(project,"plots")
        statuses.append({"stage":"quant_qc","analysis_id":"","status":"ok","reason":"","path":"00_quant_qc"})
        render(run,imported,analyses,statuses,manifest,plots)
        print(f"Quantification QC: {run/'report.html'}",flush=True)
        if cfg["strict_provenance"]:
            try:check_strict_provenance(inputs)
            except ValueError as exc:statuses.append({"stage":"provenance","analysis_id":"","status":"failed","reason":str(exc),"path":"imported"})
        for msg in imported.get("errors",[]):
            statuses.append({"stage":"import","analysis_id":"","status":"failed","reason":str(msg),"path":"imported"})
        if plot_error:
            statuses.append({"stage":"plots","analysis_id":"","status":"failed","reason":plot_error,"path":""})
        if not args.qc_only:
            if analysis_table_error:raise ValueError(analysis_table_error)
            if not analyses:raise ValueError("No analyses declared. Use --qc-only or fill analyses.tsv")
            seen=set()
            for a in analyses:
                safe_id(a["analysis_id"],"analysis_id")
                if a["analysis_id"] in seen:raise ValueError(f"duplicate analysis_id: {a['analysis_id']}")
                seen.add(a["analysis_id"])
            contrasts,contrast_error=load_deferred(project,"contrasts")
            annotations,annotation_error=load_deferred(project,"annotations")
            if contrast_error:raise ValueError(contrast_error)
            unknown={r.get("analysis_id") for r in contrasts}-seen
            if unknown:raise ValueError(f"contrasts refer to unknown analyses: {sorted(unknown)}")
            # Environment failure still leaves the complete raw QC report.
            rcmd,rversion=resolve_rscript(cfg)
            manifest["rscript"]={"argv":rcmd,"r_version":rversion}
            for analysis in analyses:
                aid=analysis["analysis_id"];aout=run/"analyses"/aid;aout.mkdir(parents=True)
                status={"stage":"analysis","analysis_id":aid,"status":"failed","reason":"","path":str(aout.relative_to(run))}
                try:
                    assay=imported["assays"].get(analysis["assay_id"])
                    if not assay:raise ValueError(f"unknown assay: {analysis['assay_id']}")
                    if assay.get("errors"):raise ValueError("; ".join(map(str,assay["errors"])))
                    if cfg["strict_provenance"]:check_strict_provenance([i for i in inputs if i["assay_id"]==analysis["assay_id"]])
                    selected=selected_samples(assay,analysis)
                    if not selected:raise ValueError("sample_set selects no included samples")
                    units={s.get("count_unit","unknown") for s in selected}
                    if len(units)>1:raise ValueError(f"analysis mixes incompatible counting units: {sorted(units)}")
                    if analysis["analysis_mode"]=="inferential" and any(str(s.get("biological_unit_status","")).lower() in {"unverified","unknown"} for s in selected):
                        raise ValueError("biological independence has not been confirmed; use exploratory mode")
                    ac=[c for c in contrasts if c["analysis_id"]==aid]
                    validate_analysis(analysis,ac,selected,allow_empty_contrasts=args.command=="validate")
                    data=aout/"input";data.mkdir()
                    sample_path=data/"samples.tsv";count_path=data/"counts.tsv"
                    write_tsv(sample_path,selected)
                    export_subset(assay["counts_path"],selected,count_path)
                    membership=[{**s,"selected_for_analysis":s["sample_id"] in {r["sample_id"] for r in selected},
                                 "selection_reason":"" if s in selected else s.get("exclusion_reason") or "outside_sample_set"}
                                for s in assay["samples"]]
                    write_tsv(data/"sample_selection.tsv",membership)
                    aa=[];deferred_annotation_error=None
                    if analysis.get("enrichment","none")!="none":
                        try:
                            if annotation_error:raise ValueError(annotation_error)
                            aa=resolve_annotations([r for r in annotations if r.get("assay_id")==analysis["assay_id"]],project)
                        except ValueError as exc:deferred_annotation_error=str(exc)
                    job={"mode":"validate" if args.command=="validate" else "run","counts_path":str(count_path),
                         "samples_path":str(sample_path),"analysis":analysis,"contrasts":ac,"annotations":aa,"annotation_error":deferred_annotation_error,"out_dir":str(aout)}
                    if assay.get("salmon"):
                        job["salmon"]=[s for s in assay["salmon"] if s["sample_id"] in {r["sample_id"] for r in selected}]
                    dump_json(data/"job.json",job)
                    env=os.environ.copy()
                    env.update({"OMP_NUM_THREADS":str(cfg["threads"]),"OPENBLAS_NUM_THREADS":str(cfg["threads"])})
                    print(f"[{aid}] {analysis['analysis_mode']} / {len(selected)} samples",flush=True)
                    with (aout/"engine.log").open("w",encoding="utf-8") as log:
                        proc=subprocess.run(rcmd+[str(REPO/"post_analysis"/"post"/"R"/"engine.R"),str(data/"job.json")],
                                            stdout=log,stderr=subprocess.STDOUT,env=env,timeout=cfg["timeout_seconds"])
                    engine_status={}
                    if (aout/"status.json").exists():
                        engine_status=json.loads((aout/"status.json").read_text(encoding="utf-8"))
                    status["engine"]=engine_status
                    if proc.returncode:
                        raise ValueError(engine_status.get("reason") or engine_status.get("error") or f"R engine exited {proc.returncode}; see engine.log")
                    status["status"]=engine_status.get("status","ok")
                    if status["status"] not in {"ok","skipped"}:raise ValueError("R engine did not complete: "+status["status"])
                    status["reason"]=engine_status.get("reason","")
                    status["engine"]=engine_status
                    manifest["analyses"].append({"analysis_id":aid,"sample_count":len(selected),"settings":analysis,
                                                 "engine":engine_status,"job_sha256":digest(data/"job.json")})
                except Exception as exc:
                    status["reason"]=str(exc)
                    if not (aout/"status.json").exists():dump_json(aout/"status.json",status)
                    print(f"[{aid}] failed: {exc}",file=sys.stderr,flush=True)
                statuses.append(status)
                dump_json(run/"manifest.json",manifest)
                write_tsv(run/"stage_status.tsv",statuses,["stage","analysis_id","status","reason","path"])
                render(run,imported,analyses,statuses,manifest,plots)
    except (Exception, KeyboardInterrupt) as exc:
        error=str(exc) or "execution interrupted"
        statuses.append({"stage":"driver","analysis_id":"","status":"failed","reason":error,"path":""})
        print(f"ERROR ~ [post] {error}",file=sys.stderr)
    finally:
        # Rendering also validates requested charts. Initial pending plots are
        # expected; only the final render participates in completion status.
        try:
            render(run,imported,analyses,statuses,manifest,plots)
            status_path=run/"plot_status.tsv"
            if status_path.exists():
                with status_path.open(encoding="utf-8",newline="") as handle:
                    for result in csv.DictReader(handle,delimiter=chr(9)):
                        failed=result.get("status")=="failed"
                        unavailable=(result.get("status")=="blocked" and args.command=="run" and not args.qc_only)
                        if failed or unavailable:
                            statuses.append({"stage":"plot","analysis_id":result.get("analysis_id",""),
                                "status":"failed","reason":result.get("plot_id","")+": "+result.get("reason",""),
                                "path":result.get("path","")})
        except Exception as exc:
            statuses.append({"stage":"report","analysis_id":"","status":"failed","reason":str(exc),"path":"report.html"})
        failures=[s for s in statuses if s["status"] in {"failed","blocked"}]
        manifest["status"]=("partial_failure" if imported.get("metrics") else "failed") if failures else "ok"
        manifest["ended_at"]=datetime.now(timezone.utc).isoformat();manifest["stages"]=statuses
        dump_json(run/"manifest.json",manifest)
        write_tsv(run/"stage_status.tsv",statuses,["stage","analysis_id","status","reason","path"])
        try:render(run,imported,analyses,statuses,manifest,plots)
        except Exception as exc:
            manifest["status"]="partial_failure";manifest["report_error"]=str(exc);failures.append({"reason":str(exc)})
            dump_json(run/"manifest.json",manifest)
        session_files=list(run.glob("analyses/*/sessionInfo.txt"))
        if session_files:shutil.copyfile(session_files[0],run/"sessionInfo.txt")
        dump_json(output/"latest_run.json",{"run_id":run_id,"path":str(run),"status":manifest["status"]})
        if not failures:dump_json(output/"latest_success.json",{"run_id":run_id,"path":str(run),"status":"ok"})
        lock.unlink(missing_ok=True)
        print(f"Report: {run/'report.html'}\nStatus: {manifest['status']}",flush=True)
    return 1 if failures else 0

def merge(args):
    from .importers import import_inputs
    from .qc import emit_qc
    inputs_file=Path(args.inputs).resolve();samples_file=Path(args.samples).resolve()
    inputs=read_table(inputs_file,"inputs");samples=read_table(samples_file,"samples")
    validate_import_tables(inputs,samples)
    # Each path is relative to the table that declares it.
    for row in inputs:
        for k in ("counts_path","stats_path","run_metadata_path","gene_map_path"):
            if row.get(k):row[k]=str(path_value(row[k],inputs_file.parent))
    for row in samples:
        if row.get("stats_source_path"):row["stats_source_path"]=str(path_value(row["stats_source_path"],samples_file.parent))
    out=Path(args.out).resolve()
    if out.exists() and any(out.iterdir()):raise ValueError(f"merge output is not empty: {out}")
    out.mkdir(parents=True,exist_ok=True)
    result=import_inputs(inputs,samples,inputs_file.parent,out)
    emit_qc(result,out)
    dump_json(out/"merge_report.json",result)
    print(out/"merge_report.json")
    return int(bool(result.get("errors")))

def composition(args):
    from .composition import run
    path=run(args.libraries,args.out,args.primary,args.signal_organism,args.signal_ppm,args.min_aligned_pairs)
    print(path);return 0

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest="family",required=True)
    m=sub.add_parser("merge",help="merge compatible count assays and preserve quantification QC")
    m.add_argument("--inputs",required=True);m.add_argument("--samples",required=True);m.add_argument("--out",required=True)
    p=sub.add_parser("post",help="table-driven post-analysis")
    actions=p.add_subparsers(dest="command",required=True)
    for name in ("init","validate","run"):
        q=actions.add_parser(name)
        g=q.add_mutually_exclusive_group(required=True)
        g.add_argument("-p","--project");g.add_argument("--project-dir")
        if name!="init":q.add_argument("--qc-only",action="store_true")
    c=actions.add_parser("composition",help="multi-organism read composition from HTSeq special counters")
    c.add_argument("--libraries",required=True,help="TSV: sample_id, organism, file (HTSeq raw output per library and genome)")
    c.add_argument("--out",required=True);c.add_argument("--primary",required=True,help="organism whose counters define the input pair total")
    c.add_argument("--signal-organism",default=None);c.add_argument("--signal-ppm",type=float,default=1000.0)
    c.add_argument("--min-aligned-pairs",type=float,default=0.0)
    args=parser.parse_args(argv)
    try:
        if args.family=="merge":return merge(args)
        if args.family=="post" and args.command=="composition":return composition(args)
        if args.command=="init":return scaffold(args)
        return execute_project(args)
    except (ValueError,OSError) as exc:
        print(f"ERROR ~ [post] {exc}",file=sys.stderr);return 1

if __name__=="__main__":
    raise SystemExit(main())
