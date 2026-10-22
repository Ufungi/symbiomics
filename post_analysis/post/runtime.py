"""Resolve R without requiring the Nextflow environment."""
from __future__ import annotations
import os
import shutil
import subprocess
from pathlib import Path

def resolve_rscript(config=None):
    config=config or {}
    explicit=config.get("rscript") or os.environ.get("SYMBIOMICS_RSCRIPT")
    candidates=[]
    if explicit:
        p=Path(str(explicit)).expanduser()
        if not p.is_file():
            raise ValueError(f"Explicit Rscript does not exist: {p}")
        candidates.append([str(p),"--vanilla"])
    else:
        r=shutil.which("Rscript")
        if r: candidates.append([r,"--vanilla"])
        conda=shutil.which("conda")
        if not conda:
            for prefix in ("miniconda3","miniforge3","anaconda3"):
                p=Path.home()/prefix/"bin"/"conda"
                if p.is_file(): conda=str(p);break
        if conda:
            candidates.append([conda,"run","--no-capture-output","-n","symbiomics-post","Rscript","--vanilla"])
    failures=[]
    for cmd in candidates:
        try:
            p=subprocess.run(cmd+["-e",'stopifnot(requireNamespace("jsonlite",quietly=TRUE),requireNamespace("edgeR",quietly=TRUE));cat(as.character(getRversion()))'],
                             capture_output=True,text=True,timeout=90)
            if p.returncode==0:
                return cmd, p.stdout.strip()
            failures.append(p.stderr.strip()[-1000:])
        except (OSError,subprocess.TimeoutExpired) as e: failures.append(str(e))
    raise ValueError("R with jsonlite and edgeR is unavailable. Create the isolated environment: "
                     "conda env create -f envs/post.yml. "+ " | ".join(failures))
