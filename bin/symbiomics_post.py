#!/usr/bin/env python3
"""Thin entrypoint for the independent post-analysis command family."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"post_analysis"))
from post.cli import main
if __name__=="__main__":
    raise SystemExit(main())
