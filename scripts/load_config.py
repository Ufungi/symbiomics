#!/usr/bin/env python3
"""Read a symbiomics config.yaml and emit flat KEY=VALUE lines for --param.

Usage:
    load_config.py PATH [--keep-profile]

Parses the YAML config with a small stdlib-only parser (flat `key: value`
maps -- all symbiomics params are scalars), keeping only keys that match a
real Nextflow `params.*` declared in nextflow.config (plus `profile`), and
prints one `KEY=VALUE` per line.

Semantics:
  * non-scalar values appear now as empty (they are not used by this subset).
  * relative file paths resolve against the directory of this config file
    itself, so a project config under project/<name>/ is self-contained;
    absolute paths pass through unchanged.
  * unknown keys are silently ignored (never leak a bogus --param).
  * `project` (not a Nextflow param) maps to
    outdir=project/<project>/output PlaceTax-style; an explicit `outdir` in
    the config wins over `project`.

Exit codes: 0 on success, 1 if the config cannot be read/parsed.
"""
import json
import os
import re
import sys


def _known_params():
    """Return the set of params.* declared in nextflow.config."""
    nf = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "nextflow.config")
    keys = set()
    try:
        with open(nf, encoding="utf-8") as fh:
            in_params = False
            for line in fh:
                stripped = line.strip()
                if stripped == "params {":
                    in_params = True
                    continue
                if in_params:
                    if stripped == "}":
                        break
                    if stripped.startswith(("//", "/*")) or not stripped:
                        continue
                    key = stripped.split("=")[0].strip()
                    if key:
                        keys.add(key)
    except OSError:
        pass
    return keys


KNOWN = _known_params()
RESERVED = {"profile"}

_LINE = re.compile(r'^\s*([A-Za-z0-9_]+)\s*:\s*(.*?)\s*$')


def _scalar(raw):
    if raw.startswith(('"', "'")) and raw.endswith(raw[0]):
        return raw[1:-1]
    low = raw.lower()
    if low in ("true", "yes"):
        return json.dumps(True)
    if low in ("false", "no"):
        return json.dumps(False)
    if low in ("null", "~", "none"):
        return json.dumps(None)
    if re.fullmatch(r"[-+]?\d+", raw):
        return str(int(raw))
    if re.fullmatch(r"[-+]?(\d+\.\d*|\.\d+)([eE][-+]?\d+)?", raw):
        return raw
    return raw


def _render(scalar, repo_root):
    if repo_root and not scalar.startswith(("/", "$", "${")):
        candidate = os.path.join(repo_root, scalar)
        if os.path.exists(candidate) or re.search(r"[\\/]", scalar):
            return os.path.normpath(candidate)
    return scalar


def main():
    if len(sys.argv) < 2:
        print(__doc__.strip())
        return 0
    if sys.argv[1] in ("-h", "--help"):
        print(__doc__.strip())
        return 0

    config_path = sys.argv[1]
    keep_profile = len(sys.argv) > 2 and sys.argv[2] == "--keep-profile"

    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    cfg = {}
    with open(config_path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            m = _LINE.match(line)
            if not m:
                print("WARN: %s:%d: cannot parse line (skipped): %s"
                      % (config_path, lineno, stripped), file=sys.stderr)
                continue
            value = m.group(2)
            value = re.sub(r'\s+#.*$', '', value).strip()
            if not value:
                value = ''
            cfg[m.group(1)] = _scalar(value)
    # Not a Nextflow param: outdir=project/<project>/output, explicit outdir wins.
    if "outdir" not in cfg and cfg.get("project"):
        cfg["outdir"] = os.path.join(repo_root, "project", cfg["project"], "output")

    if not cfg:
        print("ERROR: no config entries parsed from %s" % config_path, file=sys.stderr)
        return 1

    config_dir = os.path.dirname(os.path.abspath(config_path))
    for key, value in cfg.items():
        if key not in KNOWN and key not in RESERVED:
            continue
        if key == "profile" and not keep_profile:
            continue
        print("%s=%s" % (key, _render(value, config_dir)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
