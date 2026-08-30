#!/usr/bin/env bash
# symbiomics setup -- provision the conda environment that runs Nextflow.
#
# The pipeline itself can run tools two ways:
#   * -profile singularity|docker  (recommended -- full tool coverage)
#   * -profile conda               (per-tool envs in envs/*.yml, auto-created by
#                                   Nextflow at runtime; some tools have no conda
#                                   recipe and are skipped -- see conf/conda.config)
#
# Either way you need a conda with a "nextflow" env that ships Nextflow + a
# Java 17-24 JDK (the scripts/symbiomics launcher resolves both). Nothing else
# in the repo creates that env -- that is exactly what this script is for.
#
#   ./setup.sh                     create the runner env if missing (gap-fill)
#   ./setup.sh --update            sync the runner env to envs/runner.yml
#   ./setup.sh --pull-images       also pre-pull container images (singularity)
#   ./setup.sh --check             only verify an existing setup, never install
#
# Exit codes: 0 = all good, 1 = hard failure (conda missing / env broken).
set -uo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNNER_YML="${REPO_DIR}/envs/runner.yml"
HOST_CONF="${REPO_DIR}/conf/host.sh"

NF_ENV="${SYMBIOMICS_NF_ENV:-nextflow}"

DO_UPDATE=false
DO_PULL=false
DO_CHECK=false
for arg in "$@"; do
    case "$arg" in
        --update)      DO_UPDATE=true ;;
        --pull-images) DO_PULL=true ;;
        --check)       DO_CHECK=true ;;
        -h|--help)
            sed -n '2,14p' "$0" | sed 's/^# \{0,1\}//'
            exit 0 ;;
        *) echo "unknown option: $arg" >&2; exit 2 ;;
    esac
done

ok()   { printf '[ ok ] %s\n' "$*"; }
warn() { printf '[warn] %s\n' "$*"; }
fail() { printf '[FAIL] %s\n' "$*"; FAILED=1; }

# fail() references FAILED, so initialise once up front (also the final guard).
FAILED=0

echo "== symbiomics setup =="

# ------------------------------------------------------------ locate conda
CONDA_BASE=""
if [[ -n "${CONDA_EXE:-}" ]]; then
    CONDA_BASE="$( "$CONDA_EXE" info --base 2>/dev/null || true )"
fi
if [[ -z "$CONDA_BASE" ]]; then
    CONDA_BIN="$(command -v conda || true)"
    if [[ -n "$CONDA_BIN" ]]; then
        CONDA_BASE="$( "$CONDA_BIN" info --base 2>/dev/null || true )"
    fi
fi
CONDA_SH=""
if [[ -n "$CONDA_BASE" && -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
    CONDA_SH="$CONDA_BASE/etc/profile.d/conda.sh"
fi

if [[ -z "$CONDA_SH" ]]; then
    echo "[FAIL] could not locate a conda installation." >&2
    echo "       Install Miniforge/Miniconda first:" >&2
    echo "         https://github.com/conda-forge/miniforge" >&2
    echo "       Then re-run this script from a shell where 'conda' is on PATH," >&2
    echo "       or set CONDA_EXE to your conda binary." >&2
    exit 1
fi
echo "conda base: $CONDA_BASE"
# conda.sh is not `set -u` safe; source it with nounset off.
set +u
# shellcheck disable=SC1090
source "$CONDA_SH"
set -u

if ! command -v conda >/dev/null 2>&1; then
    echo "[FAIL] conda not usable after sourcing $CONDA_SH" >&2
    exit 1
fi

# --------------------------------------------------------- runner env yml
if [[ ! -f "$RUNNER_YML" ]]; then
    echo "[FAIL] missing $RUNNER_YML -- is this a symbiomics checkout?" >&2
    exit 1
fi

env_exists() { conda env list 2>/dev/null | awk '$1=="'"$NF_ENV"'" && NF>=1{found=1} END{exit !found}'; }

if [[ "$DO_CHECK" == true ]]; then
    echo "mode: CHECK (verify only)"
else
    if [[ "$DO_UPDATE" == true ]]; then
        echo "mode: UPDATE (sync runner env to $RUNNER_YML)"
    else
        echo "mode: INSTALL (create if missing, gap-fill only)"
    fi
fi

# ------------------------------------------------------- create/update/env
if env_exists; then
    if [[ "$DO_CHECK" == true ]]; then
        echo "[runner] env '$NF_ENV' exists (skipping checks only applies to install -- updating is not run)"
    elif [[ "$DO_UPDATE" == true ]]; then
        echo "[runner] env '$NF_ENV' exists -- updating from $RUNNER_YML"
        conda env update -n "$NF_ENV" -f "$RUNNER_YML" --prune
    else
        echo "[runner] env '$NF_ENV' exists -- skipping (use --update to sync to $RUNNER_YML)"
    fi
else
    if [[ "$DO_CHECK" == true ]]; then
        echo "[FAIL] env '$NF_ENV' does not exist (run ./setup.sh without --check to create it)" >&2
        exit 1
    fi
    echo "[runner] creating env '$NF_ENV' from $RUNNER_YML"
    conda env create -n "$NF_ENV" -f "$RUNNER_YML"
fi

# ------------------------------------------------------------- verify env
echo
echo "== Verification =="

JAVA_OUT="$( conda run -n "$NF_ENV" bash -c 'command -v java && java -version 2>&1 | head -1' 2>/dev/null || true )"
JAVA_VER="$( printf '%s' "$JAVA_OUT" | grep -oE 'version "[0-9]+' | tr -dc '0-9' )"
if [[ "$JAVA_VER" =~ ^[0-9]+$ ]] && (( JAVA_VER >= 17 && JAVA_VER <= 24 )); then
    ok "java $JAVA_VER in env '$NF_ENV'"
else
    fail "no java 17..24 in env '$NF_ENV' (found: ${JAVA_VER:-none})"
fi

NFV="$( conda run -n "$NF_ENV" nextflow -v 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true )"
if [[ -n "$NFV" ]]; then
    ok "nextflow $NFV in env '$NF_ENV'"
else
    fail "nextflow not present in env '$NF_ENV'"
fi

# -------------------------------------------------- write per-host override
if [[ "$DO_CHECK" != true ]] && [[ "$CONDA_SH" != "/home/genome/anaconda3/etc/profile.d/conda.sh" ]]; then
    mkdir -p "$(dirname "$HOST_CONF")"
    {
        echo "# Written by setup.sh -- do not edit by hand."
        echo "# Tells scripts/symbiomics where conda lives on THIS host."
        echo "SYMBIOMICS_CONDA_SH=\"$CONDA_SH\""
    } > "$HOST_CONF"
    echo "[host] wrote $HOST_CONF (conda.sh override for scripts/symbiomics)"
elif [[ "$DO_CHECK" == true && -f "$HOST_CONF" ]]; then
    echo "[host] override present: $(cat "$HOST_CONF" | head -1)"
fi

# ------------------------------------------------ optional: pull images
if [[ "$DO_PULL" == true ]]; then
    if command -v singularity >/dev/null 2>&1 || command -v apptainer >/dev/null 2>&1; then
        echo
        echo "== Pre-pulling container images (singularity) =="
        "${REPO_DIR}/scripts/check_containers.sh"
        while IFS= read -r url; do
            [[ -z "$url" || "$url" != docker://* ]] && continue
            img="${url#docker://}"
            echo "[pull] $img"
            singularity pull "docker://${img}" || warn "pull failed for $img (continue)"
        done < <(grep -oE 'docker://[^"'"'"']+' "${REPO_DIR}/conf/containers.config" | sort -u)
    else
        warn "--pull-images needs singularity/apptainer; skipping (use -profile conda instead)"
    fi
fi

echo
if (( FAILED )); then
    echo "setup FAILED -- see the [FAIL] lines above."
    exit 1
fi

echo "== Next steps =="
echo "  1. scripts/symbiomics preflight            # environmental health check"
echo "  2. (optional) download functional-annotation DBs:"
echo "       scripts/symbiomics run . -entry download_dbs --db_dir /data/db/eukannot"
echo "  3. scripts/symbiomics run . -profile singularity,local64 --input samplesheet.tsv \\"
echo "       --genome genome.fasta --taxon plant --outdir results"
echo "  See docs/installation.md for the full manual."
