#!/usr/bin/env bash
# symbiomics preflight -- environment checks with an exact remedy for each failure.
# Never exits non-zero on warnings; only hard failures (Java, Nextflow) do that.
set -uo pipefail

: "${SYMBIOMICS_DB_DIR:=/data/db/eukannot}"
DB_DIR="${1:-$SYMBIOMICS_DB_DIR}"

FAILED=0
ok()   { printf '[ ok ] %s\n' "$*"; }
warn() { printf '[warn] %s\n' "$*"; }
fail() { printf '[FAIL] %s\n' "$*"; FAILED=1; }
fix()  { printf '       fix: %s\n' "$*"; }

echo "=== symbiomics preflight ==="

# ---------------------------------------------------------------- Java
JAVA_BIN="${JAVA_CMD:-$(command -v java || true)}"
if [[ -z "$JAVA_BIN" ]]; then
    fail "no java found"
    fix "conda create -n nextflow -c conda-forge -c bioconda nextflow openjdk=23"
else
    JV=$("$JAVA_BIN" -version 2>&1 | head -1 | sed -E 's/.*version "([0-9]+).*/\1/')
    if [[ "$JV" =~ ^[0-9]+$ ]] && (( JV >= 17 && JV <= 24 )); then
        ok "java $JV via $JAVA_BIN"
    else
        fail "java $JV via $JAVA_BIN -- Nextflow needs 17..24"
        fix "use scripts/symbiomics (sets JAVA_CMD to the conda env JDK)"
    fi
fi

# ------------------------------------------------------------ Nextflow
if command -v nextflow >/dev/null 2>&1; then
    NFV=$(nextflow -v 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
    ok "nextflow ${NFV:-unknown}"
else
    fail "nextflow not on PATH"
    fix "conda activate nextflow   (or use scripts/symbiomics)"
fi

# --------------------------------------------------------- Containers
HAVE_RUNTIME=0
if command -v singularity >/dev/null 2>&1; then
    ok "singularity $(singularity --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"
    HAVE_RUNTIME=1
elif command -v apptainer >/dev/null 2>&1; then
    ok "apptainer $(apptainer --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"
    HAVE_RUNTIME=1
fi
if command -v docker >/dev/null 2>&1; then
    if docker ps >/dev/null 2>&1; then
        ok "docker $(docker --version | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1) (sudoless)"
    else
        warn "docker present but not usable without sudo"
    fi
    HAVE_RUNTIME=1
fi
if (( HAVE_RUNTIME == 0 )); then
    warn "no container runtime -- you must use -profile conda (see docs: incomplete tool coverage)"
fi

# --------------------------------------------------------------- GPU
if command -v nvidia-smi >/dev/null 2>&1; then
    GPU=$(nvidia-smi --query-gpu=name,memory.total,driver_version \
          --format=csv,noheader 2>/dev/null | head -1)
    if [[ -n "$GPU" ]]; then
        ok "GPU: $GPU"
        if command -v singularity >/dev/null 2>&1; then
            warn "GPU passthrough not exercised here; verify once with:"
            fix "singularity exec --nv docker://nvidia/cuda:12.2.0-base-ubuntu22.04 nvidia-smi"
        fi
    else
        warn "nvidia-smi present but reported no GPU"
    fi
else
    warn "no GPU -- Helixer/TMbed tracks will be unavailable"
fi

# ----------------------------------------------------------- Resources
NCPU=$(nproc 2>/dev/null || echo '?')
MEMG=$(free -g 2>/dev/null | awk '/^Mem:/{print $2}')
ok "cpus ${NCPU}, mem ${MEMG:-?} GB"

for D in "${NXF_WORK:-$DB_DIR/work}" "$DB_DIR"; do
    P="$D"
    while [[ ! -d "$P" && "$P" != "/" ]]; do P=$(dirname "$P"); done
    AVAIL=$(df -BG --output=avail "$P" 2>/dev/null | tail -1 | tr -dc '0-9')
    if [[ -z "$AVAIL" ]]; then
        warn "cannot stat free space for $D"
    elif (( AVAIL < 100 )); then
        warn "$D (on $P): only ${AVAIL} GB free -- large genomes need 3x genome size"
    else
        ok "$D (on $P): ${AVAIL} GB free"
    fi
done

# ----------------------------------------------------------- samtools
if command -v samtools >/dev/null 2>&1; then
    if samtools index 2>&1 | grep -q -- '-c'; then
        ok "samtools $(samtools --version 2>/dev/null | head -1 | awk '{print $2}') supports 'index -c' (CSI)"
    else
        warn "host samtools lacks 'index -c'; containers provide their own, so this is informational"
    fi
fi

# ----------------------------------------------------------- Databases
if [[ -d "$DB_DIR" ]]; then
    ok "db_dir $DB_DIR exists"
else
    warn "db_dir $DB_DIR does not exist yet"
    fix "scripts/symbiomics run . -entry download_dbs --db_dir $DB_DIR"
fi

# --------------------------------------------------------------- misc
NOFILE=$(ulimit -n)
if (( NOFILE < 4096 )); then
    warn "ulimit -n = $NOFILE; recommend >= 4096 for many-sample runs"
    fix "ulimit -n 4096"
else
    ok "ulimit -n = $NOFILE"
fi

if [[ "${LC_ALL:-${LANG:-}}" != C* ]]; then
    warn "locale is ${LC_ALL:-${LANG:-unset}}; some Perl tools (RepeatMasker) warn"
    fix "export LC_ALL=C"
else
    ok "locale ${LC_ALL:-$LANG}"
fi

command -v git >/dev/null 2>&1 && ok "git $(git --version | awk '{print $3}')" || warn "git not found"
command -v gh  >/dev/null 2>&1 && ok "gh CLI present" \
    || warn "gh CLI not installed -- create the GitHub repo in a browser and push over HTTPS"

echo
if (( FAILED )); then
    echo "preflight FAILED -- fix the [FAIL] lines above before running the pipeline."
    exit 1
fi
echo "preflight OK"
