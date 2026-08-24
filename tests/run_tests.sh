#!/usr/bin/env bash
# All fast checks. Nothing here needs containers, a GPU, or real data.
#   tests/run_tests.sh
set -uo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"

FAILED=0
run() {
    echo
    echo "=============================================================="
    echo "  $1"
    echo "=============================================================="
    shift
    "$@" || FAILED=1
}

run "unit: samplesheet parsing and SE/PE inference" python3 tests/test_infer_layout.py
run "unit: policy engine"                           python3 tests/test_decide_strategy.py
run "unit: multi-genome mapping summary"            python3 tests/test_summarize_multi_genome_mapping.py
run "containers resolve"                            scripts/check_containers.sh

echo
echo "=============================================================="
echo "  workflow: -stub-run (every process wired and declaring outputs)"
echo "=============================================================="
NXF_WORK="${NXF_WORK:-/tmp/symbiomics_stub_work_$$}" \
    scripts/symbiomics run . -profile test -stub-run --outdir "/tmp/symbiomics_stub_out_$$" \
    > /tmp/symbiomics_stub_$$.log 2>&1
if grep -q "symbiomics completed" /tmp/symbiomics_stub_$$.log; then
    echo "  ok   stub run completed"
else
    echo "  FAIL stub run"
    tail -30 /tmp/symbiomics_stub_$$.log
    FAILED=1
fi
rm -rf "/tmp/symbiomics_stub_work_$$" "/tmp/symbiomics_stub_out_$$" "/tmp/symbiomics_stub_$$.log"

echo
echo "=============================================================="
echo "  workflow: -entry symbiont_mapping -stub-run"
echo "=============================================================="
NXF_WORK="${NXF_WORK:-/tmp/symbiomics_mapstub_work_$$}" \
    scripts/symbiomics run . -entry symbiont_mapping -profile test_mapping -stub-run \
    --outdir "/tmp/symbiomics_mapstub_out_$$" \
    > /tmp/symbiomics_mapstub_$$.log 2>&1
if grep -q "symbiomics completed" /tmp/symbiomics_mapstub_$$.log; then
    echo "  ok   stub run completed"
else
    echo "  FAIL stub run"
    tail -30 /tmp/symbiomics_mapstub_$$.log
    FAILED=1
fi
rm -rf "/tmp/symbiomics_mapstub_work_$$" "/tmp/symbiomics_mapstub_out_$$" "/tmp/symbiomics_mapstub_$$.log"

echo
if (( FAILED )); then
    echo "FAILED"
    exit 1
fi
echo "all checks passed"
echo
echo "Next, with containers and real data:"
echo "  scripts/symbiomics run . -profile test,singularity --outdir results_test"
