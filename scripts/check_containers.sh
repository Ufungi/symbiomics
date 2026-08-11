#!/usr/bin/env bash
# Verify that every container referenced in conf/containers.config actually
# resolves. This exists because a tag that looks plausible can simply not be
# published -- two of the original pins in this repo 404'd, and the failure only
# surfaced mid-run after several images had already been pulled.
#
# Run it in CI (weekly) and before tagging a release.
#   scripts/check_containers.sh          # check
#   scripts/check_containers.sh --lock   # also write assets/container_digests.lock.tsv
set -uo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="${REPO_DIR}/conf/containers.config"
LOCK="${REPO_DIR}/assets/container_digests.lock.tsv"
WRITE_LOCK=0
[[ "${1:-}" == "--lock" ]] && WRITE_LOCK=1

FAILED=0
declare -a ROWS

while IFS= read -r url; do
    [[ -z "$url" ]] && continue
    if [[ "$url" == https://* ]]; then
        code=$(curl -s -o /dev/null -w '%{http_code}' -I --max-time 30 "$url")
        if [[ "$code" == "200" ]]; then
            printf '[ ok ] %s\n' "$url"
            ROWS+=("${url}"$'\t'"http-200"$'\t'"$(date -u +%Y-%m-%d)")
        else
            printf '[FAIL] %s  -> HTTP %s\n' "$url" "$code"
            FAILED=1
        fi
    elif [[ "$url" == docker://* ]]; then
        ref="${url#docker://}"
        if [[ "$ref" == *"@sha256:"* ]]; then
            printf '[ ok ] %s  (digest-pinned)\n' "$url"
            ROWS+=("${url}"$'\t'"digest"$'\t'"$(date -u +%Y-%m-%d)")
            continue
        fi
        printf '[warn] %s  is tag-pinned, not digest-pinned\n' "$url"
        printf '       Images whose tag provenance is untrustworthy should be\n'
        printf '       digest-pinned. Resolve with: docker buildx imagetools inspect %s\n' "$ref"
        ROWS+=("${url}"$'\t'"tag-only"$'\t'"$(date -u +%Y-%m-%d)")
    fi
done < <(grep -oE "(https|docker)://[^'\"]+" "$CONFIG" | sort -u)

if (( WRITE_LOCK )); then
    {
        printf '# container\tpin_type\tchecked\n'
        printf '%s\n' "${ROWS[@]}"
    } > "$LOCK"
    echo "wrote $LOCK"
fi

if (( FAILED )); then
    echo
    echo "One or more containers do not resolve. Fix conf/containers.config before running."
    exit 1
fi
echo
echo "all containers resolve"
