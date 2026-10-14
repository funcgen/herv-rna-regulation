#!/usr/bin/env bash
# =============================================================================
# 02_rbp / sensitivity_strict : build the six strict p <= 1e-5 datasets (Stage 11)
# =============================================================================
#   accepted primary fimo.tsv (Stage 3 observed + Stage 7 shuffles)
#     -> streaming row filter, float(field 8) <= 1e-5
#     -> results/sensitivity_strict/<dataset>/fimo.tsv
#
# FIMO is NOT rerun. One filter implementation, invoked six times.
#
# Usage:  ./run_filter_fimo_strict.sh [--out-dir DIR]
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../../../config/lib/common.sh"

MODULE_DIR="${RELEASE_ROOT}/02_rbp"
MANIFEST="${MODULE_DIR}/results/null_model/fimo_dataset_manifest.tsv"
OUT_ROOT="${MODULE_DIR}/results/sensitivity_strict"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) OUT_ROOT="$2"; shift 2 ;;
    *) die "unknown argument: $1" ;;
  esac
done

mkdir -p "${MODULE_DIR}/logs" "$OUT_ROOT"
export NCRNA_LOG_FILE="${MODULE_DIR}/logs/filter_fimo_strict.$(date +%Y%m%dT%H%M%S).log"

[[ -f "$MANIFEST" ]] || die "Stage 7 dataset manifest missing: ${MANIFEST}"

# Acceptance expectations. Hard-coded ONLY as validation targets; the values
# written to the outputs are always derived from the files themselves.
declare -A EXPECT_PRIMARY=( [observed]=716075 [shuffle01]=616810 [shuffle02]=616790
                            [shuffle03]=615771 [shuffle04]=616507 [shuffle05]=614666 )
declare -A EXPECT_STRICT=(  [observed]=13702  [shuffle01]=12224  [shuffle02]=12138
                            [shuffle03]=12352  [shuffle04]=12204  [shuffle05]=12235 )
declare -A EXPECT_BOUNDARY=([observed]=34     [shuffle01]=30     [shuffle02]=13
                            [shuffle03]=24     [shuffle04]=28     [shuffle05]=19 )

log "system=${NCRNA_SYSTEM} release_root=${RELEASE_ROOT}"
log "manifest = ${MANIFEST}"
log "out root = ${OUT_ROOT}"
log "threshold: float(field 8) <= 1e-5  (inclusive; '<' is wrong, see D016)"

for dataset in observed shuffle01 shuffle02 shuffle03 shuffle04 shuffle05; do
  rel="$(awk -F'\t' -v d="$dataset" 'NR>1 && $1==d {print $10}' "$MANIFEST")"
  [[ -n "$rel" ]] || die "dataset ${dataset} has no fimo_output entry in the manifest"
  src="${MODULE_DIR}/${rel}"
  [[ -f "$src" ]] || die "primary FIMO output missing: ${src}"

  # The input must be the accepted artifact, by sha256 from the module manifest.
  want_sha="$(awk -v p="./${rel}" '$2==p{print $1}' "${MODULE_DIR}/checksums.sha256")"
  if [[ -n "$want_sha" ]]; then
    got_sha="$(sha256sum "$src" | cut -d' ' -f1)"
    [[ "$got_sha" == "$want_sha" ]] \
      || die "primary input checksum mismatch for ${dataset}
       expected ${want_sha}
       observed ${got_sha}"
  fi

  out="${OUT_ROOT}/${dataset}/fimo.tsv"
  log "${dataset}: ${rel} -> ${out#${MODULE_DIR}/}"

  python3 "${HERE}/filter_fimo_strict.py" \
    --input  "$src" \
    --output "$out" \
    --expect-rows     "${EXPECT_PRIMARY[$dataset]}" \
    --expect-strict   "${EXPECT_STRICT[$dataset]}" \
    --expect-boundary "${EXPECT_BOUNDARY[$dataset]}" 2>&1 | tee -a "$NCRNA_LOG_FILE"

  log "${dataset}: sha256 $(sha256sum "$out" | cut -d' ' -f1)"
  log "${dataset}: bytes  $(stat -c%s "$out")  mode $(stat -c%a "$out")"
done

log "done"
