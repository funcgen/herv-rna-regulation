#!/usr/bin/env bash
# =============================================================================
# 02_rbp / scanning : parse the primary observed FIMO scan into genomic coordinates
# =============================================================================
#   results/fimo_primary/fimo.tsv                 716,075 raw calls
#     -> decode sequence names, reconstruct GRCh38 coordinates
#     -> results/scanning/primary_fimo_parsed.tsv  716,075 records, 15 columns
#     -> results/scanning/primary_fimo_parsed.bed  716,075 records, BED6
#
# Coordinate reconstruction only. No filtering, deduplication, re-scoring or
# integration with any annotation.
#
# Usage:  ./run_parse_fimo_primary.sh [--out-dir DIR]
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../../../config/lib/common.sh"

MODULE_DIR="${RELEASE_ROOT}/02_rbp"
IN_TSV="${MODULE_DIR}/results/fimo_primary/fimo.tsv"
OUT_DIR="${MODULE_DIR}/results/scanning"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    *) die "unknown argument: $1" ;;
  esac
done

mkdir -p "${MODULE_DIR}/logs" "${OUT_DIR}"
export NCRNA_LOG_FILE="${MODULE_DIR}/logs/parse_fimo_primary.$(date +%Y%m%dT%H%M%S).log"

# Stage 3 regression reference: the parser refuses to run on anything else.
EXPECTED_INPUT_SHA="3b2c524b9264af583906ac7d7abf5fc411659d993e3d2b0367dbf27a8c59cb00"
EXPECTED_ROWS=716075

[[ -f "$IN_TSV" ]] || die "primary FIMO output missing: ${IN_TSV}
     Run run_fimo_primary.sh (Stage 3) first."

OUT_TSV="${OUT_DIR}/primary_fimo_parsed.tsv"
OUT_BED="${OUT_DIR}/primary_fimo_parsed.bed"

log "system=${NCRNA_SYSTEM} release_root=${RELEASE_ROOT}"
log "input   = ${IN_TSV}"
log "outputs = ${OUT_TSV}"
log "          ${OUT_BED}"

# Idempotent: leave a correct existing result alone.
if [[ -f "$OUT_TSV" && -f "$OUT_BED" ]]; then
  rows="$(( $(wc -l < "$OUT_TSV") - 1 ))"
  bed_rows="$(wc -l < "$OUT_BED")"
  if [[ "$rows" == "$EXPECTED_ROWS" && "$bed_rows" == "$EXPECTED_ROWS" ]]; then
    log "parsed outputs already present with ${EXPECTED_ROWS} records; nothing to do"
    exit 0
  fi
  log "existing outputs have ${rows}/${bed_rows} records; rebuilding"
fi

python3 "${HERE}/parse_fimo_to_bed.py" \
  --input "$IN_TSV" \
  --out-tsv "$OUT_TSV" \
  --out-bed "$OUT_BED" \
  --expect-input-sha256 "$EXPECTED_INPUT_SHA" \
  --expect-rows "$EXPECTED_ROWS" 2>&1 | tee -a "$NCRNA_LOG_FILE"

log "tsv sha256 : $(sha256sum "$OUT_TSV" | cut -d' ' -f1)"
log "bed sha256 : $(sha256sum "$OUT_BED" | cut -d' ' -f1)"
log "tsv bytes  : $(stat -c%s "$OUT_TSV")"
log "bed bytes  : $(stat -c%s "$OUT_BED")"
log "done"
