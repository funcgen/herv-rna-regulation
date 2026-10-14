#!/usr/bin/env bash
# =============================================================================
# 01_sources : source-level validation (VAL01, VAL02, VAL03)
# =============================================================================
# VAL01  every published upstream file present on this system matches its md5
# VAL02  HERV internal FASTA contains 53,176 sequence records
# VAL03  HERV_ltr.bed contains 560,594 records
#
# Exits non-zero on any failure. This is the gate every downstream module runs
# against before it may execute.
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../../config/lib/common.sh"

MANIFEST="${SOURCES_DIR}/manifests/upstream_releases.tsv"
mkdir -p "${LOG_DIR}" "${SOURCES_DIR}/validation"
export NCRNA_LOG_FILE="${LOG_DIR}/verify_sources.$(date +%Y%m%dT%H%M%S).log"
REPORT="${SOURCES_DIR}/validation/validation_report.tsv"
printf 'validation_id\tsource_id\tcheck\texpected\tobserved\tstatus\n' > "$REPORT"

fail=0
emit() { printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5" "$6" >> "$REPORT"; }

# ---- VAL01 -----------------------------------------------------------------
while IFS=$'\t' read -r sid rec doi fname md5 size alias target; do
  [[ "$sid" == "source_id" || -z "$sid" ]] && continue
  case "$sid" in
    SRC03|SRC04|SRC05|SRC06) d="${SOURCES_DIR}/frozen/herv" ;;
    SRC07|SRC08|SRC09|SRC10) d="${SOURCES_DIR}/frozen/ltr" ;;
    SRC12a|SRC12b)           d="${SOURCES_DIR}/frozen/attract" ;;
    SRC13)                   d="${SOURCES_DIR}/derived" ;;
  esac
  f="${d}/${fname}"
  if [[ "$target" != "$NCRNA_SYSTEM" ]]; then
    emit VAL01 "$sid" md5 "$md5" "n/a" NOT_REQUIRED_ON_THIS_SYSTEM; continue
  fi
  if [[ ! -f "$f" ]]; then
    emit VAL01 "$sid" md5 "$md5" MISSING FAILED; fail=1; log "FAIL ${sid} missing: $f"; continue
  fi
  got="$(md5sum "$f" | cut -d' ' -f1)"
  if [[ "$got" == "$md5" ]]; then emit VAL01 "$sid" md5 "$md5" "$got" PASS; log "PASS ${sid} ${fname}"
  else emit VAL01 "$sid" md5 "$md5" "$got" FAILED; fail=1; log "FAIL ${sid} checksum mismatch"; fi
done < "$MANIFEST"

# ---- VAL02 -----------------------------------------------------------------
FA="${SOURCES_DIR}/frozen/herv/HERV_internal_sequences_v5.fasta"
if [[ -f "$FA" ]]; then
  n="$(grep -c '^>' "$FA")"
  if [[ "$n" == "53176" ]]; then emit VAL02 SRC03 record_count 53176 "$n" PASS; log "PASS VAL02 53,176 records"
  else emit VAL02 SRC03 record_count 53176 "$n" FAILED; fail=1; log "FAIL VAL02 got $n"; fi
else emit VAL02 SRC03 record_count 53176 ABSENT NOT_REQUIRED_ON_THIS_SYSTEM; fi

# ---- VAL03 -----------------------------------------------------------------
BED="${SOURCES_DIR}/frozen/ltr/HERV_ltr.bed"
if [[ -f "$BED" ]]; then
  n="$(wc -l < "$BED")"
  if [[ "$n" == "560594" ]]; then emit VAL03 SRC07 record_count 560594 "$n" PASS; log "PASS VAL03 560,594 records"
  else emit VAL03 SRC07 record_count 560594 "$n" FAILED; fail=1; log "FAIL VAL03 got $n"; fi
else emit VAL03 SRC07 record_count 560594 ABSENT NOT_REQUIRED_ON_THIS_SYSTEM; fi

log "validation report: ${REPORT}"
[[ "$fail" == "0" ]] || die "source validation FAILED - downstream modules must not run"
log "ALL SOURCE VALIDATIONS PASSED"
