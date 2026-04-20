#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage:
  $(basename "$0") --db ATtRACT_db.txt --pwm pwm.txt --outdir OUTDIR [options]

Required:
  --db PATH        Path to ATtRACT_db.txt
  --pwm PATH       Path to pwm.txt
  --outdir PATH    Output directory

Options:
  --prefix STR     Prefix for output files (default: attract_hsa)
  --cisbp-only     Only build CisBP-RNA (Database==C) MEME
  --all-only       Only build ALL Homo_sapiens MEME
  -h, --help       Show this help

Outputs (default prefix: attract_hsa):
  OUTDIR/
    ATtRACT_hsa_all.tsv
    ATtRACT_hsa_cisbp.tsv
    attract_hsa_all.meme
    attract_hsa_cisbp.meme
    build_report.txt

Notes:
  - MEME files use DNA alphabet ACGT.
  - ATtRACT PWMs are A C G U; U is treated as T (4th column).
EOF
}

DB_TXT=""
PWM_TXT=""
OUT_DIR=""
PREFIX="attract_hsa"
MODE="both"  # both|cisbp|all

# --- parse args ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --db) DB_TXT="$2"; shift 2 ;;
    --pwm) PWM_TXT="$2"; shift 2 ;;
    --outdir) OUT_DIR="$2"; shift 2 ;;
    --prefix) PREFIX="$2"; shift 2 ;;
    --cisbp-only) MODE="cisbp"; shift ;;
    --all-only) MODE="all"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

# --- validate ---
[[ -z "${DB_TXT}" ]] && { echo "[ERROR] --db is required" >&2; usage; exit 1; }
[[ -z "${PWM_TXT}" ]] && { echo "[ERROR] --pwm is required" >&2; usage; exit 1; }
[[ -z "${OUT_DIR}" ]] && { echo "[ERROR] --outdir is required" >&2; usage; exit 1; }

[[ ! -f "${DB_TXT}" ]] && { echo "[ERROR] DB file not found: ${DB_TXT}" >&2; exit 1; }
[[ ! -f "${PWM_TXT}" ]] && { echo "[ERROR] PWM file not found: ${PWM_TXT}" >&2; exit 1; }

mkdir -p "${OUT_DIR}"

TSV_ALL="${OUT_DIR}/ATtRACT_hsa_all.tsv"
TSV_CISBP="${OUT_DIR}/ATtRACT_hsa_cisbp.tsv"
MEME_ALL="${OUT_DIR}/${PREFIX}_all.meme"
MEME_CISBP="${OUT_DIR}/${PREFIX}_cisbp.meme"
REPORT="${OUT_DIR}/build_report.txt"

# --- filter metadata (always write TSVs; cheap + useful provenance) ---
awk -F'\t' 'NR==1 || ($4=="Homo_sapiens" && $12!="N/A" && $12!="")' \
  "${DB_TXT}" > "${TSV_ALL}"

awk -F'\t' 'NR==1 || ($4=="Homo_sapiens" && $8=="C" && $12!="N/A" && $12!="")' \
  "${DB_TXT}" > "${TSV_CISBP}"

# --- build MEME(s) ---
if [[ "${MODE}" == "both" || "${MODE}" == "all" ]]; then
  python3 02_build_attract_meme.py \
    --attract-tsv "${TSV_ALL}" \
    --pwm-txt "${PWM_TXT}" \
    --out-meme "${MEME_ALL}" \
    --label "ATtRACT_HSA_ALL"
fi

if [[ "${MODE}" == "both" || "${MODE}" == "cisbp" ]]; then
  python3 02_build_attract_meme.py \
    --attract-tsv "${TSV_CISBP}" \
    --pwm-txt "${PWM_TXT}" \
    --out-meme "${MEME_CISBP}" \
    --label "ATtRACT_HSA_CISBP"
fi

# --- report ---
{
  echo "ATtRACT build report"
  echo "===================="
  echo "Input DB:  ${DB_TXT}"
  echo "Input PWM: ${PWM_TXT}"
  echo "Output dir: ${OUT_DIR}"
  echo "Prefix: ${PREFIX}"
  echo "Mode: ${MODE}"
  echo
  echo "Filtered TSVs:"
  echo "  - ${TSV_ALL}     (rows: $(($(wc -l < "${TSV_ALL}") - 1)) + header)"
  echo "  - ${TSV_CISBP}   (rows: $(($(wc -l < "${TSV_CISBP}") - 1)) + header)"
  echo
  echo "MEME outputs:"
  if [[ -f "${MEME_ALL}" ]]; then
    echo "  - ${MEME_ALL} (motifs: $(grep -c '^MOTIF ' "${MEME_ALL}"))"
  fi
  if [[ -f "${MEME_CISBP}" ]]; then
    echo "  - ${MEME_CISBP} (motifs: $(grep -c '^MOTIF ' "${MEME_CISBP}"))"
  fi
  echo
  echo "Notes:"
  echo "  - MEME files use DNA alphabet ACGT."
  echo "  - ATtRACT PWMs are A C G U; U is treated as T (4th column)."
} > "${REPORT}"

echo "Done. Report: ${REPORT}"
