#!/usr/bin/env bash
# =============================================================================
# 02_rbp / null_model : canonical QC of the dinucleotide-shuffled replicates
# =============================================================================
# Step 20 of the frozen execution plan. Recomputes QC from the Stage 5 FASTAs on
# disk and writes the canonical standalone report.
#
# Read-only with respect to Stage 5: it never regenerates or rewrites a replicate.
#
# Usage:  ./run_qc_dinuc_shuffles.sh [--out-tsv PATH]
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../../../config/lib/common.sh"

MODULE_DIR="${RELEASE_ROOT}/02_rbp"
PARAMS="${RELEASE_ROOT}/config/parameters.yaml"
IN_FASTA="${SOURCES_DIR}/frozen/herv/HERV_internal_sequences_v5.fasta"
FASTA_DIR="${MODULE_DIR}/results/null_model/fasta"
OUT_TSV="${MODULE_DIR}/results/null_model/dinuc_shuffle_qc.tsv"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-tsv) OUT_TSV="$2"; shift 2 ;;
    *) die "unknown argument: $1" ;;
  esac
done

mkdir -p "${MODULE_DIR}/logs"
export NCRNA_LOG_FILE="${MODULE_DIR}/logs/qc_dinuc_shuffles.$(date +%Y%m%dT%H%M%S).log"

N_SHUFFLES="$(python3 "${HERE}/read_config.py" --file "$PARAMS" --key n_dinuc_shuffles)"
EXP_RECORDS="$(python3 "${HERE}/read_config.py" --file "$PARAMS" --key dinuc_shuffle_expected_records)"
EXP_BASES="$(python3 "${HERE}/read_config.py" --file "$PARAMS" --key dinuc_shuffle_expected_bases)"
BASE_SEED="$(python3 "${HERE}/read_config.py" --file "$PARAMS" --key dinuc_shuffle_base_seed)"

# Build the replicate:seed list from the explicit config table, and cross-check it
# against base_seed + replicate exactly as the Stage 5 producer does.
SEEDS=""
for i in $(seq 1 "$N_SHUFFLES"); do
  rep="$(printf '%02d' "$i")"
  seed="$(python3 "${HERE}/read_config.py" --file "$PARAMS" --seed-for "$rep")"
  derived="$(( BASE_SEED + i ))"
  [[ "$seed" == "$derived" ]] \
    || die "seed mapping conflict for replicate ${rep}: config table ${seed}, base_seed+${i} ${derived}"
  SEEDS="${SEEDS}${SEEDS:+,}${rep}:${seed}"
done

# The source must still be the published file this null model was built from.
MANIFEST="${SOURCES_DIR}/manifests/upstream_releases.tsv"
EXPECTED_MD5="$(awk -F'\t' '$4 == "HERV_internal_sequences_v5.fasta" {print $5; exit}' "$MANIFEST")"
[[ -n "$EXPECTED_MD5" ]] || die "no manifest entry for HERV_internal_sequences_v5.fasta"
verify_md5 "$IN_FASTA" "$EXPECTED_MD5"

log "system=${NCRNA_SYSTEM} release_root=${RELEASE_ROOT}"
log "original  = ${IN_FASTA}"
log "fasta dir = ${FASTA_DIR}"
log "report    = ${OUT_TSV}"
log "seeds     = ${SEEDS}"

python3 "${HERE}/qc_dinuc_shuffles.py" \
  --original "$IN_FASTA" \
  --fasta-dir "$FASTA_DIR" \
  --out-tsv "$OUT_TSV" \
  --seeds "$SEEDS" \
  --expect-records "$EXP_RECORDS" \
  --expect-bases "$EXP_BASES" 2>&1 | tee -a "$NCRNA_LOG_FILE"

log "report sha256 $(sha256sum "$OUT_TSV" | cut -d' ' -f1)"
log "done"
