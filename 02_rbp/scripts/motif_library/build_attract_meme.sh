#!/usr/bin/env bash
# =============================================================================
# 02_rbp / motif_library : build the canonical RBP MEME motif library
# =============================================================================
# Consumes the frozen, checksum-verified ATtRACT source from 01_sources and
# produces the motif library every downstream RBP stage scans with.
#
#   01_sources/frozen/attract/{ATtRACT_db.txt,pwm.txt}
#       -> filter to Homo_sapiens + CisBP-RNA entries with a usable Matrix_id
#       -> 02_rbp/results/motif_library/attract_hsa_cisbp.meme
#
# Acceptance (VAL04): exactly 93 unique RBPs across 151 motifs. The converter
# refuses to write a library that does not reproduce that universe.
#
# All scientific parameters come from config/parameters.yaml. No threshold, path
# or filter value is hard-coded here.
#
# Usage:  ./build_attract_meme.sh [--out-dir DIR]
# Env:    NCRNA_SYSTEM, NCRNA_RELEASE_ROOT (see config/paths.yaml)
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../../../config/lib/common.sh"

MODULE_DIR="${RELEASE_ROOT}/02_rbp"
PARAMS="${RELEASE_ROOT}/config/parameters.yaml"
ATTRACT_DIR="${SOURCES_DIR}/frozen/attract"

OUT_DIR="${MODULE_DIR}/results/motif_library"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    *) die "unknown argument: $1" ;;
  esac
done

mkdir -p "${MODULE_DIR}/logs" "${MODULE_DIR}/intermediate" "${OUT_DIR}"
export NCRNA_LOG_FILE="${MODULE_DIR}/logs/build_attract_meme.$(date +%Y%m%dT%H%M%S).log"

# ---------------------------------------------------------------- parameters
# Minimal YAML reader for the rbp.motif_library block. Values are scalars only.
param() {
  local key="$1"
  python3 - "$PARAMS" "$key" <<'PY'
import re, sys
path, key = sys.argv[1], sys.argv[2]
text = open(path).read()
block = re.search(r"^  motif_library:\n((?:    .*\n|\n)*)", text, re.M)
if not block:
    sys.exit("motif_library block not found in parameters.yaml")
m = re.search(rf"^    {re.escape(key)}:\s*(.+?)\s*(?:#.*)?$", block.group(1), re.M)
if not m:
    sys.exit(f"parameter not found: {key}")
print(m.group(1).strip().strip('"').strip("'"))
PY
}

ORGANISM_COLUMN="$(param organism_column)"
ORGANISM_VALUE="$(param organism_value)"
DATABASE_COLUMN="$(param database_column)"
DATABASE_VALUE="$(param database_value)"
MATRIX_ID_COLUMN="$(param matrix_id_column)"
GENE_NAME_COLUMN="$(param gene_name_column)"
MEME_LABEL="$(param meme_label)"
EXPECT_RBPS="$(param expected_unique_rbps)"
EXPECT_MOTIFS="$(param expected_motifs)"

# ---------------------------------------------------------------- inputs
ATTRACT_DB="${ATTRACT_DIR}/ATtRACT_db.txt"
PWM="${ATTRACT_DIR}/pwm.txt"

for f in "$ATTRACT_DB" "$PWM"; do
  [[ -f "$f" ]] || die "required frozen input missing: ${f}
     Run 01_sources/scripts/fetch_upstream_releases.sh and verify_sources.sh first."
done

log "system=${NCRNA_SYSTEM} release_root=${RELEASE_ROOT}"
log "attract_db=${ATTRACT_DB}"
log "pwm=${PWM}"
log "filters: ${ORGANISM_COLUMN}=${ORGANISM_VALUE}, ${DATABASE_COLUMN}=${DATABASE_VALUE}"
log "acceptance: ${EXPECT_RBPS} unique RBPs across ${EXPECT_MOTIFS} motifs"

# ---------------------------------------------------------------- build
OUT_MEME="${OUT_DIR}/attract_hsa_cisbp.meme"
OUT_MANIFEST="${OUT_DIR}/motif_library_manifest.tsv"

python3 "${HERE}/attract_to_meme.py" \
  --attract-db "$ATTRACT_DB" \
  --pwm "$PWM" \
  --out-meme "$OUT_MEME" \
  --out-manifest "$OUT_MANIFEST" \
  --label "$MEME_LABEL" \
  --organism-column "$ORGANISM_COLUMN" \
  --organism-value "$ORGANISM_VALUE" \
  --database-column "$DATABASE_COLUMN" \
  --database-value "$DATABASE_VALUE" \
  --matrix-id-column "$MATRIX_ID_COLUMN" \
  --gene-name-column "$GENE_NAME_COLUMN" \
  --expect-unique-rbps "$EXPECT_RBPS" \
  --expect-motifs "$EXPECT_MOTIFS" 2>&1 | tee -a "$NCRNA_LOG_FILE"

# ---------------------------------------------------------------- report
log "motif library : ${OUT_MEME}"
log "sha256        : $(sha256sum "$OUT_MEME" | cut -d' ' -f1)"
log "motifs        : $(grep -c '^MOTIF ' "$OUT_MEME")"
log "unique RBPs   : $(grep '^MOTIF ' "$OUT_MEME" | awk '{print $3}' | cut -d'|' -f1 | sort -u | wc -l)"
log "done"
