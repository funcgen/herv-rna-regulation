#!/usr/bin/env bash
# =============================================================================
# 02_rbp / sensitivity_strict : strict p <= 1e-5 re-analysis (Stage 12)
# =============================================================================
# Runs the ACCEPTED Stage 8, 9 and 10 analysis scripts, unmodified, against the
# strict dataset branch.
#
#                      one implementation of each calculation
#                                     |
#                    +----------------+----------------+
#                    |                                 |
#            primary manifest                  strict manifest
#            p <= 1e-4 (as run)                p <= 1e-5 (post-filtered)
#
# Historically this was three whole-file script copies (08/09/10_*_1e-5.R)
# differing only in path constants -- a direct diff of the recovered pairs shows
# 2, 2 and 1 changed lines, all paths, zero scientific difference. The release
# does not fork the code: it passes different paths to the same scripts.
#
# THIS SCRIPT CONTAINS NO SCIENTIFIC LOGIC. It resolves paths, builds the strict
# dataset manifest and invokes the same three analysis scripts.
#
# FIMO is not run. The strict inputs are the accepted Stage 11 post-filter.
#
# Usage:  ./run_strict_analysis.sh [--profile strict|primary]
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../../../config/lib/common.sh"

MODULE_DIR="${RELEASE_ROOT}/02_rbp"
PROFILES="${RELEASE_ROOT}/config/analysis_profiles.yaml"
PROFILE="strict"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --profile) PROFILE="$2"; shift 2 ;;
    *) die "unknown argument: $1" ;;
  esac
done

mkdir -p "${MODULE_DIR}/logs"
export NCRNA_LOG_FILE="${MODULE_DIR}/logs/strict_analysis.$(date +%Y%m%dT%H%M%S).log"

prof() { python3 "${HERE}/read_profile.py" --file "$PROFILES" --profile "$PROFILE" --key "$1"; }

LABEL="$(prof analysis_label)"
THRESHOLD="$(prof threshold)"
SOURCE="$(prof source)"
MANIFEST_REL="$(prof dataset_manifest)"
OBS_FIMO="$(prof observed_fimo)"
SHUF_TMPL="$(prof shuffle_fimo_template)"
GLOBAL_OUT="$(prof global_output_root)"
SUBFAM_OUT="$(prof subfamily_output_root)"
FOLLOWUP_OUT="$(prof followup_output_root)"
EXPECT_OBS="$(prof expect_observed)"
EXPECT_SHUF="$(prof expect_shuffles)"

FASTA="${RELEASE_ROOT}/$(python3 "${HERE}/read_profile.py" --file "$PROFILES" --shared fasta)"
MOTIFS="${RELEASE_ROOT}/$(python3 "${HERE}/read_profile.py" --file "$PROFILES" --shared motif_library)"
MANIFEST="${RELEASE_ROOT}/${MANIFEST_REL}"

log "profile   = ${LABEL}  threshold ${THRESHOLD}  source ${SOURCE}"
log "manifest  = ${MANIFEST}"
log "fasta     = ${FASTA}"
log "motifs    = ${MOTIFS}"

RSCRIPT="${NCRNA_RSCRIPT:-}"
if [[ -z "$RSCRIPT" ]]; then
  RSCRIPT="$(python3 -c "
import re
t = open('${RELEASE_ROOT}/config/environments/r_stats.yaml').read()
m = re.search(r'^\s*rscript:\s*\"?([^\"\n]+)\"?\s*$', t, re.M)
print(m.group(1).strip() if m else '')")"
fi
[[ -x "$RSCRIPT" ]] || die "Rscript not found: ${RSCRIPT}"

# ------------------------------------------------- 1. dataset manifest
mkdir -p "$(dirname "$MANIFEST")"
log "building the ${LABEL} dataset manifest"
python3 "${MODULE_DIR}/scripts/null_model/build_fimo_dataset_manifest.py" \
  --module-dir "$MODULE_DIR" \
  --config "${RELEASE_ROOT}/config/parameters.yaml" \
  --out "$MANIFEST" \
  --checksums-out "${MODULE_DIR}/validation/fimo_${LABEL}_datarow_checksums.tsv" \
  --observed-fimo "$OBS_FIMO" \
  --shuffle-template "$SHUF_TMPL" \
  --analysis-label "$LABEL" \
  --generated-in-stage-observed "Stage 11" \
  --generated-in-stage-shuffle "Stage 11" \
  --expect-observed "$EXPECT_OBS" \
  --expect-shuffles "$EXPECT_SHUF" 2>&1 | tee -a "$NCRNA_LOG_FILE"

# ------------------------------------------------- 2. Stage 8 logic
log "global observed-versus-null summary (Stage 8 code, ${LABEL} inputs)"
"$RSCRIPT" "${MODULE_DIR}/scripts/enrichment/summarise_global_null.R" \
  --module-dir "$MODULE_DIR" \
  --manifest   "$MANIFEST" \
  --fasta      "$FASTA" \
  --out-dir    "${RELEASE_ROOT}/${GLOBAL_OUT}" 2>&1 | tee -a "$NCRNA_LOG_FILE"

# ------------------------------------------------- 3. Stage 9 logic
log "subfamily composition-controlled enrichment (Stage 9 code, ${LABEL} inputs)"
"$RSCRIPT" "${MODULE_DIR}/scripts/enrichment/subfamily_enrichment.R" \
  --module-dir "$MODULE_DIR" \
  --manifest   "$MANIFEST" \
  --fasta      "$FASTA" \
  --motifs     "$MOTIFS" \
  --out-dir    "${RELEASE_ROOT}/${SUBFAM_OUT}" 2>&1 | tee -a "$NCRNA_LOG_FILE"

# ------------------------------------------------- 4. Stage 10 logic
log "HERVH/HERVK profiles, GO ORA and pre-ranked GSEA (Stage 10 code, ${LABEL} inputs)"
"$RSCRIPT" "${MODULE_DIR}/scripts/enrichment/hervh_hervk_go_gsea.R" \
  --stage9-dir "${RELEASE_ROOT}/${SUBFAM_OUT}" \
  --out-dir    "${RELEASE_ROOT}/${FOLLOWUP_OUT}" 2>&1 | tee -a "$NCRNA_LOG_FILE"

# Record what this branch was, so no reader has to infer it from a directory name.
cat > "${RELEASE_ROOT}/${GLOBAL_OUT}/../analysis_profile.tsv" <<TSV
key	value
analysis_profile	${LABEL}
threshold	${THRESHOLD}
source	${SOURCE}
dataset_manifest	${MANIFEST_REL}
note	Post-filter of the primary FIMO output (D016). FIMO was NOT rerun at this threshold.
TSV

log "done"
