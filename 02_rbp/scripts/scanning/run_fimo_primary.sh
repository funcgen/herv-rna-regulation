#!/usr/bin/env bash
# =============================================================================
# 02_rbp / scanning : primary observed FIMO scan (Stage 3)
# =============================================================================
# Scans the HERV internal regions with the canonical RBP motif library, scored
# against the order-1 strand-specific background.
#
#   HERV_internal_sequences_v5.fasta        53,176 records
# + attract_hsa_cisbp.meme                  151 motifs / 93 RBPs
# + herv_internal_order1_background.txt     order-1, -dna -norc
#     -> FIMO 5.5.9 --thresh 1e-4 --max-stored-scores 1e6 --norc --no-pgc
#     -> results/fimo_primary/fimo.tsv      716,075 motif occurrences
#
# The flag set is reproduced exactly from the historical cmd.txt provenance file;
# nothing is added or removed. In particular:
#   --norc     scan one orientation only; sequences are transcript-oriented
#   --no-pgc   do not parse genomic coordinates from sequence names; FIMO's
#              parser does not understand this project's header format, so
#              coordinates are reconstructed later by the Stage 4 parser
#   --thresh 0.0001  PRIMARY threshold. Not to be confused with the strict
#              sensitivity filter (p <= 1e-5) applied much later, downstream.
#
# This stage stops at raw FIMO output. It does not parse, shuffle or summarise.
#
# Usage:  ./run_fimo_primary.sh [--out-dir DIR]
# Env:    NCRNA_SYSTEM, NCRNA_RELEASE_ROOT, NCRNA_MEME_BIN
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../../../config/lib/common.sh"

MODULE_DIR="${RELEASE_ROOT}/02_rbp"
PARAMS="${RELEASE_ROOT}/config/parameters.yaml"

OUT_DIR="${MODULE_DIR}/results/fimo_primary"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    *) die "unknown argument: $1" ;;
  esac
done

mkdir -p "${MODULE_DIR}/logs"
export NCRNA_LOG_FILE="${MODULE_DIR}/logs/run_fimo_primary.$(date +%Y%m%dT%H%M%S).log"

# ---------------------------------------------------------------- parameters
rbp_param() {
  local key="$1"
  python3 - "$PARAMS" "$key" <<'PY'
import re, sys
path, key = sys.argv[1], sys.argv[2]
text = open(path).read()
block = re.search(r"^rbp:\n((?:  .*\n|\n)*)", text, re.M)
if not block:
    sys.exit("rbp block not found in parameters.yaml")
m = re.search(rf"^  {re.escape(key)}:\s*(.+?)\s*(?:#.*)?$", block.group(1), re.M)
if not m:
    sys.exit(f"parameter not found: {key}")
print(m.group(1).strip().strip('"').strip("'"))
PY
}

THRESH="$(rbp_param fimo_primary_pvalue)"
MAX_STORED="$(rbp_param fimo_max_stored_scores)"
SCAN_RC="$(rbp_param fimo_scan_reverse_complement)"

# parameters.yaml carries the threshold as 1.0e-4; the historical FIMO invocation
# used the decimal form 0.0001. They are the same number, and 'g' formatting
# renders it as 0.0001, so the recorded command matches the historical one
# character for character.
THRESH_ARG="$(python3 -c "print(format(float('${THRESH}'), 'g'))")"
[[ "$THRESH_ARG" == "0.0001" ]] \
  || die "threshold rendered as '${THRESH_ARG}', expected '0.0001' (config value: ${THRESH})"

[[ "$SCAN_RC" == "false" ]] || die "config requests reverse-complement scanning; the canonical scan is --norc"

# ---------------------------------------------------------------- tool
MEME_BIN="${NCRNA_MEME_BIN:-}"
FIMO="${MEME_BIN:+${MEME_BIN}/}fimo"
command -v "$FIMO" >/dev/null 2>&1 || FIMO="$(command -v fimo || true)"
[[ -n "$FIMO" ]] || die "fimo not found. Set NCRNA_MEME_BIN or put it on PATH."

FIMO_VERSION="$("$FIMO" --version 2>&1 | head -1 | tr -d '[:space:]')"
EXPECTED_VERSION="$(python3 -c "
import re
t=open('${RELEASE_ROOT}/config/environments/motif.yaml').read()
print(re.search(r'meme_suite_version:\s*\"?([0-9.]+)',t).group(1))")"
[[ "$FIMO_VERSION" == "$EXPECTED_VERSION" ]] \
  || die "FIMO version mismatch: found '${FIMO_VERSION}', release pins '${EXPECTED_VERSION}'"

# ---------------------------------------------------------------- inputs
FASTA="${SOURCES_DIR}/frozen/herv/HERV_internal_sequences_v5.fasta"
MOTIFS="${MODULE_DIR}/results/motif_library/attract_hsa_cisbp.meme"
BGFILE="${MODULE_DIR}/results/background/herv_internal_order1_background.txt"

for f in "$FASTA" "$MOTIFS" "$BGFILE"; do
  [[ -f "$f" ]] || die "required input missing: ${f}"
done

# Verify the two upstream stage outputs and the frozen source before a 30-minute run.
SRC_MANIFEST="${SOURCES_DIR}/manifests/upstream_releases.tsv"
want_fasta_md5="$(awk -F'\t' '$4=="HERV_internal_sequences_v5.fasta"{print $5}' "$SRC_MANIFEST")"
got_fasta_md5="$(md5sum "$FASTA" | cut -d' ' -f1)"
[[ "$got_fasta_md5" == "$want_fasta_md5" ]] \
  || die "input FASTA checksum mismatch (expected ${want_fasta_md5}, got ${got_fasta_md5})"

EXPECTED_MOTIF_MD5="9ed27b363078830c6d4f29801d7be146"
got_motif_md5="$(md5sum "$MOTIFS" | cut -d' ' -f1)"
[[ "$got_motif_md5" == "$EXPECTED_MOTIF_MD5" ]] \
  || die "motif library checksum mismatch (expected ${EXPECTED_MOTIF_MD5}, got ${got_motif_md5})"

CANONICAL_BG_MODEL="${HERE}/canonical_background_model.txt"
diff -q <(grep -vE '^#' "$BGFILE") "$CANONICAL_BG_MODEL" >/dev/null \
  || die "background model does not match the canonical frequency values"

n_motifs="$(grep -c '^MOTIF ' "$MOTIFS")"
n_rbps="$(grep '^MOTIF ' "$MOTIFS" | awk '{print $3}' | cut -d'|' -f1 | sort -u | wc -l)"
[[ "$n_motifs" == "151" && "$n_rbps" == "93" ]] \
  || die "motif library universe is ${n_motifs} motifs / ${n_rbps} RBPs, expected 151 / 93"

log "system=${NCRNA_SYSTEM} release_root=${RELEASE_ROOT}"
log "fimo    = $(command -v "$FIMO") (version ${FIMO_VERSION})"
log "fasta   = ${FASTA} (md5 ${got_fasta_md5})"
log "motifs  = ${MOTIFS} (${n_motifs} motifs, ${n_rbps} RBPs)"
log "bgfile  = ${BGFILE} (model verified against canonical values)"
log "thresh  = ${THRESH_ARG}  max-stored-scores = ${MAX_STORED}  flags = --norc --no-pgc"

# ---------------------------------------------------------------- idempotency
EXPECTED_ROWS=716075
if [[ -f "${OUT_DIR}/fimo.tsv" ]]; then
  rows="$(awk 'NR>1 && !/^#/ && NF==10' "${OUT_DIR}/fimo.tsv" | wc -l)"
  if [[ "$rows" == "$EXPECTED_ROWS" ]]; then
    log "output already present with ${rows} data rows; nothing to do"
    exit 0
  fi
  log "existing output has ${rows} rows, expected ${EXPECTED_ROWS}; rebuilding"
fi

# ---------------------------------------------------------------- run
# FIMO writes a directory, so stage it and publish only on success.
STAGE="$(mktemp -d "$(dirname "$OUT_DIR")/.stage_fimo.XXXXXX")"
cleanup() { [[ -d "$STAGE" ]] && rm -rf "$STAGE"; }
trap cleanup EXIT

log "starting FIMO scan (historically ~30 minutes)"
set +e
"$FIMO" \
  --thresh "$THRESH_ARG" \
  --max-stored-scores "$MAX_STORED" \
  --bgfile "$BGFILE" \
  --oc "${STAGE}/out" \
  --no-pgc \
  --norc \
  "$MOTIFS" \
  "$FASTA" 2>&1 | tee -a "$NCRNA_LOG_FILE"
fimo_status="${PIPESTATUS[0]}"
set -e
[[ "$fimo_status" == "0" ]] || die "FIMO exited with status ${fimo_status}"

[[ -s "${STAGE}/out/fimo.tsv" ]] || die "FIMO produced no fimo.tsv"

rows="$(awk 'NR>1 && !/^#/ && NF==10' "${STAGE}/out/fimo.tsv" | wc -l)"
[[ "$rows" == "$EXPECTED_ROWS" ]] \
  || die "FIMO produced ${rows} data rows, expected ${EXPECTED_ROWS} (VAL05). Refusing to publish."
log "VAL05 ok: ${rows} motif occurrences"

# Record the command alongside the output, as the historical runs did.
printf '%s --thresh %s --max-stored-scores %s --bgfile %s --oc %s --no-pgc --norc %s %s\n' \
  "fimo" "$THRESH_ARG" "$MAX_STORED" "$BGFILE" "${OUT_DIR}" "$MOTIFS" "$FASTA" \
  > "${STAGE}/out/cmd.txt"

rm -rf "$OUT_DIR"
mkdir -p "$(dirname "$OUT_DIR")"
mv "${STAGE}/out" "$OUT_DIR"

log "published: ${OUT_DIR}"
log "fimo.tsv bytes  : $(stat -c%s "${OUT_DIR}/fimo.tsv")"
log "fimo.tsv md5    : $(md5sum "${OUT_DIR}/fimo.tsv" | cut -d' ' -f1)"
log "data-row sha256 : $(awk 'NR>1 && !/^#/ && NF==10' "${OUT_DIR}/fimo.tsv" | sha256sum | cut -d' ' -f1)"
log "done"
