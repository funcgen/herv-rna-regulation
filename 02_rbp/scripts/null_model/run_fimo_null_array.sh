#!/usr/bin/env bash
# =============================================================================
# 02_rbp / null_model : FIMO scan of one dinucleotide-shuffled replicate (Stage 7)
# =============================================================================
# The null-model twin of the Stage 3 observed scan. Everything is held fixed
# except the sequence set:
#
#   attract_hsa_cisbp.meme                 151 motifs / 93 RBPs   (Stage 1)
# + herv_internal_order1_background.txt    ONE observed-derived   (Stage 2)
# + HERV_internal_sequences_dinuc_shuffle_NN.fasta                (Stage 5)
#     -> FIMO 5.5.9 --thresh 0.0001 --max-stored-scores 1000000 --norc --no-pgc
#     -> results/null_model/fimo/shuffleNN/fimo.tsv
#
# ONE BACKGROUND, SIX DATASETS. The historical launcher
# (run_07_fimo_rbp_dinuc_null_array.sh, recovered and read directly) supplied the
# same observed-derived bg_internal_m1.sense.txt to the observed scan and to all
# five shuffles. That is a scientific property of the null: the motif model and
# the background are held fixed while only sequence ARRANGEMENT varies, under
# exact dinucleotide preservation. Fitting a background per replicate would
# absorb the very composition signal the null exists to measure.
#
# The observed dataset is NOT rescanned here: Stage 3 already produced it and its
# output is immutable. It is represented in the six-dataset manifest instead.
#
# This stage stops at raw FIMO output. No parsing, no merging, no statistics.
#
# Usage:  ./run_fimo_null_array.sh --replicate NN [--out-dir DIR]
# Env:    NCRNA_SYSTEM, NCRNA_RELEASE_ROOT, NCRNA_MEME_BIN
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../../../config/lib/common.sh"

MODULE_DIR="${RELEASE_ROOT}/02_rbp"
PARAMS="${RELEASE_ROOT}/config/parameters.yaml"
REPLICATE=""
OUT_ROOT="${MODULE_DIR}/results/null_model/fimo"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --replicate) REPLICATE="$2"; shift 2 ;;
    --out-dir)   OUT_ROOT="$2";  shift 2 ;;
    *) die "unknown argument: $1" ;;
  esac
done

[[ -n "$REPLICATE" ]] || die "--replicate NN is required (01..05)"
REPLICATE="$(printf '%02d' "$((10#$REPLICATE))")"

mkdir -p "${MODULE_DIR}/logs"
export NCRNA_LOG_FILE="${MODULE_DIR}/logs/fimo_null_${REPLICATE}.$(date +%Y%m%dT%H%M%S).log"

cfg() { python3 "${HERE}/read_config.py" --file "$PARAMS" "$@"; }

# ---------------------------------------------------------------- parameters
# Read from the SAME config keys the Stage 3 observed scan uses. There is one
# canonical parameter set; this stage does not define a second one.
THRESH="$(cfg --key fimo_primary_pvalue)"
MAX_STORED="$(cfg --key fimo_max_stored_scores)"
SCAN_RC="$(cfg --key fimo_scan_reverse_complement)"
SEED="$(cfg --seed-for "$REPLICATE")"
EXPECTED_ROWS="$(cfg --hits-for "$REPLICATE")"

THRESH_ARG="$(python3 -c "print(format(float('${THRESH}'), 'g'))")"
[[ "$THRESH_ARG" == "0.0001" ]] \
  || die "threshold rendered as '${THRESH_ARG}', expected '0.0001' (config value: ${THRESH})"
[[ "$SCAN_RC" == "false" ]] || die "config requests reverse-complement scanning; the canonical scan is --norc"
[[ "$MAX_STORED" == "1000000" ]] || die "max-stored-scores is '${MAX_STORED}', expected 1000000"
[[ -n "$SEED" && -n "$EXPECTED_ROWS" ]] || die "no seed / expected hit count configured for replicate ${REPLICATE}"

# The seed must agree with the two independent statements of the mapping, exactly
# as Stage 5 required. An array renumbering must not be able to change identity.
BASE_SEED="$(cfg --key dinuc_shuffle_base_seed)"
derived="$(( BASE_SEED + 10#$REPLICATE ))"
[[ "$SEED" == "$derived" ]] \
  || die "seed mapping conflict for replicate ${REPLICATE}: config ${SEED}, base_seed+${REPLICATE} ${derived}"

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
FASTA="${MODULE_DIR}/results/null_model/fasta/HERV_internal_sequences_dinuc_shuffle_${REPLICATE}.fasta"
MOTIFS="${MODULE_DIR}/results/motif_library/attract_hsa_cisbp.meme"
BGFILE="${MODULE_DIR}/results/background/herv_internal_order1_background.txt"

for f in "$FASTA" "$MOTIFS" "$BGFILE"; do
  [[ -f "$f" ]] || die "required input missing: ${f}"
done

EXPECTED_MOTIF_MD5="9ed27b363078830c6d4f29801d7be146"
got_motif_md5="$(md5sum "$MOTIFS" | cut -d' ' -f1)"
[[ "$got_motif_md5" == "$EXPECTED_MOTIF_MD5" ]] \
  || die "motif library checksum mismatch (expected ${EXPECTED_MOTIF_MD5}, got ${got_motif_md5})"

# The background is verified against the canonical frequency values, not a byte
# md5, because the file's comment header records the input path (see Stage 2).
CANONICAL_BG_MODEL="${HERE}/../scanning/canonical_background_model.txt"
diff -q <(grep -vE '^#' "$BGFILE") "$CANONICAL_BG_MODEL" >/dev/null \
  || die "background model does not match the canonical frequency values"

n_motifs="$(grep -c '^MOTIF ' "$MOTIFS")"
n_rbps="$(grep '^MOTIF ' "$MOTIFS" | awk '{print $3}' | cut -d'|' -f1 | sort -u | wc -l)"
[[ "$n_motifs" == "151" && "$n_rbps" == "93" ]] \
  || die "motif library universe is ${n_motifs} motifs / ${n_rbps} RBPs, expected 151 / 93"

# The shuffled FASTA must be the accepted Stage 5 artifact, by sha256 from the
# module manifest. A shuffle regenerated with a different seed would otherwise
# scan silently.
MODULE_SUMS="${MODULE_DIR}/checksums.sha256"
rel="./results/null_model/fasta/$(basename "$FASTA")"
want_fasta_sha="$(awk -v p="$rel" '$2==p{print $1}' "$MODULE_SUMS")"
[[ -n "$want_fasta_sha" ]] || die "shuffled FASTA ${rel} has no entry in ${MODULE_SUMS}"
got_fasta_sha="$(sha256sum "$FASTA" | cut -d' ' -f1)"
[[ "$got_fasta_sha" == "$want_fasta_sha" ]] \
  || die "shuffled FASTA checksum mismatch for replicate ${REPLICATE}
       expected ${want_fasta_sha}
       observed ${got_fasta_sha}"

n_records="$(grep -c '^>' "$FASTA")"
[[ "$n_records" == "$(cfg --key dinuc_shuffle_expected_records)" ]] \
  || die "shuffled FASTA has ${n_records} records, expected $(cfg --key dinuc_shuffle_expected_records)"

OUT_DIR="${OUT_ROOT}/shuffle${REPLICATE}"

log "system=${NCRNA_SYSTEM} release_root=${RELEASE_ROOT}"
log "replicate = ${REPLICATE}  seed = ${SEED}  expected hits = ${EXPECTED_ROWS}"
log "fimo    = $(command -v "$FIMO") (version ${FIMO_VERSION})"
log "fasta   = ${FASTA} (sha256 ${got_fasta_sha}, ${n_records} records)"
log "motifs  = ${MOTIFS} (${n_motifs} motifs, ${n_rbps} RBPs, md5 ${got_motif_md5})"
log "bgfile  = ${BGFILE} (observed-derived; the SAME file used by Stage 3 and by every replicate)"
log "thresh  = ${THRESH_ARG}  max-stored-scores = ${MAX_STORED}  flags = --norc --no-pgc"

# ---------------------------------------------------------------- idempotency
if [[ -f "${OUT_DIR}/fimo.tsv" ]]; then
  rows="$(awk 'NR>1 && !/^#/ && NF==10' "${OUT_DIR}/fimo.tsv" | wc -l)"
  if [[ "$rows" == "$EXPECTED_ROWS" ]]; then
    log "output already present with ${rows} data rows; nothing to do"
    exit 0
  fi
  log "existing output has ${rows} rows, expected ${EXPECTED_ROWS}; rebuilding"
fi

# ---------------------------------------------------------------- run
mkdir -p "$OUT_ROOT"
STAGE="$(mktemp -d "${OUT_ROOT}/.stage_fimo_${REPLICATE}.XXXXXX")"
cleanup() { [[ -d "$STAGE" ]] && rm -rf "$STAGE"; }
trap cleanup EXIT

log "starting FIMO scan for shuffle${REPLICATE} (Stage 3 observed took ~27 min)"
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

# Structural gate before publication.
hdr="$(head -1 "${STAGE}/out/fimo.tsv")"
[[ "$(awk -F'\t' 'NR==1{print NF}' <<<"$hdr")" == "10" ]] \
  || die "fimo.tsv header does not have 10 columns"
grep -q '^# FIMO ' "${STAGE}/out/fimo.tsv" \
  || log "note: no '# FIMO' trailer line found"

rows="$(awk 'NR>1 && !/^#/ && NF==10' "${STAGE}/out/fimo.tsv" | wc -l)"
[[ "$rows" == "$EXPECTED_ROWS" ]] \
  || die "FIMO produced ${rows} data rows, expected ${EXPECTED_ROWS} (VAL09). Refusing to publish."
log "VAL09 ok: ${rows} motif occurrences for shuffle${REPLICATE}"

printf '%s --thresh %s --max-stored-scores %s --bgfile %s --oc %s --no-pgc --norc %s %s\n' \
  "fimo" "$THRESH_ARG" "$MAX_STORED" "$BGFILE" "${OUT_DIR}" "$MOTIFS" "$FASTA" \
  > "${STAGE}/out/cmd.txt"

rm -rf "$OUT_DIR"
mkdir -p "$OUT_ROOT"
mv "${STAGE}/out" "$OUT_DIR"

loci="$(awk 'NR>1 && !/^#/ && NF==10 {print $3}' "${OUT_DIR}/fimo.tsv" | sort -u | wc -l)"
log "published: ${OUT_DIR}"
log "positive loci   : ${loci}"
log "fimo.tsv bytes  : $(stat -c%s "${OUT_DIR}/fimo.tsv")"
log "fimo.tsv md5    : $(md5sum "${OUT_DIR}/fimo.tsv" | cut -d' ' -f1)"
log "data-row sha256 : $(awk 'NR>1 && !/^#/ && NF==10' "${OUT_DIR}/fimo.tsv" | sha256sum | cut -d' ' -f1)"
log "done"
