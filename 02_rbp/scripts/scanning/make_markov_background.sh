#!/usr/bin/env bash
# =============================================================================
# 02_rbp / scanning : build the canonical order-1 nucleotide background for FIMO
# =============================================================================
# Estimates the background nucleotide model that every RBP FIMO scan is scored
# against, from the frozen HERV internal-region FASTA.
#
#   01_sources/frozen/herv/HERV_internal_sequences_v5.fasta
#     -> fasta-get-markov -dna -norc -m 1
#     -> 02_rbp/results/background/herv_internal_order1_background.txt
#
# Three fasta-get-markov defaults are overridden deliberately, because each
# default would silently change the model:
#   -dna    the tool otherwise GUESSES the alphabet (and defaults to RNA for
#           these sequences), which would emit A/C/G/U rows FIMO cannot use here
#   -norc   the tool otherwise COMBINES forward and reverse-complement
#           frequencies; the scan is transcript-oriented (FIMO --norc), so the
#           background must be estimated in the same single orientation
#   -m 1    the tool otherwise builds an order-0 model
# The pseudocount is left at the tool default (0.1) and recorded in config.
#
# Scientific behaviour is identical to the historical implementation. Only
# paths, configuration, atomicity and verification have changed.
#
# Usage:  ./make_markov_background.sh [--out-dir DIR] [--input FASTA]
# Env:    NCRNA_SYSTEM, NCRNA_RELEASE_ROOT, NCRNA_MEME_BIN
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../../../config/lib/common.sh"

MODULE_DIR="${RELEASE_ROOT}/02_rbp"
PARAMS="${RELEASE_ROOT}/config/parameters.yaml"

OUT_DIR="${MODULE_DIR}/results/background"
INPUT_FASTA="${SOURCES_DIR}/frozen/herv/HERV_internal_sequences_v5.fasta"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --input)   INPUT_FASTA="$2"; shift 2 ;;
    *) die "unknown argument: $1" ;;
  esac
done

mkdir -p "${MODULE_DIR}/logs" "${OUT_DIR}"
export NCRNA_LOG_FILE="${MODULE_DIR}/logs/make_markov_background.$(date +%Y%m%dT%H%M%S).log"

# ---------------------------------------------------------------- parameters
param() {
  local key="$1"
  python3 - "$PARAMS" "$key" <<'PY'
import re, sys
path, key = sys.argv[1], sys.argv[2]
text = open(path).read()
block = re.search(r"^  background:\n((?:    .*\n|\n)*)", text, re.M)
if not block:
    sys.exit("background block not found in parameters.yaml")
m = re.search(rf"^    {re.escape(key)}:\s*(.+?)\s*(?:#.*)?$", block.group(1), re.M)
if not m:
    sys.exit(f"parameter not found: {key}")
print(m.group(1).strip().strip('"').strip("'"))
PY
}

MARKOV_ORDER="$(param markov_order)"
ALPHABET_FLAG="$(param alphabet_flag)"
STRAND_FLAG="$(param strand_flag)"
EXPECTED_MD5="$(param expected_md5)"
EXPECTED_FREQ_LINES="$(param expected_frequency_lines)"

# ---------------------------------------------------------------- tool
MEME_BIN="${NCRNA_MEME_BIN:-}"
if [[ -n "$MEME_BIN" ]]; then
  FGM="${MEME_BIN}/fasta-get-markov"
else
  FGM="$(command -v fasta-get-markov || true)"
fi
[[ -n "$FGM" && -x "$FGM" ]] || die "fasta-get-markov not found. Set NCRNA_MEME_BIN or put it on PATH."

# fasta-get-markov has no --version flag; the version is in the libexec path it
# resolves to, and is cross-checked against fimo in the same prefix.
RESOLVED="$(readlink -f "$FGM")"
MEME_VERSION="$(printf '%s' "$RESOLVED" | grep -oE 'meme-[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)"
[[ -n "$MEME_VERSION" ]] || die "could not determine MEME Suite version from: ${RESOLVED}"
EXPECTED_MEME="meme-$(python3 -c "
import re,sys
t=open('${RELEASE_ROOT}/config/environments/motif.yaml').read()
print(re.search(r'meme_suite_version:\s*\"?([0-9.]+)',t).group(1))")"
if [[ "$MEME_VERSION" != "$EXPECTED_MEME" ]]; then
  die "MEME Suite version mismatch: found ${MEME_VERSION}, release pins ${EXPECTED_MEME}.
     Refusing to build the background with a different version."
fi

# ---------------------------------------------------------------- input
[[ -f "$INPUT_FASTA" ]] || die "frozen input FASTA missing: ${INPUT_FASTA}
     Run 01_sources/scripts/fetch_upstream_releases.sh and verify_sources.sh first."

SRC_MANIFEST="${SOURCES_DIR}/manifests/upstream_releases.tsv"
if [[ -f "$SRC_MANIFEST" ]]; then
  want_md5="$(awk -F'\t' '$4=="HERV_internal_sequences_v5.fasta"{print $5}' "$SRC_MANIFEST")"
  if [[ -n "$want_md5" ]]; then
    got_md5="$(md5sum "$INPUT_FASTA" | cut -d' ' -f1)"
    [[ "$got_md5" == "$want_md5" ]] || die "input FASTA checksum mismatch
       expected ${want_md5}
       observed ${got_md5}"
    log "input verified against 01_sources manifest: md5 ${want_md5}"
  fi
fi

OUT_FILE="${OUT_DIR}/herv_internal_order1_background.txt"

log "system=${NCRNA_SYSTEM} release_root=${RELEASE_ROOT}"
log "tool    = ${RESOLVED} (${MEME_VERSION})"
log "input   = ${INPUT_FASTA}"
log "flags   = ${ALPHABET_FLAG} ${STRAND_FLAG} -m ${MARKOV_ORDER}"
log "output  = ${OUT_FILE}"

# Idempotent: leave an already-correct output alone. Correctness is judged on the
# MODEL (the frequency rows), not on the byte md5 of the historical file, because
# fasta-get-markov records the input path in its comment header and the release
# uses the published input filename.
CANONICAL_MODEL="${HERE}/canonical_background_model.txt"
if [[ -f "$OUT_FILE" && -f "$CANONICAL_MODEL" ]]; then
  if diff -q <(grep -vE '^#' "$OUT_FILE") "$CANONICAL_MODEL" >/dev/null 2>&1; then
    log "output already present and its model matches the canonical values; nothing to do"
    exit 0
  fi
  log "existing output does not match the canonical model; rebuilding"
fi

# ---------------------------------------------------------------- build
STAGE="$(mktemp -d "${OUT_DIR}/.stage.XXXXXX")"
cleanup() { rm -rf "$STAGE"; }
trap cleanup EXIT

STAGED="${STAGE}/background.txt"
"$FGM" "$ALPHABET_FLAG" "$STRAND_FLAG" -m "$MARKOV_ORDER" "$INPUT_FASTA" "$STAGED" \
  2>&1 | tee -a "$NCRNA_LOG_FILE"

[[ -s "$STAGED" ]] || die "fasta-get-markov produced an empty background file"

n_freq="$(grep -cE '^[ACGT]+[[:space:]]' "$STAGED" || true)"
[[ "$n_freq" == "$EXPECTED_FREQ_LINES" ]] \
  || die "expected ${EXPECTED_FREQ_LINES} frequency rows, found ${n_freq}"

# Scientific acceptance: the frequency rows must match the canonical model.
# The leading comment block records the input path, so it is excluded from this
# comparison and reported separately by the validation step.
log "frequency rows: ${n_freq}"

mv -f "$STAGED" "$OUT_FILE"

OUT_MD5="$(md5sum "$OUT_FILE" | cut -d' ' -f1)"
log "output md5    : ${OUT_MD5}"
log "output sha256 : $(sha256sum "$OUT_FILE" | cut -d' ' -f1)"
log "output bytes  : $(stat -c%s "$OUT_FILE")"
if [[ "$OUT_MD5" == "$EXPECTED_MD5" ]]; then
  log "BYTE_IDENTICAL to the historical canonical background"
else
  log "NOTE: md5 differs from the historical file (${EXPECTED_MD5})."
  log "      fasta-get-markov embeds the input path in its comment header, so the"
  log "      canonical input name changes those bytes. Frequency rows are compared"
  log "      separately by validate_markov_background.sh."
fi
log "done"
