#!/usr/bin/env bash
# =============================================================================
# 03_transcript_context step 3 (effective order 10)
#   sort the canonical step-1 and step-2 assets into GRCh38 genome order
# =============================================================================
# Refactor of scripts/context/02_sort_context_assets.sh. ORDERING ONLY: no row is
# added, removed, deduplicated or modified, and no coordinate or name changes.
#
# Sort command, reproduced from history:
#     zcat <asset> | bedtools sort -g <chrom sizes> -i -
# Keys: chromosome rank in the chrom-sizes file, then start (ascending). End,
# name, score and strand are NOT keys. Within equal (chrom, start), bedtools
# 2.31.1 orders records by its internal (unstable) sort: deterministic for the
# same input bytes and binary, but not a declarative rule and not idempotent.
# Using the historical command reproduces the historical sorted content of every
# asset D022 did not change. Downstream steps must not depend on tie order.
#
# Inputs: the seven accepted canonical assets, verified against
# results/gencode_features/checksums.sha256. The step-2 files are the D022
# per-gene corrected assets; the historical step-2 reproduction under
# validation/historical_reproduction/step2/ is forensic and never read here.
#
# The historical script also sorted input/annotations/tss_gencode.transcripts.bed.
# That file is not a step-1 or step-2 output; it came from the EXCLUDE-classified
# dsrna/bidirectional branch and has no canonical producer (erratum D023, open).
# It is deliberately not sorted here.
#
# Outputs: <class>/sorted/<asset>.sorted.bed.gz beside the unsorted sources,
# which are left untouched. gzip -n: no name or mtime, so output is
# byte-reproducible. Staged on shared storage, validated by
# validate_sorted_assets.py, then renamed into place.
#
# Usage: sort_assets.sh --results-root DIR --genome FILE
# =============================================================================
set -euo pipefail
export LC_ALL=C

usage() { echo "usage: $0 --results-root DIR --genome FILE" >&2; exit 1; }
ROOT=""; GENOME=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --results-root) ROOT="$2"; shift 2 ;;
    --genome)       GENOME="$2"; shift 2 ;;
    *) echo "unknown argument: $1" >&2; usage ;;
  esac
done
[[ -n "$ROOT" && -n "$GENOME" ]] || usage
ROOT="$(cd "$ROOT" && pwd)"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

case "$ROOT" in
  *historical_reproduction*) echo "[error] refusing a historical-reproduction path as input root" >&2; exit 1 ;;
esac

# ---------------------------------------------------------------- tools
BEDTOOLS="${NCRNA_BEDTOOLS:-bedtools}"
command -v "$BEDTOOLS" >/dev/null 2>&1 || { echo "[error] bedtools not found: $BEDTOOLS" >&2; exit 1; }
BT_VERSION="$("$BEDTOOLS" --version 2>&1 | head -1 | awk '{print $2}')"
[[ "$BT_VERSION" == "v2.31.1" ]] || { echo "[error] bedtools ${BT_VERSION}, release pins v2.31.1" >&2; exit 1; }
PYTHON="${NCRNA_PYTHON:-python3}"
echo "[ok] bedtools $(command -v "$BEDTOOLS") ${BT_VERSION}; $("$PYTHON" -V 2>&1)"

# ---------------------------------------------------------------- chrom order
GENOME_MD5="f94face851cb2129ab423f54db04fa0c"
got="$(md5sum < "$GENOME" | cut -c1-32)"
[[ "$got" == "$GENOME_MD5" ]] || { echo "[error] chrom sizes md5 ${got}, expected ${GENOME_MD5}" >&2; exit 1; }
echo "[ok] chrom sizes ${GENOME} md5 ${got}, $(wc -l < "$GENOME") sequences"

# ---------------------------------------------------------------- assets
# class/asset (paths relative to ROOT)
ASSETS=(
  genes/gencode.v48.genes
  transcripts/gencode.v48.transcripts
  features/gencode.v48.exons
  features/gencode.v48.cds
  features/gencode.v48.utr
  features/gencode.v48.gene_exons_merged
  features/gencode.v48.gene_introns
)

# accepted inputs: every asset must be listed in, and verify against, the manifest
( cd "$ROOT"
  for a in "${ASSETS[@]}"; do
    grep -q "  \./${a}\.bed\.gz\$" checksums.sha256 || { echo "[error] ${a}.bed.gz not in checksums.sha256" >&2; exit 1; }
  done
  grep -E "  \./($(IFS='|'; echo "${ASSETS[*]}" | sed 's/\./\\./g'))\.bed\.gz\$" checksums.sha256 | sha256sum -c --quiet -
) || { echo "[error] input verification failed" >&2; exit 1; }
echo "[ok] ${#ASSETS[@]} accepted inputs verified against ${ROOT}/checksums.sha256"

# ---------------------------------------------------------------- stage
STAGE="$(mktemp -d "${ROOT}/.stage_sort.XXXXXX")"
cleanup() { [[ -d "$STAGE" ]] && rm -rf "$STAGE"; }
trap cleanup EXIT

for a in "${ASSETS[@]}"; do
  cls="$(dirname "$a")"; name="$(basename "$a")"
  mkdir -p "${STAGE}/${cls}/sorted"
  zcat "${ROOT}/${a}.bed.gz" \
    | "$BEDTOOLS" sort -g "$GENOME" -i - \
    | gzip -n -6 > "${STAGE}/${cls}/sorted/${name}.sorted.bed.gz"
  echo "[ok] sorted ${a}  ($(zcat "${STAGE}/${cls}/sorted/${name}.sorted.bed.gz" | wc -l) rows)"
done

# ---------------------------------------------------------------- validation gate
"$PYTHON" "${HERE}/validate_sorted_assets.py" \
  --results-root "$ROOT" --sorted-root "$STAGE" --genome "$GENOME" \
  || { echo "[error] validation of staged outputs failed; nothing published" >&2; exit 1; }

# ---------------------------------------------------------------- publish
for a in "${ASSETS[@]}"; do
  cls="$(dirname "$a")"; name="$(basename "$a")"
  mkdir -p "${ROOT}/${cls}/sorted"
  chmod 644 "${STAGE}/${cls}/sorted/${name}.sorted.bed.gz"
  mv "${STAGE}/${cls}/sorted/${name}.sorted.bed.gz" "${ROOT}/${cls}/sorted/${name}.sorted.bed.gz"
done
echo "[ok] published ${#ASSETS[@]} sorted assets under ${ROOT}/*/sorted/"
for a in "${ASSETS[@]}"; do
  f="${ROOT}/$(dirname "$a")/sorted/$(basename "$a").sorted.bed.gz"
  echo "[info] $(basename "$f")  gz $(sha256sum < "$f" | cut -c1-64)  content $(zcat "$f" | sha256sum | cut -c1-64)"
done
