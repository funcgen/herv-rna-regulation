#!/usr/bin/env bash
# =============================================================================
# 03_transcript_context: transcript TSS BED + its genome-ordered copy
#   D023 dependency repair; prerequisite of effective order 13 (nearest TSS)
# =============================================================================
#   accepted step-1 transcripts/gencode.v48.transcripts.bed.gz
#     -> build_transcript_tss.py
#     -> tss/gencode.v48.tss.transcripts.bed.gz                (GTF order, as historically)
#     -> zcat | bedtools sort -g <chrom sizes> -i - | gzip -n   (the accepted step-3 command)
#     -> tss/sorted/gencode.v48.tss.transcripts.sorted.bed.gz  (consumed by order 13)
#
# Historically the raw file was input/annotations/tss_gencode.transcripts.bed,
# made by scripts/dsrna/bidirectional/01_extract_gencode_tss_bed.py, and step 3
# (02_sort_context_assets.sh) sorted it into context/assets/tss/sorted/. Only the
# TSS asset is sorted here; the accepted step-3 outputs are not touched.
#
# Usage: build_tss.sh --results-root DIR --genome FILE
# =============================================================================
set -euo pipefail
export LC_ALL=C
usage() { echo "usage: $0 --results-root DIR --genome FILE" >&2; exit 1; }
ROOT=""; GENOME=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --results-root) ROOT="$2"; shift 2 ;;
    --genome)       GENOME="$2"; shift 2 ;;
    *) usage ;;
  esac
done
[[ -n "$ROOT" && -n "$GENOME" ]] || usage
ROOT="$(cd "$ROOT" && pwd)"; HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BEDTOOLS="${NCRNA_BEDTOOLS:-bedtools}"; PYTHON="${NCRNA_PYTHON:-python3}"
[[ "$("$BEDTOOLS" --version)" == "bedtools v2.31.1" ]] || { echo "[error] release pins bedtools v2.31.1" >&2; exit 1; }
echo "[env] $(command -v "$BEDTOOLS") $("$BEDTOOLS" --version); $("$PYTHON" -V 2>&1)"
[[ "$(md5sum < "$GENOME" | cut -c1-32)" == "f94face851cb2129ab423f54db04fa0c" ]] || { echo "[error] chrom sizes md5" >&2; exit 1; }
( cd "$ROOT" && grep -E '  \./transcripts/gencode\.v48\.transcripts\.bed\.gz$' checksums.sha256 | sha256sum -c --quiet - ) \
  || { echo "[error] accepted transcript BED does not verify" >&2; exit 1; }
echo "[ok] accepted step-1 transcript BED verified"

mkdir -p "${ROOT}/tss/sorted"
"$PYTHON" "${HERE}/build_transcript_tss.py" \
  --transcripts "${ROOT}/transcripts/gencode.v48.transcripts.bed.gz" --out "${ROOT}/tss/gencode.v48.tss.transcripts.bed.gz"

OUT="${ROOT}/tss/sorted/gencode.v48.tss.transcripts.sorted.bed.gz"
TMP="$(mktemp "${ROOT}/tss/sorted/.stage_tss_sorted.XXXXXX")"
trap 'rm -f "$TMP"' EXIT
zcat "${ROOT}/tss/gencode.v48.tss.transcripts.bed.gz" | "$BEDTOOLS" sort -g "$GENOME" -i - | gzip -n -6 > "$TMP"
[[ "$(zcat "$TMP" | wc -l)" == "$(zcat "${ROOT}/tss/gencode.v48.tss.transcripts.bed.gz" | wc -l)" ]] || { echo "[error] sort changed row count" >&2; exit 1; }
chmod 644 "$TMP"; mv "$TMP" "$OUT"; trap - EXIT
for f in "${ROOT}/tss/gencode.v48.tss.transcripts.bed.gz" "$OUT"; do
  echo "[info] $(basename "$f")  rows $(zcat "$f" | wc -l)  gz $(sha256sum < "$f" | cut -c1-64)  content $(zcat "$f" | sha256sum | cut -c1-64)"
done
