#!/usr/bin/env bash
# =============================================================================
# 03_transcript_context step 6 (effective order 13)
#   nearest gene and nearest transcript TSS for every HERV internal region
# =============================================================================
# Refactor of scripts/context/04_closest_gene_tss.sh. The four commands and flags
# are reproduced exactly:
#
#   bedtools closest -g <chrom sizes>    -d -a <HERV sorted> -b <genes sorted>
#   bedtools closest -g <chrom sizes> -s -d -a <HERV sorted> -b <genes sorted>
#   bedtools closest -g <chrom sizes>    -d -a <HERV sorted> -b <TSS sorted>
#   bedtools closest -g <chrom sizes> -s -d -a <HERV sorted> -b <TSS sorted>
#
# No -t, so bedtools' default "-t all": every equally close feature is reported
# as its own row. -d appends the distance: 0 for overlap; otherwise the gap in
# bases + 1 (so book-ended features are at distance 1); unsigned. A HERV with no
# candidate on its sequence gets one row of ". -1 -1 . -1 ." and distance -1.
# Output order is bedtools' (A order; ties in B order). Each file starts with the
# historical 13-column header. Nothing is deduplicated or re-sorted.
#
# Inputs (accepted): step-5 results/herv_internal/internal_regions.sorted.bed.gz,
# step-3 genes/sorted/gencode.v48.genes.sorted.bed.gz, and the D023 transcript TSS
# tss/sorted/gencode.v48.tss.transcripts.sorted.bed.gz. Not D022-dependent.
#
# Outputs (plain TSV, historical names) in <results>/closest/:
#   internal_closest_gene.tsv  internal_closest_gene_same_strand.tsv
#   internal_closest_tss.tsv   internal_closest_tss_same_strand.tsv
#
# Usage: closest_gene_tss.sh --results-root DIR --genome FILE
# =============================================================================
set -euo pipefail
export LC_ALL=C
usage() { echo "usage: $0 --results-root DIR --genome FILE" >&2; exit 1; }
RES=""; GENOME=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --results-root) RES="$2"; shift 2 ;;
    --genome)       GENOME="$2"; shift 2 ;;
    *) usage ;;
  esac
done
[[ -n "$RES" && -n "$GENOME" ]] || usage
RES="$(cd "$RES" && pwd)"

BEDTOOLS="${NCRNA_BEDTOOLS:-bedtools}"
echo "[env] which bedtools: $(command -v "$BEDTOOLS")"
echo "[env] bedtools --version: $("$BEDTOOLS" --version)"
[[ "$("$BEDTOOLS" --version)" == "bedtools v2.31.1" ]] || { echo "[error] release pins bedtools v2.31.1" >&2; exit 1; }
[[ "$(md5sum < "$GENOME" | cut -c1-32)" == "f94face851cb2129ab423f54db04fa0c" ]] || { echo "[error] chrom sizes md5" >&2; exit 1; }

A="${RES}/herv_internal/internal_regions.sorted.bed.gz"
GENES="${RES}/gencode_features/genes/sorted/gencode.v48.genes.sorted.bed.gz"
TSS="${RES}/gencode_features/tss/sorted/gencode.v48.tss.transcripts.sorted.bed.gz"
( cd "${RES}/herv_internal" && sha256sum -c --quiet checksums.sha256 ) || { echo "[error] HERV input does not verify" >&2; exit 1; }
# (sha256sum -c without --quiet lists each verified file. Do not route this through
# `tee /dev/stderr`: when stderr is the SLURM log, tee reopens and truncates it.)
( cd "${RES}/gencode_features" && grep -E '  \./(genes/sorted/gencode\.v48\.genes|tss/sorted/gencode\.v48\.tss\.transcripts)\.sorted\.bed\.gz$' checksums.sha256 \
    | sha256sum -c - ) || { echo "[error] gene/TSS inputs do not verify" >&2; exit 1; }
[[ "$(zcat "$A" | wc -l)" == 53176 && "$(zcat "$GENES" | wc -l)" == 78894 && "$(zcat "$TSS" | wc -l)" == 387954 ]] \
  || { echo "[error] input universe sizes" >&2; exit 1; }
echo "[ok] accepted inputs verified: 53,176 HERVs, 78,894 genes, 387,954 transcript TSSs"

mkdir -p "${RES}/closest"
STAGE="$(mktemp -d "${RES}/closest/.stage_closest.XXXXXX")"
trap '[[ -d "$STAGE" ]] && rm -rf "$STAGE"' EXIT
HDR=$'#a_chrom\ta_start\ta_end\ta_name\ta_score\ta_strand\tb_chrom\tb_start\tb_end\tb_name\tb_score\tb_strand\tdistance_bp'

run() {  # out b_file [flags...]
  local out="$1" b="$2"; shift 2
  { echo "$HDR"; "$BEDTOOLS" closest -g "$GENOME" "$@" -d -a <(gzip -dc "$A") -b <(gzip -dc "$b"); } > "${STAGE}/${out}"
  echo "[ok] ${out}: $(($(wc -l < "${STAGE}/${out}")-1)) rows"
}
run internal_closest_gene.tsv             "$GENES"
run internal_closest_gene_same_strand.tsv "$GENES" -s
run internal_closest_tss.tsv              "$TSS"
run internal_closest_tss_same_strand.tsv  "$TSS"   -s

for f in "${STAGE}"/internal_closest_*.tsv; do
  [[ "$(head -1 "$f")" == "$HDR" ]] || { echo "[error] header $f" >&2; exit 1; }
  chmod 644 "$f"; mv "$f" "${RES}/closest/"
done
( cd "${RES}/closest" && sha256sum internal_closest_*.tsv > checksums.sha256 )
echo "[ok] published 4 closest tables to ${RES}/closest/"
