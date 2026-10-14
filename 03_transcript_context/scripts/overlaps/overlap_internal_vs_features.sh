#!/usr/bin/env bash
# =============================================================================
# 03_transcript_context step 5 (effective order 12)
#   HERV internal regions x gene features, any-strand and same-strand
# =============================================================================
# Refactor of scripts/context/03_overlap_internal_regions.sh. Every command and
# flag is reproduced exactly:
#
#   bedtools sort -g <chrom sizes> -i <HERV internal BED>  | gzip   -> internal_regions.sorted
#   any-strand : bedtools intersect      -wa -wb -a <HERV sorted> -b <feature sorted>
#   same-strand: bedtools intersect -s   -wa -wb -a <HERV sorted> -b <feature sorted>
#
# No -sorted, -g, -f, -F, -r, -u, -v or -wo. Every overlapping (HERV, feature)
# pair is written once per overlapping feature interval, in bedtools' output
# order (A-file order, then B order); nothing is deduplicated or re-sorted.
# Each file starts with the historical 12-column header.
#
# Feature classes (accepted step-3 sorted assets):
#   genes         genes/sorted/gencode.v48.genes.sorted.bed.gz
#   gene_exons    features/sorted/gencode.v48.gene_exons_merged.sorted.bed.gz  (D022 per-gene)
#   gene_introns  features/sorted/gencode.v48.gene_introns.sorted.bed.gz       (D022 per-gene)
#   cds           features/sorted/gencode.v48.cds.sorted.bed.gz
#   utr           features/sorted/gencode.v48.utr.sorted.bed.gz               (generic UTR, as historically)
#
# D022: the gene-level exon/intron inputs are the corrected per-gene assets;
# the historical genome-wide files in validation/historical_reproduction/ are
# forensic only and are refused. D023: no TSS input. The step-4 UTR split and
# terminal exons are not inputs of this step.
#
# Outputs (plain TSV, historical names) in <results>/overlap/:
#   internal_vs_{genes,gene_exons,gene_introns,cds,utr}.tsv
#   internal_vs_{genes,gene_exons,gene_introns,cds,utr}.same_strand.tsv
# plus <results>/herv_internal/internal_regions.sorted.bed.gz (also consumed by
# effective order 13). Staged, then renamed into place.
#
# Usage: overlap_internal_vs_features.sh --results-root DIR --herv-bed FILE --genome FILE
# =============================================================================
set -euo pipefail
export LC_ALL=C

usage() { echo "usage: $0 --results-root DIR --herv-bed FILE --genome FILE" >&2; exit 1; }
RES=""; HERV=""; GENOME=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --results-root) RES="$2"; shift 2 ;;
    --herv-bed)     HERV="$2"; shift 2 ;;
    --genome)       GENOME="$2"; shift 2 ;;
    *) echo "unknown argument: $1" >&2; usage ;;
  esac
done
[[ -n "$RES" && -n "$HERV" && -n "$GENOME" ]] || usage
RES="$(cd "$RES" && pwd)"
F="${RES}/gencode_features"

BEDTOOLS="${NCRNA_BEDTOOLS:-bedtools}"
command -v "$BEDTOOLS" >/dev/null 2>&1 || { echo "[error] bedtools not found: $BEDTOOLS" >&2; exit 1; }
echo "[env] which bedtools: $(command -v "$BEDTOOLS")"
echo "[env] bedtools --version: $("$BEDTOOLS" --version)"
[[ "$("$BEDTOOLS" --version)" == "bedtools v2.31.1" ]] || { echo "[error] release pins bedtools v2.31.1" >&2; exit 1; }

# ---------------------------------------------------------------- accepted inputs
[[ "$(md5sum < "$HERV" | cut -c1-32)" == "dc0c77396229b1616e39a2cb13b294c0" ]] \
  || { echo "[error] HERV internal BED is not the accepted SRC04 release file" >&2; exit 1; }
[[ "$(wc -l < "$HERV")" == 53176 ]] || { echo "[error] HERV BED row count" >&2; exit 1; }
[[ "$(md5sum < "$GENOME" | cut -c1-32)" == "f94face851cb2129ab423f54db04fa0c" ]] \
  || { echo "[error] chrom sizes md5" >&2; exit 1; }

declare -A B=(
  [genes]="genes/sorted/gencode.v48.genes.sorted.bed.gz"
  [gene_exons]="features/sorted/gencode.v48.gene_exons_merged.sorted.bed.gz"
  [gene_introns]="features/sorted/gencode.v48.gene_introns.sorted.bed.gz"
  [cds]="features/sorted/gencode.v48.cds.sorted.bed.gz"
  [utr]="features/sorted/gencode.v48.utr.sorted.bed.gz"
)
ORDER=(genes gene_exons gene_introns cds utr)
( cd "$F"
  for t in "${ORDER[@]}"; do
    case "${B[$t]}" in *historical_reproduction*|*tss*|*utr_5p3p*|*last_exon*) echo "[error] refused input ${B[$t]}" >&2; exit 1 ;; esac
    grep -E "  \./${B[$t]//./\\.}\$" checksums.sha256
  done | sha256sum -c --quiet - ) || { echo "[error] feature inputs do not verify against ${F}/checksums.sha256" >&2; exit 1; }
echo "[ok] HERV BED (SRC04, 53,176 loci), chrom sizes and 5 step-3 feature BEDs verified"

# ---------------------------------------------------------------- stage
mkdir -p "${RES}/overlap" "${RES}/herv_internal"
STAGE="$(mktemp -d "${RES}/overlap/.stage_overlap.XXXXXX")"
trap '[[ -d "$STAGE" ]] && rm -rf "$STAGE"' EXIT

"$BEDTOOLS" sort -g "$GENOME" -i "$HERV" | gzip -n -6 > "${STAGE}/internal_regions.sorted.bed.gz"
A="${STAGE}/internal_regions.sorted.bed.gz"
HDR=$'#a_chrom\ta_start\ta_end\ta_name\ta_score\ta_strand\tb_chrom\tb_start\tb_end\tb_name\tb_score\tb_strand'

for t in "${ORDER[@]}"; do
  b="${F}/${B[$t]}"
  { echo "$HDR"; "$BEDTOOLS" intersect -wa -wb -a <(gzip -dc "$A") -b <(gzip -dc "$b"); } \
    > "${STAGE}/internal_vs_${t}.tsv"
  { echo "$HDR"; "$BEDTOOLS" intersect -s -wa -wb -a <(gzip -dc "$A") -b <(gzip -dc "$b"); } \
    > "${STAGE}/internal_vs_${t}.same_strand.tsv"
  echo "[ok] internal_vs_${t}: $(($(wc -l < "${STAGE}/internal_vs_${t}.tsv")-1)) any-strand, $(($(wc -l < "${STAGE}/internal_vs_${t}.same_strand.tsv")-1)) same-strand rows"
done

# ---------------------------------------------------------------- publish
for f in "${STAGE}"/internal_vs_*.tsv; do
  [[ "$(head -1 "$f")" == "$HDR" ]] || { echo "[error] header $f" >&2; exit 1; }
  chmod 644 "$f"; mv "$f" "${RES}/overlap/"
done
chmod 644 "$A"; mv "$A" "${RES}/herv_internal/internal_regions.sorted.bed.gz"
echo "[ok] published 10 overlap tables to ${RES}/overlap/ and the sorted HERV BED to ${RES}/herv_internal/"
