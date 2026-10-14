#!/usr/bin/env bash
# =============================================================================
# 03_transcript_context step 8 (effective order 15)
#   HERV LTRs x transcript-resolved features, OPPOSITE strand, with bp overlaps
# =============================================================================
# Refactor of scripts/context/07_overlap_ltrs_vs_features.sh.
#
# Historical 07 did two things. It rebuilt three feature BEDs inline (the UTR
# 5'/3' split, the 3'-UTR subset and the last exon per transcript) -- that logic
# is now step 4, whose outputs reproduce 07's assets byte for byte, so the shared
# step-4 assets are consumed here instead. And it ran, exactly as reproduced here:
#
#   bedtools sort -g <chrom sizes> -i <LTR BED> | gzip                     -> ltrs.sorted
#   bedtools intersect -sorted -g <chrom sizes> -S -wo -a <LTRs sorted> -b <feature sorted>
#
# -S keeps only pairs on OPPOSITE strands. -wo writes the LTR record, the feature
# record and the overlap length in bp. One row per overlapping (LTR, feature
# interval): an LTR overlapping one exon in several isoforms gives several rows.
# Nothing is collapsed here; the SPARCS candidate builder (07_sparcs) collapses
# per LTR x transcript x feature class. An opposite-strand overlap is a genomic
# configuration, not evidence of antisense transcription or dsRNA.
# Historical 07 wrote no same-strand or unstranded LTR table.
#
# Features (accepted, genome-sorted):
#   tx_exons   step-3 features/sorted/gencode.v48.exons.sorted.bed.gz
#   cds        step-3 features/sorted/gencode.v48.cds.sorted.bed.gz
#   utr_5p3p   step-4 features/sorted/gencode.v48.utr_5p3p.sorted.bed.gz
#   utr3       step-4 features/sorted/gencode.v48.utr3.sorted.bed.gz
#   last_exon  step-4 features/sorted/gencode.v48.last_exon_per_tx.sorted.bed.gz
# No gene-level merged exons/introns (D022) and no TSS (D023).
#
# Outputs in <results>/overlap_ltrs/ (historical names and 13-column header):
#   ltrs_vs_gencode_{tx_exons,cds,utr_5p3p,utr3,last_exon}.antisense.wo.tsv
# and <results>/herv_ltrs/ltrs.sorted.bed.gz.
#
# Usage: overlap_ltrs_vs_features.sh --results-root DIR --ltr-bed FILE --genome FILE
# =============================================================================
set -euo pipefail
export LC_ALL=C
usage() { echo "usage: $0 --results-root DIR --ltr-bed FILE --genome FILE" >&2; exit 1; }
RES=""; LTR=""; GENOME=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --results-root) RES="$2"; shift 2 ;;
    --ltr-bed)      LTR="$2"; shift 2 ;;
    --genome)       GENOME="$2"; shift 2 ;;
    *) usage ;;
  esac
done
[[ -n "$RES" && -n "$LTR" && -n "$GENOME" ]] || usage
RES="$(cd "$RES" && pwd)"; F="${RES}/gencode_features"

BEDTOOLS="${NCRNA_BEDTOOLS:-bedtools}"
echo "[env] which bedtools: $(command -v "$BEDTOOLS")"
echo "[env] bedtools --version: $("$BEDTOOLS" --version)"
[[ "$("$BEDTOOLS" --version)" == "bedtools v2.31.1" ]] || { echo "[error] release pins bedtools v2.31.1" >&2; exit 1; }
[[ "$(md5sum < "$GENOME" | cut -c1-32)" == "f94face851cb2129ab423f54db04fa0c" ]] || { echo "[error] chrom sizes md5" >&2; exit 1; }
[[ "$(md5sum < "$LTR" | cut -c1-32)" == "ef318953ea35e2821431f943c5e20007" ]] \
  || { echo "[error] LTR BED is not the accepted SRC07 release file" >&2; exit 1; }
bad_strand="$(awk -F'\t' '$6!="+" && $6!="-"' "$LTR" | wc -l)"
[[ "$bad_strand" == 0 ]] || { echo "[error] ${bad_strand} LTRs without strand + or -" >&2; exit 1; }
echo "[ok] SRC07 LTR BED: $(wc -l < "$LTR") LTRs, all with strand + or -"

declare -A B=(
  [tx_exons]="features/sorted/gencode.v48.exons.sorted.bed.gz"
  [cds]="features/sorted/gencode.v48.cds.sorted.bed.gz"
  [utr_5p3p]="features/sorted/gencode.v48.utr_5p3p.sorted.bed.gz"
  [utr3]="features/sorted/gencode.v48.utr3.sorted.bed.gz"
  [last_exon]="features/sorted/gencode.v48.last_exon_per_tx.sorted.bed.gz"
)
ORDER=(tx_exons cds utr_5p3p utr3 last_exon)
( cd "$F"
  for t in "${ORDER[@]}"; do grep -E "  \./${B[$t]//./\\.}\$" checksums.sha256; done | sha256sum -c - ) \
  || { echo "[error] feature inputs do not verify against ${F}/checksums.sha256" >&2; exit 1; }

mkdir -p "${RES}/overlap_ltrs" "${RES}/herv_ltrs"
STAGE="$(mktemp -d "${RES}/overlap_ltrs/.stage_ltrs.XXXXXX")"
trap '[[ -d "$STAGE" ]] && rm -rf "$STAGE"' EXIT
"$BEDTOOLS" sort -g "$GENOME" -i "$LTR" | gzip -n -6 > "${STAGE}/ltrs.sorted.bed.gz"
A="${STAGE}/ltrs.sorted.bed.gz"
HDR=$'#a_chrom\ta_start\ta_end\ta_name\ta_score\ta_strand\tb_chrom\tb_start\tb_end\tb_name\tb_score\tb_strand\toverlap_bp'
for t in "${ORDER[@]}"; do
  out="${STAGE}/ltrs_vs_gencode_${t}.antisense.wo.tsv"
  { echo "$HDR"; "$BEDTOOLS" intersect -sorted -g "$GENOME" -S -wo -a <(gzip -dc "$A") -b <(gzip -dc "${F}/${B[$t]}"); } > "$out"
  echo "[ok] ltrs_vs_gencode_${t}: $(($(wc -l < "$out")-1)) rows"
done

for f in "${STAGE}"/ltrs_vs_gencode_*.tsv; do
  [[ "$(head -1 "$f")" == "$HDR" ]] || { echo "[error] header $f" >&2; exit 1; }
  chmod 644 "$f"; mv "$f" "${RES}/overlap_ltrs/"
done
chmod 644 "$A"; mv "$A" "${RES}/herv_ltrs/ltrs.sorted.bed.gz"
( cd "${RES}/overlap_ltrs" && sha256sum ltrs_vs_gencode_*.tsv > checksums.sha256 )
( cd "${RES}/herv_ltrs" && sha256sum ltrs.sorted.bed.gz > checksums.sha256 )
echo "[ok] published 5 antisense LTR overlap tables and the sorted LTR BED"
