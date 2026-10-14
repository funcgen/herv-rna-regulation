#!/usr/bin/env bash
# =============================================================================
# 03_transcript_context step 7 (effective order 14)
#   HERV protein-domain predictions x gene/transcript features, same strand, with
#   base-pair overlap lengths
# =============================================================================
# Merges scripts/context/06_overlap_domains_vs_features.sh and
# scripts/context/05_overlap_domains_vs_exons.sh. 05's only output,
# domains_vs_gene_exons.same_strand.wo.tsv, was byte-identical to 06's
# domains_vs_gencode_gene_exons_merged.same_strand.wo.tsv (same relationship,
# re-run 2026-04-15), so it is produced once, under 06's name.
#
# Commands, reproduced exactly:
#   bedtools sort -g <chrom sizes> -i <domain BED> | gzip                  -> domains.sorted
#   bedtools intersect -sorted -g <chrom sizes> -s -wo -a <domains sorted> -b <feature sorted>
# -wo writes the A record, the B record and the overlap length in bp; every
# overlapping (domain, feature interval) pair is one row -- a domain overlapping
# the same exon in several isoforms, or several exons, yields several rows.
# Nothing is merged, reduced or deduplicated here. No coverage filter.
#
# Feature classes (accepted, genome-sorted):
#   tx_exons           step-3 features/sorted/gencode.v48.exons.sorted.bed.gz        (transcript-resolved)
#   cds                step-3 features/sorted/gencode.v48.cds.sorted.bed.gz
#   utr                step-3 features/sorted/gencode.v48.utr.sorted.bed.gz          (generic, unsplit)
#   utr_5p3p           step-4 features/sorted/gencode.v48.utr_5p3p.sorted.bed.gz     (historically rebuilt inline
#                      by 06 with the identical rule; now the shared step-4 asset)
#   gene_exons_merged  step-3 features/sorted/gencode.v48.gene_exons_merged.sorted.bed.gz (D022 per-gene)
# Terminal exons and TSSs are not inputs of this step.
#
# Outputs in <results>/overlap/ (plain TSV, historical 13-column header):
#   domains_vs_gencode_{tx_exons,cds,utr,utr_5p3p,gene_exons_merged}.same_strand.wo.tsv
# and <results>/herv_domains/domains.sorted.bed.gz.
#
# Usage: overlap_domains_vs_features.sh --results-root DIR --domains-bed FILE --genome FILE
# =============================================================================
set -euo pipefail
export LC_ALL=C
usage() { echo "usage: $0 --results-root DIR --domains-bed FILE --genome FILE" >&2; exit 1; }
RES=""; DOM=""; GENOME=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --results-root) RES="$2"; shift 2 ;;
    --domains-bed)  DOM="$2"; shift 2 ;;
    --genome)       GENOME="$2"; shift 2 ;;
    *) usage ;;
  esac
done
[[ -n "$RES" && -n "$DOM" && -n "$GENOME" ]] || usage
RES="$(cd "$RES" && pwd)"; F="${RES}/gencode_features"

BEDTOOLS="${NCRNA_BEDTOOLS:-bedtools}"
echo "[env] which bedtools: $(command -v "$BEDTOOLS")"
echo "[env] bedtools --version: $("$BEDTOOLS" --version)"
[[ "$("$BEDTOOLS" --version)" == "bedtools v2.31.1" ]] || { echo "[error] release pins bedtools v2.31.1" >&2; exit 1; }
[[ "$(md5sum < "$GENOME" | cut -c1-32)" == "f94face851cb2129ab423f54db04fa0c" ]] || { echo "[error] chrom sizes md5" >&2; exit 1; }
[[ "$(md5sum < "$DOM" | cut -c1-32)" == "f02b8bb89f19a0a2a3ac45efcd900898" ]] \
  || { echo "[error] domain BED is not the accepted SRC05 release file" >&2; exit 1; }
echo "[ok] SRC05 domain BED: $(wc -l < "$DOM") predictions"

declare -A B=(
  [tx_exons]="features/sorted/gencode.v48.exons.sorted.bed.gz"
  [cds]="features/sorted/gencode.v48.cds.sorted.bed.gz"
  [utr]="features/sorted/gencode.v48.utr.sorted.bed.gz"
  [utr_5p3p]="features/sorted/gencode.v48.utr_5p3p.sorted.bed.gz"
  [gene_exons_merged]="features/sorted/gencode.v48.gene_exons_merged.sorted.bed.gz"
)
ORDER=(tx_exons cds utr utr_5p3p gene_exons_merged)
# (no `tee /dev/stderr`: it truncates the SLURM log)
( cd "$F"
  for t in "${ORDER[@]}"; do grep -E "  \./${B[$t]//./\\.}\$" checksums.sha256; done | sha256sum -c - ) \
  || { echo "[error] feature inputs do not verify against ${F}/checksums.sha256" >&2; exit 1; }

mkdir -p "${RES}/overlap" "${RES}/herv_domains"
STAGE="$(mktemp -d "${RES}/overlap/.stage_domains.XXXXXX")"
trap '[[ -d "$STAGE" ]] && rm -rf "$STAGE"' EXIT
"$BEDTOOLS" sort -g "$GENOME" -i "$DOM" | gzip -n -6 > "${STAGE}/domains.sorted.bed.gz"
A="${STAGE}/domains.sorted.bed.gz"
HDR=$'#a_chrom\ta_start\ta_end\ta_name\ta_score\ta_strand\tb_chrom\tb_start\tb_end\tb_name\tb_score\tb_strand\toverlap_bp'
for t in "${ORDER[@]}"; do
  out="${STAGE}/domains_vs_gencode_${t}.same_strand.wo.tsv"
  { echo "$HDR"; "$BEDTOOLS" intersect -sorted -g "$GENOME" -s -wo -a <(gzip -dc "$A") -b <(gzip -dc "${F}/${B[$t]}"); } > "$out"
  echo "[ok] domains_vs_gencode_${t}: $(($(wc -l < "$out")-1)) rows"
done

for f in "${STAGE}"/domains_vs_gencode_*.tsv; do
  [[ "$(head -1 "$f")" == "$HDR" ]] || { echo "[error] header $f" >&2; exit 1; }
  chmod 644 "$f"; mv "$f" "${RES}/overlap/"
done
chmod 644 "$A"; mv "$A" "${RES}/herv_domains/domains.sorted.bed.gz"
( cd "${RES}/overlap" && sha256sum internal_vs_*.tsv domains_vs_*.tsv > checksums.sha256 )
( cd "${RES}/herv_domains" && sha256sum domains.sorted.bed.gz > checksums.sha256 )
echo "[ok] published 5 domain overlap tables and the sorted domain BED"
