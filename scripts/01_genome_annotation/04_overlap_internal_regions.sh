#!/usr/bin/env bash
set -euo pipefail

# ---- SET INTERNAL BED FILE ----
INTERNAL_BED="bed/internals_merged.bed"
# -------------------------------

GENOME="input/genome/grch38.primary.genome"
LOG="context/logs/03_overlap_internal_regions.log"
: > "$LOG"

mkdir -p context/assets context/results/overlap context/logs

if [[ ! -s "$GENOME" ]]; then
  echo "[ERROR] Genome file not found or empty: $GENOME" | tee -a "$LOG"
  echo "Create it from a .fai: cut -f1,2 GRCh38.fa.fai > $GENOME" | tee -a "$LOG"
  exit 1
fi

[[ -f "$INTERNAL_BED" ]] || { echo "[ERROR] Missing $INTERNAL_BED" | tee -a "$LOG"; exit 1; }

# Ensure internal regions are sorted using the SAME genome order as other assets
INTERNAL_SORTED="context/assets/internal_regions.sorted.bed.gz"
echo "[INFO] bedtools sort internal regions: $INTERNAL_BED -> $INTERNAL_SORTED" | tee -a "$LOG"
bedtools sort -g "$GENOME" -i "$INTERNAL_BED" | gzip -c > "$INTERNAL_SORTED"

run_intersect () {
  local tag="$1"
  local bfile="$2"
  local out="context/results/overlap/${tag}.tsv"

  [[ -f "$bfile" ]] || { echo "[ERROR] Missing B file: $bfile" | tee -a "$LOG"; exit 1; }

  echo "[INFO] bedtools intersect: $tag" | tee -a "$LOG"

  {
    echo -e "#a_chrom\ta_start\ta_end\ta_name\ta_score\ta_strand\tb_chrom\tb_start\tb_end\tb_name\tb_score\tb_strand"
    bedtools intersect -wa -wb \
      -a <(gzip -dc "$INTERNAL_SORTED") \
      -b <(gzip -dc "$bfile")
  } > "$out"

  echo "[OK] Wrote $out" | tee -a "$LOG"
}

run_intersect_same () {
  local tag="$1"
  local bfile="$2"
  local out="context/results/overlap/${tag}.same_strand.tsv"

  [[ -f "$bfile" ]] || { echo "[ERROR] Missing B file: $bfile" | tee -a "$LOG"; exit 1; }

  echo "[INFO] bedtools intersect (same strand -s): $tag" | tee -a "$LOG"

  {
    echo -e "#a_chrom\ta_start\ta_end\ta_name\ta_score\ta_strand\tb_chrom\tb_start\tb_end\tb_name\tb_score\tb_strand"
    bedtools intersect -s -wa -wb \
      -a <(gzip -dc "$INTERNAL_SORTED") \
      -b <(gzip -dc "$bfile")
  } > "$out"

  echo "[OK] Wrote $out" | tee -a "$LOG"
}

run_pair () {
  local tag="$1"
  local bfile="$2"
  run_intersect "$tag" "$bfile"
  run_intersect_same "$tag" "$bfile"
}

# Gene metadata / sanity
run_pair "internal_vs_genes" "context/assets/genes/sorted/gencode.v48.genes.sorted.bed.gz"

# Core context classification
run_pair "internal_vs_gene_exons"   "context/assets/features/sorted/gencode.v48.gene_exons_merged.sorted.bed.gz"
run_pair "internal_vs_gene_introns" "context/assets/features/sorted/gencode.v48.gene_introns.sorted.bed.gz"

# Extra feature overlaps
run_pair "internal_vs_cds" "context/assets/features/sorted/gencode.v48.cds.sorted.bed.gz"
run_pair "internal_vs_utr" "context/assets/features/sorted/gencode.v48.utr.sorted.bed.gz"

echo "[OK] All overlap intersections done (any-strand + same-strand)." | tee -a "$LOG"
echo "[INFO] Same-strand outputs are: context/results/overlap/*.same_strand.tsv" | tee -a "$LOG"