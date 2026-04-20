#!/usr/bin/env bash
set -euo pipefail

GENOME="input/genome/grch38.primary.genome"
LOG="context/logs/04_closest_gene_tss.log"
: > "$LOG"

mkdir -p context/results/closest context/logs

if [[ ! -s "$GENOME" ]]; then
  echo "[ERROR] Genome file not found or empty: $GENOME" | tee -a "$LOG"
  echo "Create it from a .fai: cut -f1,2 GRCh38.fa.fai > $GENOME" | tee -a "$LOG"
  exit 1
fi

INTERNAL_SORTED="context/assets/internal_regions.sorted.bed.gz"
GENES_B="context/assets/genes/sorted/gencode.v48.genes.sorted.bed.gz"
TSS_B="context/assets/tss/sorted/gencode.v48.tss.transcripts.sorted.bed.gz"

[[ -f "$INTERNAL_SORTED" ]] || { echo "[ERROR] Missing $INTERNAL_SORTED" | tee -a "$LOG"; exit 1; }
[[ -f "$GENES_B" ]] || { echo "[ERROR] Missing $GENES_B" | tee -a "$LOG"; exit 1; }
[[ -f "$TSS_B" ]]   || { echo "[ERROR] Missing $TSS_B" | tee -a "$LOG"; exit 1; }

echo "[INFO] bedtools closest: internal -> closest gene" | tee -a "$LOG"
{
  echo -e "#a_chrom\ta_start\ta_end\ta_name\ta_score\ta_strand\tb_chrom\tb_start\tb_end\tb_name\tb_score\tb_strand\tdistance_bp"
  bedtools closest -g "$GENOME" -d \
    -a <(gzip -dc "$INTERNAL_SORTED") \
    -b <(gzip -dc "$GENES_B")
} > context/results/closest/internal_closest_gene.tsv

echo "[OK] Wrote context/results/closest/internal_closest_gene.tsv" | tee -a "$LOG"

echo "[INFO] bedtools closest (same strand): internal -> closest gene" | tee -a "$LOG"
{
  echo -e "#a_chrom\ta_start\ta_end\ta_name\ta_score\ta_strand\tb_chrom\tb_start\tb_end\tb_name\tb_score\tb_strand\tdistance_bp"
  bedtools closest -g "$GENOME" -s -d \
    -a <(gzip -dc "$INTERNAL_SORTED") \
    -b <(gzip -dc "$GENES_B")
} > context/results/closest/internal_closest_gene_same_strand.tsv

echo "[OK] Wrote context/results/closest/internal_closest_gene_same_strand.tsv" | tee -a "$LOG"


echo "[INFO] bedtools closest: internal -> closest TSS" | tee -a "$LOG"
{
  echo -e "#a_chrom\ta_start\ta_end\ta_name\ta_score\ta_strand\tb_chrom\tb_start\tb_end\tb_name\tb_score\tb_strand\tdistance_bp"
  bedtools closest -g "$GENOME" -d \
    -a <(gzip -dc "$INTERNAL_SORTED") \
    -b <(gzip -dc "$TSS_B")
} > context/results/closest/internal_closest_tss.tsv

echo "[OK] Wrote context/results/closest/internal_closest_tss.tsv" | tee -a "$LOG"

echo "[INFO] bedtools closest (same strand): internal -> closest TSS" | tee -a "$LOG"
{
  echo -e "#a_chrom\ta_start\ta_end\ta_name\ta_score\ta_strand\tb_chrom\tb_start\tb_end\tb_name\tb_score\tb_strand\tdistance_bp"
  bedtools closest -g "$GENOME" -s -d \
    -a <(gzip -dc "$INTERNAL_SORTED") \
    -b <(gzip -dc "$TSS_B")
} > context/results/closest/internal_closest_tss_same_strand.tsv

echo "[OK] Wrote context/results/closest/internal_closest_tss_same_strand.tsv" | tee -a "$LOG"



echo "[OK] Closest done." | tee -a "$LOG"
