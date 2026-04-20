#!/usr/bin/env bash
set -euo pipefail

# ---- INPUTS ----
GENOME="input/genome/grch38.primary.genome"

# Your actual file (plain BED, not gz)
DOMAINS_BED="input/annotations/ERV_GyDB_v6_domains.bed"

# Exons are gzipped in your setup
EXONS_BED="context/assets/features/sorted/gencode.v48.gene_exons_merged.sorted.bed.gz"
# ----------------

LOG="context/logs/05_overlap_domains_vs_exons.log"
: > "$LOG"
mkdir -p context/results/overlap context/logs

[[ -s "$GENOME" ]] || { echo "[ERROR] Missing/empty genome: $GENOME" | tee -a "$LOG"; exit 1; }
[[ -s "$DOMAINS_BED" ]] || { echo "[ERROR] Missing/empty domains BED: $DOMAINS_BED" | tee -a "$LOG"; exit 1; }
[[ -s "$EXONS_BED" ]] || { echo "[ERROR] Missing/empty exons BED: $EXONS_BED" | tee -a "$LOG"; exit 1; }

# Helper: stream BED whether gz or plain
stream_bed() {
  local f="$1"
  if [[ "$f" == *.gz ]]; then
    gzip -dc "$f"
  else
    cat "$f"
  fi
}

# Optional: enforce sorting for speed/consistency (recommended for big files)
DOMAINS_SORTED="context/assets/domains.sorted.bed.gz"
EXONS_SORTED="context/assets/exons.sorted.bed.gz"

echo "[INFO] Sorting domains with genome order -> $DOMAINS_SORTED" | tee -a "$LOG"
bedtools sort -g "$GENOME" -i <(stream_bed "$DOMAINS_BED") | gzip -c > "$DOMAINS_SORTED"

echo "[INFO] Sorting exons with genome order -> $EXONS_SORTED" | tee -a "$LOG"
bedtools sort -g "$GENOME" -i <(stream_bed "$EXONS_BED") | gzip -c > "$EXONS_SORTED"

OUT="context/results/overlap/domains_vs_gene_exons.same_strand.wo.tsv"

echo "[INFO] bedtools intersect domains vs exons (same strand, -wo, sorted chromsweep)" | tee -a "$LOG"
echo "       A=$DOMAINS_SORTED" | tee -a "$LOG"
echo "       B=$EXONS_SORTED" | tee -a "$LOG"
echo "       OUT=$OUT" | tee -a "$LOG"

{
  echo -e "#a_chrom\ta_start\ta_end\ta_name\ta_score\ta_strand\tb_chrom\tb_start\tb_end\tb_name\tb_score\tb_strand\toverlap_bp"
  # IMPORTANT: -wo already writes A and B; do NOT add -wa
  bedtools intersect -sorted -g "$GENOME" -s -wo \
    -a <(gzip -dc "$DOMAINS_SORTED") \
    -b <(gzip -dc "$EXONS_SORTED")
} > "$OUT"

echo "[OK] Wrote $OUT" | tee -a "$LOG"