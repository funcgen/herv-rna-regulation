#!/usr/bin/env bash
set -euo pipefail

GENOME="input/genome/grch38.primary.genome"
LOG="context/logs/02_sort_context_assets.log"
: > "$LOG"

mkdir -p \
  context/assets/genes/sorted \
  context/assets/features/sorted \
  context/assets/transcripts/sorted \
  context/assets/tss/sorted

if [[ ! -s "$GENOME" ]]; then
  echo "[ERROR] Genome file not found or empty: $GENOME" | tee -a "$LOG"
  echo "Create it from a .fai: cut -f1,2 GRCh38.fa.fai > $GENOME" | tee -a "$LOG"
  exit 1
fi

sort_bed_gz () {
  local in_gz="$1"
  local out_gz="$2"
  echo "[INFO] bedtools sort $in_gz -> $out_gz" | tee -a "$LOG"
  zcat "$in_gz" \
    | bedtools sort -g "$GENOME" -i - \
    | gzip -c > "$out_gz"
}

# genes/features
sort_bed_gz context/assets/genes/gencode.v48.genes.bed.gz \
            context/assets/genes/sorted/gencode.v48.genes.sorted.bed.gz

sort_bed_gz context/assets/transcripts/gencode.v48.transcripts.bed.gz \
            context/assets/transcripts/sorted/gencode.v48.transcripts.sorted.bed.gz

sort_bed_gz context/assets/features/gencode.v48.exons.bed.gz \
            context/assets/features/sorted/gencode.v48.exons.sorted.bed.gz

sort_bed_gz context/assets/features/gencode.v48.cds.bed.gz \
            context/assets/features/sorted/gencode.v48.cds.sorted.bed.gz

sort_bed_gz context/assets/features/gencode.v48.utr.bed.gz \
            context/assets/features/sorted/gencode.v48.utr.sorted.bed.gz

# Add these two if you want gene-level exon union available sorted
if [[ -f context/assets/features/gencode.v48.gene_exons_merged.bed.gz ]]; then
  sort_bed_gz context/assets/features/gencode.v48.gene_exons_merged.bed.gz \
              context/assets/features/sorted/gencode.v48.gene_exons_merged.sorted.bed.gz
fi

if [[ -f context/assets/features/gencode.v48.gene_introns.bed.gz ]]; then
  sort_bed_gz context/assets/features/gencode.v48.gene_introns.bed.gz \
              context/assets/features/sorted/gencode.v48.gene_introns.sorted.bed.gz
fi

# TSS (plain bed)
echo "[INFO] bedtools sort TSS bed..." | tee -a "$LOG"
bedtools sort -g "$GENOME" -i input/annotations/tss_gencode.transcripts.bed \
  | gzip -c > context/assets/tss/sorted/gencode.v48.tss.transcripts.sorted.bed.gz

echo "[OK] Done sorting assets." | tee -a "$LOG"
