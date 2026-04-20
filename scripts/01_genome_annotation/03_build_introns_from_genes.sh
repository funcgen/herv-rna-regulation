#!/usr/bin/env bash
set -euo pipefail

LOG="context/logs/02b_build_introns_from_genes.log"
: > "$LOG"

GENOME="input/genome/grch38.primary.genome"
if [[ ! -s "$GENOME" ]]; then
  echo "[ERROR] Genome file not found or empty: $GENOME" | tee -a "$LOG"
  exit 1
fi

echo "[INFO] Building gene-level exons (merged) and introns" | tee -a "$LOG"

TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

GENES="context/assets/genes/gencode.v48.genes.bed.gz"
EXONS="context/assets/features/gencode.v48.exons.bed.gz"

# 1) exons → gene_id + gene_name + gene_type + strand)
zcat "$EXONS" \
  | awk 'BEGIN{OFS="\t"}{
      split($4,a,"|");
      # exon name is: transcript_id|gene_id|gene_name|exon|gene_type|transcript_type
      gene_id=a[2];
      gene_name=a[3];
      gene_type=a[5];
      strand=$6;
      name=gene_id"|"gene_name"|"gene_type;
      print $1,$2,$3,name,$5,strand
    }' \
  > "$TMPDIR/exons.gene.bed6"

# Sort with bedtools + genome order
bedtools sort -g "$GENOME" -i "$TMPDIR/exons.gene.bed6" > "$TMPDIR/exons.gene.sorted.bed6"

# 2) Merge exons per gene (+ strand)
bedtools merge \
  -i "$TMPDIR/exons.gene.sorted.bed6" \
  -c 4,6 -o distinct,distinct \
  > "$TMPDIR/gene_exons.merged.bed6"

awk 'BEGIN{OFS="\t"}{ print $1,$2,$3,$4,0,$5 }' "$TMPDIR/gene_exons.merged.bed6" \
  > "$TMPDIR/gene_exons.merged.bed"

gzip -c "$TMPDIR/gene_exons.merged.bed" > context/assets/features/gencode.v48.gene_exons_merged.bed.gz
echo "[OK] Wrote merged gene exons: context/assets/features/gencode.v48.gene_exons_merged.bed.gz" | tee -a "$LOG"

mkdir -p context/assets/features/sorted
zcat context/assets/features/gencode.v48.gene_exons_merged.bed.gz \
  | bedtools sort -g "$GENOME" -i - \
  | gzip -c > context/assets/features/sorted/gencode.v48.gene_exons_merged.sorted.bed.gz

# 3) Gene-level introns: gene bodies - merged exons
zcat "$GENES" > "$TMPDIR/genes.bed"
bedtools sort -g "$GENOME" -i "$TMPDIR/genes.bed" > "$TMPDIR/genes.sorted.bed"

bedtools subtract \
  -a "$TMPDIR/genes.sorted.bed" \
  -b "$TMPDIR/gene_exons.merged.bed" \
  > "$TMPDIR/gene_introns.raw.bed"

awk 'BEGIN{OFS="\t"}{
  split($4,a,"|");
  newname=a[1]"|"a[2]"|"a[3]"|intron";
  print $1,$2,$3,newname,$5,$6
}' "$TMPDIR/gene_introns.raw.bed" > "$TMPDIR/gene_introns.bed"

bedtools sort -g "$GENOME" -i "$TMPDIR/gene_introns.bed" \
  | gzip -c > context/assets/features/gencode.v48.gene_introns.bed.gz

zcat context/assets/features/gencode.v48.gene_introns.bed.gz \
  | bedtools sort -g "$GENOME" -i - \
  | gzip -c > context/assets/features/sorted/gencode.v48.gene_introns.sorted.bed.gz

N=$(zcat context/assets/features/gencode.v48.gene_introns.bed.gz | wc -l)
echo "[OK] Wrote gene introns: context/assets/features/gencode.v48.gene_introns.bed.gz ($N intervals)" | tee -a "$LOG"
