#!/usr/bin/env bash
#SBATCH --job-name=ncrna_introns
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
#SBATCH --output=/path/to/herv-rna-regulation/03_transcript_context/logs/introns_%j.out
# 03_transcript_context step 2 (effective order 9): per-gene merged exons and
# per-gene introns (D022 scientific implementation correction).
#
# Pure Python (standard library); bedtools is not used. The chrom-sizes file is
# still required: it defines the genome order of the outputs.
#
# Resources: 2.17M exon tuples grouped by gene in Python, ~1-2 GB. 8 GB is ample.
set -euo pipefail

HERE="$(pwd)"
ROOT="$(cd "${HERE}/../../.." && pwd)"
[[ -x "${HERE}/build_introns.sh" ]] \
  || { echo "ERROR: submit with --chdir set to the gencode_features directory" >&2; exit 1; }

echo "[env] python3: $(command -v python3) $(python3 -V 2>&1)"
echo "[env] host: $(hostname)"

M="${ROOT}/03_transcript_context"
FEATURES="${M}/results/gencode_features/features"
GENOME="${M}/intermediate/grch38.primary.genome"

# Chrom sizes: a deterministic projection of the declared SRC01 reference,
# cut -f1,2 of its .fai. Derived here rather than copied from the historical
# tree, and checked against the historical file's md5.
if [[ ! -s "$GENOME" ]]; then
  FAI="${NCRNA_GRCH38_FAI:-/path/to/references/GRCh38/genome/GRCh38.primary_assembly.genome.fa.fai}"
  [[ -s "$FAI" ]] || { echo "ERROR: GRCh38 .fai not found: $FAI" >&2; exit 1; }
  mkdir -p "$(dirname "$GENOME")"
  cut -f1,2 "$FAI" > "${GENOME}.tmp"
  got="$(md5sum "${GENOME}.tmp" | cut -d' ' -f1)"
  want="f94face851cb2129ab423f54db04fa0c"
  [[ "$got" == "$want" ]] \
    || { echo "ERROR: derived chrom sizes md5 ${got}, expected ${want}" >&2; rm -f "${GENOME}.tmp"; exit 1; }
  mv "${GENOME}.tmp" "$GENOME"
  echo "[ok] derived chrom sizes from the .fai, md5 ${got}"
fi

exec "${HERE}/build_introns.sh" \
  --genes  "${M}/results/gencode_features/genes/gencode.v48.genes.bed.gz" \
  --exons  "${FEATURES}/gencode.v48.exons.bed.gz" \
  --genome "$GENOME" \
  --outdir "$FEATURES"
