#!/usr/bin/env bash
#SBATCH --job-name=ncrna_closest
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
#SBATCH --output=/path/to/herv-rna-regulation/03_transcript_context/logs/closest_%j.out
# 03_transcript_context step 6 (effective order 13): nearest gene and TSS.
# Resources: four bedtools closest runs; B is at most 388k one-base TSS rows
# (~0.1 GB in memory). 4 GB is ample.
set -euo pipefail
HERE="$(pwd)"; ROOT="$(cd "${HERE}/../../.." && pwd)"; M="${ROOT}/03_transcript_context"
[[ -x "${HERE}/closest_gene_tss.sh" ]] || { echo "ERROR: submit with --chdir set to scripts/overlaps" >&2; exit 1; }
export NCRNA_BEDTOOLS=/path/to/conda/envs/hervarium-ncrna/bin/bedtools
echo "[env] host $(hostname)"; echo "[env] which bedtools: $(which "$NCRNA_BEDTOOLS")"; echo "[env] bedtools --version: $("$NCRNA_BEDTOOLS" --version)"
exec "${HERE}/closest_gene_tss.sh" --results-root "${M}/results" --genome "${M}/intermediate/grch38.primary.genome"
