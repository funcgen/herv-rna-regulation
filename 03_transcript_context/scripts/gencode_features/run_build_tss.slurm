#!/usr/bin/env bash
#SBATCH --job-name=ncrna_tss
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
#SBATCH --output=/path/to/herv-rna-regulation/03_transcript_context/logs/tss_%j.out
# D023: transcript TSS BED and its genome-ordered copy. Streaming Python over
# 388k transcripts plus one bedtools sort of 388k one-base rows (< 0.5 GB).
set -euo pipefail
HERE="$(pwd)"; ROOT="$(cd "${HERE}/../../.." && pwd)"; M="${ROOT}/03_transcript_context"
[[ -x "${HERE}/build_tss.sh" ]] || { echo "ERROR: submit with --chdir set to the gencode_features directory" >&2; exit 1; }
export NCRNA_BEDTOOLS=/path/to/conda/envs/hervarium-ncrna/bin/bedtools
echo "[env] host $(hostname)"
exec "${HERE}/build_tss.sh" --results-root "${M}/results/gencode_features" --genome "${M}/intermediate/grch38.primary.genome"
