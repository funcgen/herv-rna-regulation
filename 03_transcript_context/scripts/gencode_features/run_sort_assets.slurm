#!/usr/bin/env bash
#SBATCH --job-name=ncrna_sort_assets
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
#SBATCH --output=/path/to/herv-rna-regulation/03_transcript_context/logs/sort_assets_%j.out
# 03_transcript_context step 3 (effective order 10): genome-order sort, ordering only.
#
# Resources: bedtools sort holds one asset in memory (largest: 2.17M exons,
# ~0.5 GB). The validation gate holds an asset's input and output lines plus
# counters (~5 GB peak on the exons). 16 GB is comfortable; the probe job took
# 2.5 min for the same seven sorts.
set -euo pipefail
HERE="$(pwd)"
ROOT="$(cd "${HERE}/../../.." && pwd)"
[[ -x "${HERE}/sort_assets.sh" ]] || { echo "ERROR: submit with --chdir set to the gencode_features directory" >&2; exit 1; }
export NCRNA_BEDTOOLS="${NCRNA_BEDTOOLS:-/path/to/conda/envs/hervarium-ncrna/bin/bedtools}"
echo "[env] host $(hostname)"
M="${ROOT}/03_transcript_context"
exec "${HERE}/sort_assets.sh" --results-root "${M}/results/gencode_features" --genome "${M}/intermediate/grch38.primary.genome"
