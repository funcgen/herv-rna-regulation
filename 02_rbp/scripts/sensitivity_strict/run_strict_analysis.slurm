#!/usr/bin/env bash
#SBATCH --job-name=ncrna_strict_analysis
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
#SBATCH --output=/path/to/herv-rna-regulation/02_rbp/logs/strict_analysis_%j.out
# Stage 12: strict p <= 1e-5 re-analysis, reusing the Stage 8-10 code.
#
# The strict inputs are ~1.7 MB each against ~85 MB primary, so this is far
# lighter than the primary runs (which peaked at 0.8, 0.9 and 2.0 GB). The
# 2.0 GB Stage 10 peak was the Bioconductor annotation stack, which is
# dataset-independent, so 8 GB covers the whole chain.
set -euo pipefail
HERE="$(pwd)"
[[ -x "${HERE}/run_strict_analysis.sh" ]] \
  || { echo "ERROR: submit with --chdir set to the sensitivity_strict directory" >&2; exit 1; }
export NCRNA_SYSTEM="${NCRNA_SYSTEM:-cnag}"
exec ./run_strict_analysis.sh "$@"
