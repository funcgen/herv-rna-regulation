#!/usr/bin/env bash
#SBATCH --job-name=ncrna_strict_filter
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
#SBATCH --output=/path/to/herv-rna-regulation/02_rbp/logs/strict_filter_%j.out
# Stage 11: build the six strict p <= 1e-5 datasets.
#
# Pure streaming I/O over ~500 MB of FIMO text, one line at a time, so memory is
# constant and tiny regardless of file size. 2 GB is already generous; the
# implementation never holds a file in memory.
#
# Submit with --chdir set to this script's directory.
set -euo pipefail
HERE="$(pwd)"
[[ -x "${HERE}/run_filter_fimo_strict.sh" ]] \
  || { echo "ERROR: submit with --chdir set to the sensitivity_strict directory" >&2; exit 1; }
export NCRNA_SYSTEM="${NCRNA_SYSTEM:-cnag}"
exec ./run_filter_fimo_strict.sh "$@"
