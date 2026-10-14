#!/usr/bin/env bash
#SBATCH --job-name=ncrna_motif_library
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
# Batch wrapper for the motif-library build.
#
# Submit with --chdir set to this script's directory. SLURM copies the batch
# script into a spool directory, so the job's working directory is used to
# locate the producer rather than BASH_SOURCE.
set -euo pipefail

[[ -x "./build_attract_meme.sh" ]] \
  || { echo "ERROR: submit with --chdir set to the motif_library script directory" >&2; exit 1; }

export NCRNA_SYSTEM="${NCRNA_SYSTEM:-cnag}"
exec ./build_attract_meme.sh "$@"
