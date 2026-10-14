#!/usr/bin/env bash
#SBATCH --job-name=ncrna_parse_fimo
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
# Batch wrapper for the primary FIMO coordinate parser.
#
# Measured: ~15 s, MaxRSS 1.4 GB. 8 GB / 20 min is generous and backfills
# immediately. The parser needs no MEME environment — Python standard library only.
#
# Submit with --chdir set to this script's directory. SLURM copies the batch
# script into a spool directory, so the job's working directory is used to
# locate the producer rather than BASH_SOURCE.
set -euo pipefail

HERE="$(pwd)"
[[ -x "${HERE}/run_parse_fimo_primary.sh" ]] \
  || { echo "ERROR: submit with --chdir set to the scanning script directory" >&2; exit 1; }

export NCRNA_SYSTEM="${NCRNA_SYSTEM:-cnag}"

exec ./run_parse_fimo_primary.sh "$@"
