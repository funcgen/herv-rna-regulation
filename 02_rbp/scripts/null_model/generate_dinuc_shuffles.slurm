#!/usr/bin/env bash
#SBATCH --job-name=ncrna_dinuc_shuf
#SBATCH --array=1-5
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
# Array wrapper: one task per replicate, mirroring the historical execution
# model (run_06_make_dinuc_shuffles_array.sh, array 1-5, job 1085376).
#
# The array index only SELECTS a replicate. The seed comes from the explicit
# table in config/parameters.yaml, and generate_dinuc_shuffles.sh refuses to run
# if that table and base_seed + replicate disagree. Renumbering this array can
# therefore not silently reassign a seed.
#
# Submit with --chdir set to this script's directory: SLURM copies the batch
# script into a spool directory, so BASH_SOURCE cannot locate the producer.
set -euo pipefail

HERE="$(pwd)"
[[ -x "${HERE}/generate_dinuc_shuffles.sh" ]] \
  || { echo "ERROR: submit with --chdir set to the null_model script directory" >&2; exit 1; }

export NCRNA_SYSTEM="${NCRNA_SYSTEM:-cnag}"

REPLICATE="$(printf '%02d' "${SLURM_ARRAY_TASK_ID}")"
echo "array task ${SLURM_ARRAY_TASK_ID} -> replicate ${REPLICATE}"

exec ./generate_dinuc_shuffles.sh --replicate "${REPLICATE}" "$@"
