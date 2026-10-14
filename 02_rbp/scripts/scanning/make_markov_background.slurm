#!/usr/bin/env bash
#SBATCH --job-name=ncrna_markov_background
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
# Batch wrapper for the order-1 background build.
#
# Submit with --chdir set to this script's directory, e.g.
#   sbatch --chdir="$(dirname "$0")" make_markov_background.slurm
# SLURM copies the batch script into a spool directory, so the script cannot
# locate itself via BASH_SOURCE; the job's working directory is used instead.
#
# The MEME location comes from config/environments/motif.yaml (env_path). No
# site- or user-specific path is hard-coded here. Override with NCRNA_MEME_BIN.
set -euo pipefail

HERE="$(pwd)"
[[ -x "${HERE}/make_markov_background.sh" ]] \
  || { echo "ERROR: submit with --chdir set to the scanning script directory" >&2; exit 1; }

export NCRNA_SYSTEM="${NCRNA_SYSTEM:-cnag}"

if [[ -z "${NCRNA_MEME_BIN:-}" ]]; then
  MOTIF_YAML="${HERE}/../../../config/environments/motif.yaml"
  if [[ -f "$MOTIF_YAML" ]]; then
    ENV_PATH="$(python3 -c "
import re
text = open('${MOTIF_YAML}').read()
m = re.search(r'^\s*env_path:\s*\"?([^\"\n]+)\"?\s*$', text, re.M)
print(m.group(1).strip() if m else '')")"
    [[ -n "$ENV_PATH" ]] && export NCRNA_MEME_BIN="${ENV_PATH}/bin"
  fi
fi

exec ./make_markov_background.sh "$@"
