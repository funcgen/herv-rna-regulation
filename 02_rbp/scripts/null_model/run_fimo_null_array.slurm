#!/usr/bin/env bash
#SBATCH --job-name=ncrna_fimo_null
#SBATCH --array=1-5
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
# Array wrapper: one task per shuffled replicate (Stage 7).
#
# Tasks 1-5 map to shuffle01..05. There is deliberately NO task 0: the observed
# dataset was produced by Stage 3 and is immutable, so rescanning it would only
# risk overwriting an accepted artifact. It is represented in the six-dataset
# manifest with generated_in_stage = "Stage 3".
#
# The array index only SELECTS a replicate. Seed, expected hit count and the
# FASTA checksum all come from config, and the producer aborts if the index and
# the configured seed disagree.
#
# Resources: Stage 3 measured MaxRSS 3.88 GB over 26m43s for 716,075 hits on the
# same motif library and background. The shuffles produce slightly FEWER hits
# (~616k), so 8 GB / 2 h is comfortable. The historical 50 GB request was ~13x
# over-provisioned and backfilled poorly.
#
# Submit with --chdir set to this script's directory: SLURM copies the batch
# script into a spool directory, so BASH_SOURCE cannot locate the producer.
set -euo pipefail

HERE="$(pwd)"
[[ -x "${HERE}/run_fimo_null_array.sh" ]] \
  || { echo "ERROR: submit with --chdir set to the null_model script directory" >&2; exit 1; }

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

REPLICATE="$(printf '%02d' "${SLURM_ARRAY_TASK_ID}")"
echo "array task ${SLURM_ARRAY_TASK_ID} -> shuffle${REPLICATE}"

exec ./run_fimo_null_array.sh --replicate "${REPLICATE}"
