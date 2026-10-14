#!/usr/bin/env bash
#SBATCH --job-name=ncrna_fimo_primary
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
# Batch wrapper for the primary observed FIMO scan.
#
# Resources match the historical run exactly: 1 CPU, 50 GB, 2 h wall
# (~30 min actual). FIMO is single-threaded for this workload, and a tight
# request backfills sooner than an oversized one.
#
# Submit with --chdir set to this script's directory. SLURM copies the batch
# script into a spool directory, so the job's working directory is used to
# locate the producer rather than BASH_SOURCE.
set -euo pipefail

HERE="$(pwd)"
[[ -x "${HERE}/run_fimo_primary.sh" ]] \
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

exec ./run_fimo_primary.sh "$@"
