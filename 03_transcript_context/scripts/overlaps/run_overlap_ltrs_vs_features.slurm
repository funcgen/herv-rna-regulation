#!/usr/bin/env bash
#SBATCH --job-name=ncrna_ltr_ov
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
#SBATCH --output=/path/to/herv-rna-regulation/03_transcript_context/logs/ltr_overlaps_%j.out
# 03_transcript_context step 8 (effective order 15): LTRs x transcript features, -sorted -g -S -wo.
# Resources: bedtools sort holds 560,594 LTRs (~0.1 GB); -sorted chromsweep then streams. 4 GB ample.
set -euo pipefail
HERE="$(pwd)"; ROOT="$(cd "${HERE}/../../.." && pwd)"; M="${ROOT}/03_transcript_context"
[[ -x "${HERE}/overlap_ltrs_vs_features.sh" ]] || { echo "ERROR: submit with --chdir set to scripts/overlaps" >&2; exit 1; }
export NCRNA_BEDTOOLS=/path/to/conda/envs/hervarium-ncrna/bin/bedtools
echo "[env] host $(hostname)"
exec "${HERE}/overlap_ltrs_vs_features.sh" --results-root "${M}/results" \
  --ltr-bed "${ROOT}/01_sources/frozen/ltr/HERV_ltr.bed" --genome "${M}/intermediate/grch38.primary.genome"
