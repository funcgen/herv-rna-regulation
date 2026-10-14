#!/usr/bin/env bash
#SBATCH --job-name=ncrna_internal_ov
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
#SBATCH --output=/path/to/herv-rna-regulation/03_transcript_context/logs/internal_overlaps_%j.out
# 03_transcript_context step 5 (effective order 12): HERV internal x gene features.
# Resources: ten bedtools intersects (no -sorted) of 53,176 HERVs against at most
# 902k feature rows (CDS); bedtools holds B in memory (~0.3 GB). The semantic-audit
# job ran four of these within minutes. 8 GB is ample.
set -euo pipefail
HERE="$(pwd)"
ROOT="$(cd "${HERE}/../../.." && pwd)"
[[ -x "${HERE}/overlap_internal_vs_features.sh" ]] || { echo "ERROR: submit with --chdir set to scripts/overlaps" >&2; exit 1; }
export NCRNA_BEDTOOLS=/path/to/conda/envs/hervarium-ncrna/bin/bedtools
echo "[env] host $(hostname)"
echo "[env] which bedtools: $(which "$NCRNA_BEDTOOLS")"
echo "[env] bedtools --version: $("$NCRNA_BEDTOOLS" --version)"
M="${ROOT}/03_transcript_context"
exec "${HERE}/overlap_internal_vs_features.sh" --results-root "${M}/results" \
  --herv-bed "${ROOT}/01_sources/frozen/herv/HERV_internal_v5.bed" --genome "${M}/intermediate/grch38.primary.genome"
