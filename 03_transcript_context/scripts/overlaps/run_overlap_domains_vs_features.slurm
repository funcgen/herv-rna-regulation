#!/usr/bin/env bash
#SBATCH --job-name=ncrna_domain_ov
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
#SBATCH --output=/path/to/herv-rna-regulation/03_transcript_context/logs/domain_overlaps_%j.out
# 03_transcript_context step 7 (effective order 14): domains x features, -sorted -g -s -wo.
# Resources: chromsweep (-sorted) streams both inputs; 17,787 domains. 4 GB is ample.
set -euo pipefail
HERE="$(pwd)"; ROOT="$(cd "${HERE}/../../.." && pwd)"; M="${ROOT}/03_transcript_context"
[[ -x "${HERE}/overlap_domains_vs_features.sh" ]] || { echo "ERROR: submit with --chdir set to scripts/overlaps" >&2; exit 1; }
export NCRNA_BEDTOOLS=/path/to/conda/envs/hervarium-ncrna/bin/bedtools
echo "[env] host $(hostname)"
exec "${HERE}/overlap_domains_vs_features.sh" --results-root "${M}/results" \
  --domains-bed "${ROOT}/01_sources/frozen/herv/HERV_domains_v5.bed" --genome "${M}/intermediate/grch38.primary.genome"
