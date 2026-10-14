#!/usr/bin/env bash
#SBATCH --job-name=ncrna_thresh_cmp
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
#SBATCH --output=/path/to/herv-rna-regulation/02_rbp/logs/thresh_cmp_%j.out
# Stage 13: threshold comparison and sensitivity summary.
#
# Table comparison only: the largest input is the 186-row HERVH/HERVK enrichment
# table. No Bioconductor stack is loaded, so this is far lighter than stages
# 10 and 12 (2.0-2.3 GB); 2 GB is ample.
set -euo pipefail
HERE="$(pwd)"
ROOT="$(cd "${HERE}/../../.." && pwd)"
[[ -f "${HERE}/compare_thresholds.R" ]] \
  || { echo "ERROR: submit with --chdir set to the sensitivity_strict directory" >&2; exit 1; }
if [[ -z "${NCRNA_RSCRIPT:-}" ]]; then
  NCRNA_RSCRIPT="$(python3 -c "
import re
t = open('${ROOT}/config/environments/r_stats.yaml').read()
m = re.search(r'^\s*rscript:\s*\"?([^\"\n]+)\"?\s*$', t, re.M)
print(m.group(1).strip() if m else '')")"
fi
[[ -x "$NCRNA_RSCRIPT" ]] || { echo "ERROR: Rscript not found: ${NCRNA_RSCRIPT}" >&2; exit 1; }
M="${ROOT}/02_rbp"
exec "$NCRNA_RSCRIPT" "${HERE}/compare_thresholds.R" \
  --primary-subfamily "${M}/results/enrichment/subfamily" \
  --strict-subfamily  "${M}/results/sensitivity_strict/analysis/subfamily" \
  --primary-global    "${M}/results/enrichment" \
  --strict-global     "${M}/results/sensitivity_strict/analysis/global" \
  --out-root          "${M}/results/sensitivity_strict/analysis"
