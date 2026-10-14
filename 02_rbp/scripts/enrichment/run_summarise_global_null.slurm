#!/usr/bin/env bash
#SBATCH --job-name=ncrna_global_null
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
#SBATCH --output=/path/to/herv-rna-regulation/02_rbp/logs/global_null_%j.out
# Stage 8: global observed-versus-null summary.
#
# Resources: the six FIMO tables total ~500 MB and readr holds them one at a time,
# but bind_rows keeps all six per-locus tables (6 x 53,176 rows) in memory at once.
# That is small; 16 GB is generous for a tidyverse read of ~4M FIMO rows. Far
# lighter than FIMO itself.
#
# Submit with --chdir set to this script's directory.
set -euo pipefail

HERE="$(pwd)"
ROOT="$(cd "${HERE}/../../.." && pwd)"
[[ -f "${HERE}/summarise_global_null.R" ]] \
  || { echo "ERROR: submit with --chdir set to the enrichment script directory" >&2; exit 1; }

# Rscript from config/environments/r_stats.yaml unless overridden.
if [[ -z "${NCRNA_RSCRIPT:-}" ]]; then
  R_YAML="${ROOT}/config/environments/r_stats.yaml"
  NCRNA_RSCRIPT="$(python3 -c "
import re
t = open('${R_YAML}').read()
m = re.search(r'^\s*rscript:\s*\"?([^\"\n]+)\"?\s*$', t, re.M)
print(m.group(1).strip() if m else '')")"
fi
[[ -x "$NCRNA_RSCRIPT" ]] || { echo "ERROR: Rscript not found: ${NCRNA_RSCRIPT}" >&2; exit 1; }

MODULE="${ROOT}/02_rbp"
exec "$NCRNA_RSCRIPT" "${HERE}/summarise_global_null.R" \
  --module-dir "${MODULE}" \
  --manifest "${MODULE}/results/null_model/fimo_dataset_manifest.tsv" \
  --fasta    "${ROOT}/01_sources/frozen/herv/HERV_internal_sequences_v5.fasta" \
  --out-dir  "${MODULE}/results/enrichment"
