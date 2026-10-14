#!/usr/bin/env bash
#SBATCH --job-name=ncrna_go_gsea
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
#SBATCH --output=/path/to/herv-rna-regulation/02_rbp/logs/go_gsea_%j.out
# Stage 10: HERVH/HERVK profiles, GO ORA and pre-ranked GSEA.
#
# Resources: the analysis operates on 186 Stage 9 rows and a 93-RBP universe, so
# the data are tiny. The memory goes to the Bioconductor annotation stack -
# org.Hs.eg.db and GO.db are SQLite-backed and clusterProfiler/fgsea hold GO
# term membership in memory. Stages 8 and 9 peaked below 900 MB on far larger
# tables; 8 GB is ample headroom for the annotation objects without the
# unevidenced 32-50 GB requests this release has been trimming.
#
# Submit with --chdir set to this script's directory.
set -euo pipefail

HERE="$(pwd)"
ROOT="$(cd "${HERE}/../../.." && pwd)"
[[ -f "${HERE}/hervh_hervk_go_gsea.R" ]] \
  || { echo "ERROR: submit with --chdir set to the enrichment script directory" >&2; exit 1; }

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
exec "$NCRNA_RSCRIPT" "${HERE}/hervh_hervk_go_gsea.R" \
  --stage9-dir "${MODULE}/results/enrichment/subfamily" \
  --out-dir    "${MODULE}/results/enrichment/hervh_hervk_followup"
