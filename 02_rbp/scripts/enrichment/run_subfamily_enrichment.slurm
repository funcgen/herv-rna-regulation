#!/usr/bin/env bash
#SBATCH --job-name=ncrna_subfam_enrich
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
#SBATCH --output=/path/to/herv-rna-regulation/02_rbp/logs/subfamily_enrichment_%j.out
# Stage 9: subfamily composition-controlled enrichment.
#
# Heavier than Stage 8 (which used 816 MB / 40 s). The dense locus x RBP grid is
# built per subfamily and the largest subfamily has ~9,000 loci x 93 RBPs, and
# the sparse joined count table spans all six datasets over ~1.5M locus/RBP
# combinations. The historical run produced a 28 MB gzipped dense table. 32 GB
# is comfortable headroom without the old 50 GB overprovisioning; the loop frees
# and gc()s each subfamily as it goes.
#
# Submit with --chdir set to this script's directory.
set -euo pipefail

HERE="$(pwd)"
ROOT="$(cd "${HERE}/../../.." && pwd)"
[[ -f "${HERE}/subfamily_enrichment.R" ]] \
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
exec "$NCRNA_RSCRIPT" "${HERE}/subfamily_enrichment.R" \
  --module-dir "${MODULE}" \
  --manifest   "${MODULE}/results/null_model/fimo_dataset_manifest.tsv" \
  --fasta      "${ROOT}/01_sources/frozen/herv/HERV_internal_sequences_v5.fasta" \
  --motifs     "${MODULE}/results/motif_library/attract_hsa_cisbp.meme" \
  --out-dir    "${MODULE}/results/enrichment/subfamily"
