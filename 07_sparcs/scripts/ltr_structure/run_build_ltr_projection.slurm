#!/usr/bin/env bash
#SBATCH --job-name=sparcs_ltr_projection
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
#SBATCH --output=/path/to/herv-rna-regulation/07_sparcs/logs/ltr_projection_%j.out
# 07_sparcs step 3 (frozen order 44), run as effective order 16 (D017).
# Resources: one streaming pass over the 6.27 GB FIMO member plus md5 of all
# three inputs; memory is one counter per LTR (560,594) and the catalogue status
# map, well under 1 GB of Python objects. The first run (4 GB) peaked at MaxRSS
# 4.19 GB, most plausibly page cache of the streamed file charged to the job;
# 8 GB requested for headroom.
set -euo pipefail
HERE="$(pwd)"; ROOT="$(cd "${HERE}/../../.." && pwd)"
[[ -f "${HERE}/build_ltr_projection.py" ]] || { echo "ERROR: submit with --chdir set to 07_sparcs/scripts/ltr_structure" >&2; exit 1; }
L="${ROOT}/01_sources/frozen/ltr"
echo "[env] host $(hostname); $(command -v python3) $(python3 -V 2>&1)"
exec python3 "${HERE}/build_ltr_projection.py" --ltr-bed "${L}/HERV_ltr.bed" \
  --catalogue "${L}/HERV_LTR_U3_R_U5_catalogue.tsv" --fimo "${L}/extracted/merged_fimo_parsed_v4.tsv" \
  --out "${ROOT}/07_sparcs/results/ltr_structure/herv_ltr_minimal_projection.tsv.gz"
