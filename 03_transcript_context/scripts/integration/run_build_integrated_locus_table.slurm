#!/usr/bin/env bash
#SBATCH --job-name=ncrna_integrate
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
#SBATCH --output=/path/to/herv-rna-regulation/03_transcript_context/logs/integrate_%j.out
# 03_transcript_context step 9 (effective order 17): integrated HERV-locus table.
# Resources: pandas holds the raw FIMO (716k rows), four closest tables and ten
# overlap tables at once; a few GB at most. 16 GB leaves wide headroom.
set -euo pipefail
HERE="$(pwd)"; ROOT="$(cd "${HERE}/../../.." && pwd)"; M="${ROOT}/03_transcript_context"
[[ -f "${HERE}/build_integrated_locus_table.py" ]] || { echo "ERROR: submit with --chdir set to scripts/integration" >&2; exit 1; }
PY=/path/to/conda/envs/hervarium-ncrna/bin/python
echo "[env] host $(hostname)"; echo "[env] python: ${PY} $("$PY" --version 2>&1)"
echo "[env] pandas $("$PY" -c 'import pandas; print(pandas.__version__)')"
OUT="${M}/results/integration"
mkdir -p "$OUT"
"$PY" "${HERE}/build_integrated_locus_table.py" \
  --internal-bed   "${M}/results/herv_internal/internal_regions.sorted.bed.gz" \
  --domains-tsv    "${ROOT}/01_sources/frozen/herv/HERV_loci_annotated_domains_v5.tsv" \
  --map-tsv        "${ROOT}/01_sources/frozen/ltr/HERV_internal_LTR_map.tsv" \
  --ltr-projection "${ROOT}/07_sparcs/results/ltr_structure/herv_ltr_minimal_projection.tsv.gz" \
  --rbp-fimo       "${ROOT}/02_rbp/results/fimo_primary/fimo.tsv" \
  --closest-dir    "${M}/results/closest" \
  --overlap-dir    "${M}/results/overlap" \
  --rbp-qval 1 \
  --out "${OUT}/transcript_context.tsv"
( cd "$OUT" && sha256sum transcript_context.tsv > checksums.sha256 && cat checksums.sha256 )
echo "[ok] rows: $(($(wc -l < "${OUT}/transcript_context.tsv")-1)); columns: $(head -1 "${OUT}/transcript_context.tsv" | awk -F'\t' '{print NF}')"
