#!/usr/bin/env bash
#SBATCH --job-name=ncrna_herv_id
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
#SBATCH --output=/path/to/herv-rna-regulation/03_transcript_context/logs/herv_id_%j.out
# 03_transcript_context step 10 (effective order 18): attach HERV_id.
# Resources: a text-level pass over a 27 MB table and a 30 MB GTF; well under 1 GB.
set -euo pipefail
HERE="$(pwd)"; ROOT="$(cd "${HERE}/../../.." && pwd)"; M="${ROOT}/03_transcript_context"
[[ -f "${HERE}/add_herv_locus_id.py" ]] || { echo "ERROR: submit with --chdir set to scripts/integration" >&2; exit 1; }
PY=/path/to/conda/envs/hervarium-ncrna/bin/python
echo "[env] host $(hostname); ${PY} $("$PY" --version 2>&1)"
OUT="${M}/results/integration"
( cd "$OUT" && sha256sum -c checksums.sha256 )
"$PY" "${HERE}/add_herv_locus_id.py" --table "${OUT}/transcript_context.tsv" \
  --gtf "${ROOT}/01_sources/derived/hervs_internal_and_ltrs.gtf" --out "${OUT}/transcript_context.with_herv_id.tsv"
( cd "$OUT" && sha256sum transcript_context.tsv transcript_context.with_herv_id.tsv > checksums.sha256 && cat checksums.sha256 )
echo "[ok] rows: $(($(wc -l < "${OUT}/transcript_context.with_herv_id.tsv")-1)); columns: $(head -1 "${OUT}/transcript_context.with_herv_id.tsv" | awk -F'\t' '{print NF}'); md5 $(md5sum < "${OUT}/transcript_context.with_herv_id.tsv" | cut -c1-32); bytes $(stat -c%s "${OUT}/transcript_context.with_herv_id.tsv")"
