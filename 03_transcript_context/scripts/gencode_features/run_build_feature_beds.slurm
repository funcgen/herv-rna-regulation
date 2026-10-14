#!/usr/bin/env bash
#SBATCH --job-name=ncrna_gencode_beds
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
#SBATCH --output=/path/to/herv-rna-regulation/03_transcript_context/logs/gencode_beds_%j.out
# 03_transcript_context step 1 (effective order 8): extract GENCODE feature BEDs.
#
# The producer STREAMS the 1.8 GB GTF and writes each record as it is read, so
# memory is constant and small regardless of input size. The historical version
# buffered ~4M records in memory and needed several GB; 4 GB here is already
# generous. The cost is dominated by gzip compression, not RAM.
#
# Submit with --chdir set to this script's directory.
set -euo pipefail

HERE="$(pwd)"
ROOT="$(cd "${HERE}/../../.." && pwd)"
[[ -f "${HERE}/build_feature_beds.py" ]] \
  || { echo "ERROR: submit with --chdir set to the gencode_features directory" >&2; exit 1; }

GTF="${NCRNA_GENCODE_V48_GTF:-/path/to/references/GRCh38/annotations/gencode.v48.primary_assembly.annotation.gtf}"
EXPECT_SHA="0bc33e92c494f98f0ad91b0c4cde56f7caafe8e62d270a9f8554d0991765a155"

exec python3 "${HERE}/build_feature_beds.py" \
  --gtf "$GTF" \
  --outdir "${ROOT}/03_transcript_context/results/gencode_features" \
  --prefix "gencode.v48" \
  --expect-gtf-sha256 "$EXPECT_SHA"
