#!/usr/bin/env bash
#SBATCH --job-name=ncrna_terminal_exons
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
#SBATCH --output=/path/to/herv-rna-regulation/03_transcript_context/logs/terminal_exons_%j.out
# 03_transcript_context step 4 (effective order 11): transcript-aware UTR split
# and terminal exon per transcript.
#
# Resources: CDS bounds for ~100k coding transcripts, 392k UTR rows and one best
# exon per transcript (388k) in Python; bedtools sort holds at most the 392k UTR
# rows. Well under 2 GB; 8 GB requested.
set -euo pipefail
HERE="$(pwd)"
ROOT="$(cd "${HERE}/../../.." && pwd)"
[[ -f "${HERE}/build_terminal_exons.py" ]] || { echo "ERROR: submit with --chdir set to the gencode_features directory" >&2; exit 1; }
M="${ROOT}/03_transcript_context"; R="${M}/results/gencode_features"
BT=/path/to/conda/envs/hervarium-ncrna/bin/bedtools
G="${M}/intermediate/grch38.primary.genome"
echo "[env] host $(hostname); $(python3 -V 2>&1); $($BT --version)"
[[ "$(md5sum < "$G" | cut -c1-32)" == f94face851cb2129ab423f54db04fa0c ]] || { echo "ERROR: chrom sizes md5" >&2; exit 1; }
# (no `tee /dev/stderr` here: it reopens and truncates the SLURM log)
( cd "$R" && grep -E '  \./features/sorted/gencode\.v48\.(exons|cds|utr)\.sorted\.bed\.gz$' checksums.sha256 | sha256sum -c - ) \
  || { echo "ERROR: accepted step-3 inputs do not verify" >&2; exit 1; }
echo "[ok] 3 accepted step-3 inputs verified"
exec python3 "${HERE}/build_terminal_exons.py" \
  --exons "${R}/features/sorted/gencode.v48.exons.sorted.bed.gz" \
  --cds   "${R}/features/sorted/gencode.v48.cds.sorted.bed.gz" \
  --utr   "${R}/features/sorted/gencode.v48.utr.sorted.bed.gz" \
  --genome "$G" --outdir "${R}/features/sorted" --bedtools "$BT"
