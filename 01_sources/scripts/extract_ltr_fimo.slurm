#!/usr/bin/env bash
#SBATCH --job-name=ncrna_extract_ltr_fimo
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=your_partition
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --account=your_slurm_account
# Submit the SRC09 archive extraction as a batch job so it survives session loss.
# Usage:  sbatch --chdir <01_sources/scripts> extract_ltr_fimo.slurm
set -euo pipefail
export NCRNA_SYSTEM=cnag
exec ./extract_ltr_fimo.sh
