#!/usr/bin/env bash
#SBATCH --job-name=locus_rep
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
#SBATCH --output=/path/to/herv-rna-regulation/04_domains_terminal_exons/logs/locus_rep_%j.out
# 04_domains_terminal_exons step 4 (effective order 33): locus-representative sensitivity.
# Runs twice (second into scratch) and requires byte-identical outputs. 425 rows: < 0.5 GB.
set -euo pipefail
HERE="$(pwd)"; M="$(cd "$HERE/../.." && pwd)"; ROOT="$(cd "$M/.." && pwd)"
[[ -f "$HERE/compare_locus_representative.py" ]] || { echo "ERROR: submit with --chdir set to scripts/conservation_comparison" >&2; exit 1; }
PY=/path/to/conda/envs/hervarium-ncrna/bin/python
COV=$($PY -c "import re;t=open('$ROOT/config/parameters.yaml').read();b=re.search(r'^domains:\s*\n(.*?)(?=^\S|\Z)',t,re.M|re.S).group(1);print(re.search(r'^\s+min_hmm_profile_coverage:\s*([0-9.eE+-]+)',b,re.M).group(1))")
IN=$M/results/domain_transcript_table/all_coverage/domain_centered_unique_domains.tsv
OUT=$M/results/locus_representative
echo "[env] host $(hostname); $($PY -V 2>&1); scipy $($PY -c 'import scipy;print(scipy.__version__)'); pandas $($PY -c 'import pandas;print(pandas.__version__)')"
echo "[params] coverage_threshold=$COV (domains.min_hmm_profile_coverage)"
$PY "$HERE/compare_locus_representative.py" --domains $IN --outdir $OUT --coverage-threshold $COV
RR=$M/validation/conservation_comparison/.rerun_locus; rm -rf $RR
$PY "$HERE/compare_locus_representative.py" --domains $IN --outdir $RR --coverage-threshold $COV > /dev/null
diff <(cd $OUT && sha256sum *) <(cd $RR && sha256sum *) && echo "[ok] re-run byte-identical (4 files)"
rm -rf $RR
( cd $OUT && sha256sum * > checksums.sha256 && cat checksums.sha256 )
