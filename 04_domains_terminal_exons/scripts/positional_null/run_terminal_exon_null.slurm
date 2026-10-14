#!/usr/bin/env bash
#SBATCH --job-name=pos_null
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
#SBATCH --output=/path/to/herv-rna-regulation/04_domains_terminal_exons/logs/pos_null_%j.out
# 04_domains_terminal_exons step 2 (effective order 31): transcript-aware positional null.
# Runs the canonical script twice (the second run into a scratch directory) and requires
# byte-identical outputs. Resources: < 0.3 GB measured for the historical v2; 4 GB.
set -euo pipefail
HERE="$(pwd)"; M="$(cd "$HERE/../.." && pwd)"; ROOT="$(cd "$M/.." && pwd)"
[[ -f "$HERE/terminal_exon_null.py" ]] || { echo "ERROR: submit with --chdir set to scripts/positional_null" >&2; exit 1; }
PY=/path/to/conda/envs/hervarium-ncrna/bin/python
param() {  # value of domains.<key> in config/parameters.yaml (no YAML module on CNAG)
  $PY -c "import re,sys;t=open('$ROOT/config/parameters.yaml').read();b=re.search(r'^domains:\s*\n(.*?)(?=^\S|\Z)',t,re.M|re.S).group(1);m=re.search(r'^\s+$1:\s*([0-9.eE+-]+)',b,re.M);sys.exit('missing domains.$1') if not m else print(m.group(1))"
}
OV=$(param min_exonic_overlap_fraction); SEED=$(param positional_null_seed); NP=$(param positional_null_n_permutations)
IN=$M/results/domain_transcript_table/conserved
EX=$ROOT/03_transcript_context/results/gencode_features/features/sorted/gencode.v48.exons.sorted.bed.gz
OUT=$M/results/positional_null
echo "[env] host $(hostname); $($PY -V 2>&1); numpy $($PY -c 'import numpy;print(numpy.__version__)'); pandas $($PY -c 'import pandas;print(pandas.__version__)')"
echo "[params] min_exon_overlap=$OV seed=$SEED n_permutations=$NP scope=both"
run() {
  $PY "$HERE/terminal_exon_null.py" --domains $IN/domain_centered_unique_domains.tsv \
    --associations $IN/domain_centered_qualifying_domain_transcript_associations.tsv \
    --tx-exons-bed $EX --outdir "$1" --min-exon-overlap $OV --n-permutations $NP --seed $SEED --scope both
}
t0=$SECONDS; run $OUT; echo "[run] canonical elapsed $((SECONDS - t0)) s"
RR=$M/validation/positional_null/.rerun; rm -rf $RR
t0=$SECONDS; run $RR > /dev/null; echo "[run] rerun elapsed $((SECONDS - t0)) s"
if diff <(cd $OUT && sha256sum *) <(cd $RR && sha256sum *); then echo "[ok] re-run byte-identical (5 files)"; else echo "[error] re-run differs"; exit 1; fi
rm -rf $RR
( cd $OUT && sha256sum * > checksums.sha256 && cat checksums.sha256 )
