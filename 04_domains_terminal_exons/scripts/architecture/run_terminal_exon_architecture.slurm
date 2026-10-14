#!/usr/bin/env bash
#SBATCH --job-name=te_arch
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
#SBATCH --output=/path/to/herv-rna-regulation/04_domains_terminal_exons/logs/te_arch_%j.out
# 04_domains_terminal_exons step 5 (effective order 34): terminal-exon architecture (Figure 2D
# source tables, VAL37). Runs twice (second into scratch) and requires byte-identical outputs.
# Resources: the historical run on the same inputs peaked at 1.86 GB MaxRSS in 95 s; 16 GB headroom.
set -euo pipefail
HERE="$(pwd)"; M="$(cd "$HERE/../.." && pwd)"; ROOT="$(cd "$M/.." && pwd)"
[[ -f "$HERE/terminal_exon_architecture.py" ]] || { echo "ERROR: submit with --chdir set to scripts/architecture" >&2; exit 1; }
PY=/path/to/conda/envs/hervarium-ncrna/bin/python
F=$ROOT/03_transcript_context/results/gencode_features/features/sorted
A=$M/results/domain_transcript_table/all_coverage/domain_centered_qualifying_domain_transcript_associations.tsv
ARGS="--associations $A --tx-exons-bed $F/gencode.v48.exons.sorted.bed.gz --last-exons-bed $F/gencode.v48.last_exon_per_tx.sorted.bed.gz --params $ROOT/config/parameters.yaml"
OUT=$M/results/terminal_exon_architecture
echo "[env] host $(hostname); $($PY -V 2>&1); pandas $($PY -c 'import pandas;print(pandas.__version__)'); numpy $($PY -c 'import numpy;print(numpy.__version__)'); scipy $($PY -c 'import scipy;print(scipy.__version__)')"
t0=$SECONDS
$PY "$HERE/terminal_exon_architecture.py" $ARGS --outdir $OUT
echo "[run] elapsed $((SECONDS - t0)) s"
RR=$M/validation/architecture/.rerun; rm -rf $RR
$PY "$HERE/terminal_exon_architecture.py" $ARGS --outdir $RR > /dev/null
diff <(cd $OUT && sha256sum *) <(cd $RR && sha256sum *) && echo "[ok] re-run byte-identical (9 files)"
rm -rf $RR
( cd $OUT && sha256sum *.tsv > checksums.sha256 && cat checksums.sha256 )
