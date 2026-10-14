#!/usr/bin/env bash
#SBATCH --job-name=dom_table
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
#SBATCH --output=/path/to/herv-rna-regulation/04_domains_terminal_exons/logs/dom_table_%j.out
# 04_domains_terminal_exons step 1 (effective order 30): the domain-centred transcript table,
# both parameter profiles (conserved, all_coverage) with the same producer.
# Resources: one pass over the 2.17 M-row exon BED into a 388 k-transcript dict and a
# 1,315-row overlap table (~0.5 GB). 4 GB is ample.
set -euo pipefail
HERE="$(pwd)"; M="$(cd "$HERE/../.." && pwd)"; ROOT="$(cd "$M/.." && pwd)"
[[ -f "$HERE/build_domain_table.py" ]] || { echo "ERROR: submit with --chdir set to scripts/domain_transcript_table" >&2; exit 1; }
PY=/path/to/conda/envs/hervarium-ncrna/bin/python
C=$ROOT/03_transcript_context/results; F=$C/gencode_features/features/sorted
OUT=$M/results/domain_transcript_table
echo "[env] host $(hostname); $($PY -V 2>&1); pandas $($PY -c 'import pandas;print(pandas.__version__)')"
for profile in conserved all_coverage; do
  t0=$SECONDS
  $PY "$HERE/build_domain_table.py" --domains-bed $ROOT/01_sources/frozen/herv/HERV_domains_v5.bed \
    --domains-tsv $ROOT/01_sources/frozen/herv/HERV_loci_annotated_domains_v5.tsv \
    --domain-exon-overlaps $C/overlap/domains_vs_gencode_tx_exons.same_strand.wo.tsv \
    --tx-exons-bed $F/gencode.v48.exons.sorted.bed.gz --last-exons-bed $F/gencode.v48.last_exon_per_tx.sorted.bed.gz \
    --params $ROOT/config/parameters.yaml --profile $profile --outdir $OUT
  echo "[run] $profile elapsed $((SECONDS - t0)) s"
done
( cd "$OUT" && find . -type f ! -name checksums.sha256 | LC_ALL=C sort | xargs sha256sum > checksums.sha256 && cat checksums.sha256 )
