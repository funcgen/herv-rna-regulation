#!/usr/bin/env bash
#SBATCH --job-name=sparcs_cand
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=your_partition
#SBATCH --account=your_slurm_account
#SBATCH --output=/path/to/herv-rna-regulation/07_sparcs/logs/sparcs_cand_%j.out
# 07_sparcs step 1 (effective order 43): antisense terminal-LTR candidates, offline (unshare -rn).
# Runs twice (second into scratch) and requires byte-identical outputs.
# Resources: the unchanged historical script on the same tables ran in 38 s (< 1 GB); 8 GB headroom.
set -euo pipefail
HERE="$(pwd)"; M="$(cd "$HERE/../.." && pwd)"; ROOT="$(cd "$M/.." && pwd)"
[[ -f "$HERE/build_sparcs_candidates.py" ]] || { echo "ERROR: submit with --chdir set to scripts/candidates" >&2; exit 1; }
PY=/path/to/conda/envs/hervarium-ncrna/bin/python
R3=$ROOT/03_transcript_context/results; G=$R3/gencode_features
ARGS="--overlap-dir $R3/overlap_ltrs --ltr-bed $R3/herv_ltrs/ltrs.sorted.bed.gz --transcripts-bed $G/transcripts/sorted/gencode.v48.transcripts.sorted.bed.gz --genes-bed $G/genes/sorted/gencode.v48.genes.sorted.bed.gz --params $ROOT/config/parameters.yaml"
OUT=$M/results/candidates
echo "[env] host $(hostname); $($PY -V 2>&1); pandas $($PY -c 'import pandas;print(pandas.__version__)'); numpy $($PY -c 'import numpy;print(numpy.__version__)'); network: none (unshare -rn)"
echo "[in] producer sha256 $(sha256sum < "$HERE/build_sparcs_candidates.py" | cut -c1-64); parameters.yaml sha256 $(sha256sum < $ROOT/config/parameters.yaml | cut -c1-64)"
t0=$SECONDS
unshare -rn $PY "$HERE/build_sparcs_candidates.py" $ARGS --outdir $OUT
echo "[time] producer elapsed $((SECONDS - t0)) s (MaxRSS: sacct)"
RR=$M/validation/candidates/.rerun; rm -rf $RR
unshare -rn $PY "$HERE/build_sparcs_candidates.py" $ARGS --outdir $RR > /dev/null
diff <(cd $OUT && sha256sum *) <(cd $RR && sha256sum *) && echo "[ok] re-run byte-identical ($(ls $OUT | wc -l) files)"
rm -rf $RR "$HERE/__pycache__"
cat $OUT/checksums.sha256
( cd $OUT && sha256sum -c checksums.sha256 )
