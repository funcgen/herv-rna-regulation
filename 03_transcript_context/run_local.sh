#!/usr/bin/env bash
# 03_transcript_context: local step (effective order 19, HERV x gene overlap summaries; see docs/WORKFLOW.md). All
# other 03_transcript_context steps run on SLURM.
# The step is the producer stage of the validated driver named in its "from:" line, copied line for line; the
# validation and historical-comparison stages of that driver are not part of this release.
# Requires: the 03 integrated locus table (effective order 18), NCRNA_RSCRIPT (optional).
# Usage: 03_transcript_context/run_local.sh [step ...]        (no argument: all steps in order)
set -euo pipefail
M="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(cd "$M/.." && pwd)"
STEPS=(gene_overlap_tables)

step_gene_overlap_tables() (
  # from: 03_transcript_context/scripts/integration/run_build_gene_overlap_tables.sh
  HERE="$M/scripts/integration"
  RS="${NCRNA_RSCRIPT:-Rscript}"
  IN="$M/results/integration/transcript_context.with_herv_id.tsv"
  OUT="$M/results/gene_overlap"
  [[ "$(sha256sum < "$IN" | cut -c1-64)" == cbaead1273c916fbcfe4337c7602a9d12c1d5b8bdacfd2a0dd05c0b51f662995 ]] \
    || { echo "[error] input is not the accepted order-18 table" >&2; exit 1; }
  /usr/bin/time -f "[time] producer elapsed %e s, MaxRSS %M KB, exit %x" \
    "$RS" "$HERE/build_gene_overlap_tables.R" --input "$IN" --outdir "$OUT"
  ( cd "$OUT" && find . -type f ! -name checksums.sha256 | LC_ALL=C sort | xargs sha256sum > checksums.sha256 && cat checksums.sha256 )
)

[[ $# -gt 0 ]] || set -- "${STEPS[@]}"
for s in "$@"; do
  declare -F "step_$s" > /dev/null || { echo "unknown step: $s (steps: ${STEPS[*]})" >&2; exit 1; }
  echo "== 03_transcript_context: $s"; "step_$s"
done
