#!/usr/bin/env bash
# 06_lncrna_function: local steps, in workflow order (effective orders 39-42; see docs/WORKFLOW.md).
# Each step is the producer stage of the validated driver named in its "from:" line, copied line for line; the
# validation, determinism, fixture and historical-comparison stages of those drivers are not part of this release.
# Requires: the frozen PubMed corpus (SRC16) under 01_sources/frozen/pubmed/, NCRNA_NLP_PYTHON (the Python of the NLP
#           environment pinned in config/environments/python_nlp.lock.txt), NCRNA_RSCRIPT (optional), unshare (util-linux)
#           and the 05_lncrna results.
# Usage: 06_lncrna_function/run_local.sh [step ...]        (no argument: all steps in order)
set -euo pipefail
M="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(cd "$M/.." && pwd)"
STEPS=(pubmed_corpus function_statements literature_classification integrated_annotation)

step_pubmed_corpus() (
  # from: 06_lncrna_function/scripts/retrieval/run_fetch_pubmed_abstracts.sh
  HERE="$M/scripts/retrieval"
  SRC=$ROOT/01_sources; CORPUS=$SRC/frozen/pubmed/lncrna_abstracts_dataset_full_no_function_filter.csv
  SUMS=$SRC/manifests/manual_snapshot_checksums.tsv; SNAPS=$SRC/manifests/manual_snapshots.tsv
  OUT=$M/results/pubmed_corpus
  ARGS=(--corpus "$CORPUS" --snapshot-checksums "$SUMS" --snapshots "$SNAPS")
  /usr/bin/time -f "[time] producer elapsed %e s, MaxRSS %M KB, exit %x" python3 "$HERE/fetch_pubmed_abstracts.py" "${ARGS[@]}" --outdir "$OUT"
  ( cd "$OUT" && sha256sum *.tsv > checksums.sha256 && cat checksums.sha256 )
)

step_function_statements() (
  # from: 06_lncrna_function/scripts/text_mining/run_extract_functions.sh
  HERE="$M/scripts/text_mining"
  PY="${NCRNA_NLP_PYTHON:?set NCRNA_NLP_PYTHON to the python of the validated NLP environment}"
  export PYTHONNOUSERSITE=1
  OFF="unshare -rn"
  SRC=$ROOT/01_sources; CORPUS=$SRC/frozen/pubmed/lncrna_abstracts_dataset_full_no_function_filter.csv
  SUMS=$SRC/manifests/manual_snapshot_checksums.tsv; SNAPS=$SRC/manifests/manual_snapshots.tsv
  OUT=$M/results/function_statements
  ARGS=(--corpus "$CORPUS" --snapshot-checksums "$SUMS" --snapshots "$SNAPS")
  LOCK=$ROOT/config/environments/python_nlp.lock.txt
  diff <($PY -m pip freeze 2>/dev/null | LC_ALL=C sort) <(grep -v '^#' "$LOCK" | LC_ALL=C sort) > /dev/null \
    && echo "[ok] interpreter matches $LOCK" || { echo "[error] interpreter does not match the pinned lock $LOCK" >&2; exit 1; }
  /usr/bin/time -f "[time] producer elapsed %e s, MaxRSS %M KB, exit %x" $OFF "$PY" "$HERE/extract_functions.py" "${ARGS[@]}" --outdir "$OUT" 2>&1 | grep -v "Extracting Functions:" || true
  [[ -s "$OUT/lncrna_functions_extracted_v5.csv" ]] || { echo "[error] producer output missing" >&2; exit 1; }
  ( cd "$OUT" && sha256sum lncrna_functions_extracted_v5.csv > checksums.sha256 && cat checksums.sha256 )
)

step_literature_classification() (
  # from: 06_lncrna_function/scripts/text_mining/run_classify_lncrna.sh
  HERE="$M/scripts/text_mining"
  WORK="${TMPDIR:-/tmp}/literature_classification_work.$$"
  mkdir -p "$WORK"
  RS="${NCRNA_RSCRIPT:-Rscript}"; OFF="unshare -rn"
  V5=$M/results/function_statements/lncrna_functions_extracted_v5.csv
  AN=$ROOT/05_lncrna/results/lncrna_annotation/lncRNA_herv_coverage.annotated_with_lncrnawiki.full.tsv
  OUT=$M/results/literature_classification
  /usr/bin/time -f "[time] producer elapsed %e s, MaxRSS %M KB, exit %x" $OFF "$RS" "$HERE/classify_lncrna.R" --v5 "$V5" --annotation "$AN" --outdir "$OUT" > "$WORK/producer.stdout" 2> >(grep -v "^\s*$" >&2)
  grep -E "^\[" "$WORK/producer.stdout" || true
  ( cd "$OUT" && sha256sum *.tsv > checksums.sha256 && cat checksums.sha256 )
  rm -rf "$WORK"
)

step_integrated_annotation() (
  # from: 06_lncrna_function/scripts/integration/run_integrate_annotation.sh
  HERE="$M/scripts/integration"
  RS="${NCRNA_RSCRIPT:-Rscript}"; OFF="unshare -rn"
  AN=$ROOT/05_lncrna/results/lncrna_annotation/lncRNA_herv_coverage.annotated_with_lncrnawiki.full.tsv
  LS=$M/results/literature_classification/lit_summary_per_gene.tsv
  CA=$ROOT/05_lncrna/results/coverage_categories/lncRNA_herv_coverage_categories.tsv
  OUT=$M/results/integrated_annotation
  ARGS=(--annotation "$AN" --literature-summary "$LS" --categories "$CA")
  /usr/bin/time -f "[time] producer elapsed %e s, MaxRSS %M KB, exit %x" $OFF "$RS" "$HERE/integrate_annotation.R" "${ARGS[@]}" --outdir "$OUT" 2>&1 | grep -E "^\[|time|Error"
  ( cd "$OUT" && sha256sum *.tsv > checksums.sha256 && cat checksums.sha256 )
)

[[ $# -gt 0 ]] || set -- "${STEPS[@]}"
for s in "$@"; do
  declare -F "step_$s" > /dev/null || { echo "unknown step: $s (steps: ${STEPS[*]})" >&2; exit 1; }
  echo "== 06_lncrna_function: $s"; "step_$s"
done
