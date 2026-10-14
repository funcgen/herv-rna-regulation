#!/usr/bin/env bash
# 05_lncrna: local steps, in workflow order (effective orders 35-38; see docs/WORKFLOW.md).
# Each step is the producer stage of the validated driver named in its "from:" line, copied line for line; the
# validation, determinism, fixture and historical-comparison stages of those drivers are not part of this release.
# Requires: NCRNA_GENCODE_V48_GTF (steps lncrna_coverage, pc_coverage), NCRNA_RSCRIPT (optional, default Rscript),
#           the 03_transcript_context results and the dated snapshots SRC14/SRC17 (step lncrna_annotation).
# Usage: 05_lncrna/run_local.sh [step ...]        (no argument: all steps in order)
set -euo pipefail
M="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(cd "$M/.." && pwd)"
STEPS=(lncrna_coverage pc_coverage coverage_categories lncrna_annotation)

step_lncrna_coverage() (
  # from: 05_lncrna/scripts/run_compute_lncrna_herv_coverage.sh
  HERE="$M/scripts"
  GTF="${NCRNA_GENCODE_V48_GTF:?set NCRNA_GENCODE_V48_GTF to the SRC02 GTF}"
  GTF_SHA="$(awk -F'\t' '$1=="SRC02"' "$ROOT/01_sources/SOURCE_STATUS.tsv" | grep -oE 'sha256 of the uncompressed GTF[^:]*: [0-9a-f]{64}' | grep -oE '[0-9a-f]{64}')"
  [[ ${#GTF_SHA} == 64 ]] || { echo "[error] SRC02 sha256 not found in SOURCE_STATUS.tsv" >&2; exit 1; }
  GO=$ROOT/03_transcript_context/results/gene_overlap
  IN=$GO/tables/lncRNA_overlapping_HERVs.gene_level.tsv
  OUT=$M/results/herv_coverage
  /usr/bin/time -f "[time] producer elapsed %e s, MaxRSS %M KB, exit %x" \
    python3 "$HERE/compute_lncrna_herv_coverage.py" --overlap-table "$IN" --overlap-manifest "$GO/checksums.sha256" \
      --gtf "$GTF" --gtf-sha256 "$GTF_SHA" --outdir "$OUT"
  ( cd "$OUT" && sha256sum lncRNA_herv_coverage.tsv > checksums.sha256 && cat checksums.sha256 )
)

step_pc_coverage() (
  # from: 05_lncrna/scripts/run_compute_pc_herv_coverage.sh
  HERE="$M/scripts"
  GTF="${NCRNA_GENCODE_V48_GTF:?set NCRNA_GENCODE_V48_GTF to the SRC02 GTF}"
  GTF_SHA="$(awk -F'\t' '$1=="SRC02"' "$ROOT/01_sources/SOURCE_STATUS.tsv" | grep -oE 'sha256 of the uncompressed GTF[^:]*: [0-9a-f]{64}' | grep -oE '[0-9a-f]{64}')"
  [[ ${#GTF_SHA} == 64 ]] || { echo "[error] SRC02 sha256 not found in SOURCE_STATUS.tsv" >&2; exit 1; }
  GO=$ROOT/03_transcript_context/results/gene_overlap
  IN=$GO/tables/protein_coding_overlapping_HERVs.gene_level.tsv
  OUT=$M/results/pc_herv_coverage
  TABLE=protein_coding_overlapping_herv_coverage.tsv
  /usr/bin/time -f "[time] producer elapsed %e s, MaxRSS %M KB, exit %x" \
    python3 "$HERE/compute_pc_herv_coverage.py" --overlap-table "$IN" --overlap-manifest "$GO/checksums.sha256" \
      --gtf "$GTF" --gtf-sha256 "$GTF_SHA" --outdir "$OUT"
  ( cd "$OUT" && sha256sum $TABLE > checksums.sha256 && cat checksums.sha256 )
)

step_coverage_categories() (
  # from: 05_lncrna/scripts/run_assign_coverage_categories.sh
  HERE="$M/scripts"
  RS="${NCRNA_RSCRIPT:-Rscript}"
  COV=$M/results/herv_coverage/lncRNA_herv_coverage.tsv; PARAMS=$ROOT/config/parameters.yaml
  OUT=$M/results/coverage_categories
  /usr/bin/time -f "[time] producer elapsed %e s, MaxRSS %M KB, exit %x" \
    "$RS" "$HERE/assign_coverage_categories.R" --coverage "$COV" --coverage-manifest "$M/results/herv_coverage/checksums.sha256" \
      --params "$PARAMS" --outdir "$OUT"
  ( cd "$OUT" && sha256sum lncRNA_herv_coverage_categories.tsv lncRNA_herv_coverage_category_counts.tsv > checksums.sha256 && cat checksums.sha256 )
)

step_lncrna_annotation() (
  # from: 05_lncrna/scripts/run_annotate_lncrna_tables.sh
  HERE="$M/scripts"
  RS="${NCRNA_RSCRIPT:-Rscript}"
  SRC=$ROOT/01_sources; JSON=$SRC/frozen/lncrnawiki/lncrnawiki.json
  GO=$SRC/frozen/ensembl_go/ensembl_biomart_go_per_gene.recovered_2026-04-16.tsv; SUMS=$SRC/manifests/manual_snapshot_checksums.tsv
  COV=$M/results/herv_coverage/lncRNA_herv_coverage.tsv; CAT=$M/results/coverage_categories/lncRNA_herv_coverage_categories.tsv
  OUT=$M/results/lncrna_annotation
  ARGS=(--coverage "$COV" --categories "$CAT" --lncrnawiki-json "$JSON" --go-snapshot "$GO" --snapshot-checksums "$SUMS")
  /usr/bin/time -f "[time] producer elapsed %e s, MaxRSS %M KB, exit %x" "$RS" "$HERE/annotate_lncrna_tables.R" "${ARGS[@]}" --outdir "$OUT" 2>&1 | grep -v "^  "
  ( cd "$OUT" && sha256sum *.tsv > checksums.sha256 && cat checksums.sha256 )
)

[[ $# -gt 0 ]] || set -- "${STEPS[@]}"
for s in "$@"; do
  declare -F "step_$s" > /dev/null || { echo "unknown step: $s (steps: ${STEPS[*]})" >&2; exit 1; }
  echo "== 05_lncrna: $s"; "step_$s"
done
