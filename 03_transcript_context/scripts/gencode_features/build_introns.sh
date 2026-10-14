#!/usr/bin/env bash
# =============================================================================
# 03_transcript_context step 2 (effective order 9)
#   per-gene merged exons and per-gene introns
# =============================================================================
# Canonical definition (D022, scientific implementation correction):
#   for each gene independently, the union of that gene's own exons across all
#   its transcripts (overlapping or touching intervals merge), and that gene's
#   introns = its gene body minus its own exon union.
#
#   accepted step-1 genes + exons  ->  build_gene_introns.py  (grouped by gene_id)
#     -> results/.../gencode.v48.gene_exons_merged.bed.gz
#     -> results/.../gencode.v48.gene_introns.bed.gz
#
# bedtools is NOT used. The historical workflow ran `bedtools merge` and
# `bedtools subtract` over all exons genome-wide, with no gene or strand
# grouping, which fused exons of different genes and let one gene's exons mask
# another gene's introns (D022). `bedtools merge -s` would not fix that:
# same-strand overlapping genes would still fuse. Grouping must be by gene.
# The historical behaviour is preserved as a forensic reproduction only, in
# validation/historical_reproduction/step2/.
#
# Usage:
#   build_introns.sh --genes F --exons F --genome F --outdir D
# =============================================================================
set -euo pipefail

usage() { echo "usage: $0 --genes F --exons F --genome F --outdir D" >&2; exit 1; }

GENES=""; EXONS=""; GENOME=""; OUTDIR=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --genes)  GENES="$2";  shift 2 ;;
    --exons)  EXONS="$2";  shift 2 ;;
    --genome) GENOME="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    *) echo "unknown argument: $1" >&2; usage ;;
  esac
done
[[ -n "$GENES" && -n "$EXONS" && -n "$GENOME" && -n "$OUTDIR" ]] || usage

for f in "$GENES" "$EXONS" "$GENOME"; do
  [[ -s "$f" ]] || { echo "[error] missing or empty input: $f" >&2; exit 1; }
done

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${NCRNA_PYTHON:-python3}"
echo "[ok] $("$PYTHON" -V 2>&1) at $(command -v "$PYTHON")"

# ---------------------------------------------------------------- accepted inputs
EXPECT_GENES_ROWS=78894
EXPECT_EXONS_ROWS=2165096
n_genes="$(zcat "$GENES" | wc -l)"
n_exons="$(zcat "$EXONS" | wc -l)"
[[ "$n_genes" == "$EXPECT_GENES_ROWS" ]] \
  || { echo "[error] genes BED has ${n_genes} rows, expected ${EXPECT_GENES_ROWS}" >&2; exit 1; }
[[ "$n_exons" == "$EXPECT_EXONS_ROWS" ]] \
  || { echo "[error] exons BED has ${n_exons} rows, expected ${EXPECT_EXONS_ROWS}" >&2; exit 1; }
echo "[ok] step-1 inputs: ${n_genes} genes, ${n_exons} exons"

# ---------------------------------------------------------------- build + publish
# build_gene_introns.py stages both files beside the outputs and renames them
# into place only after both are written.
exec "$PYTHON" "${HERE}/build_gene_introns.py" \
  --genes "$GENES" --exons "$EXONS" --genome "$GENOME" --outdir "$OUTDIR"
