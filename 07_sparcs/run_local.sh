#!/usr/bin/env bash
# 07_sparcs: local steps, in workflow order (effective orders 44-46; see docs/WORKFLOW.md). Steps 16 and 43 of this
# module run on SLURM (scripts/ltr_structure/run_build_ltr_projection.slurm, scripts/candidates/run_build_sparcs_candidates.slurm).
# Each step is the producer stage of the validated driver named in its "from:" line, copied line for line; the
# validation, determinism, fixture and historical-comparison stages of those drivers are not part of this release.
# Inputs passed to the drivers as options are read from these variables (defaults: the 01_sources locations):
#   NCRNA_LTR_FIMO        extracted SRC09 member merged_fimo_parsed_v4.tsv (01_sources/scripts/extract_ltr_fimo.sh)
#   NCRNA_U3R_CATALOGUE   SRC08 HERV_LTR_U3_R_U5_catalogue.tsv
# Requires: NCRNA_GENCODE_V48_GTF (step enrichment), NCRNA_RSCRIPT (optional), unshare (util-linux), the frozen MSigDB C7
#           snapshot SRC18 (step enrichment; 01_sources/scripts/freeze_msigdb_c7.R).
# Usage: 07_sparcs/run_local.sh [step ...]        (no argument: all steps in order)
set -euo pipefail
M="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(cd "$M/.." && pwd)"
STEPS=(interferon_motifs enrichment u3ru5_integrity)

step_interferon_motifs() (
  # from: 07_sparcs/scripts/interferon_motifs/run_summarise_ifn_motifs.sh
  HERE="$M/scripts/interferon_motifs"
  FIMO="${NCRNA_LTR_FIMO:-$ROOT/01_sources/frozen/ltr/extracted/merged_fimo_parsed_v4.tsv}"
  RS="${NCRNA_RSCRIPT:-Rscript}"; OFF="unshare -rn"
  ARGS=(--fimo "$FIMO" --archive-members "$ROOT/01_sources/manifests/archive_members.tsv" --params "$ROOT/config/parameters.yaml")
  OUT=$M/results/interferon_motifs
  /usr/bin/time -f "[time] producer elapsed %e s, MaxRSS %M KB, exit %x" $OFF "$RS" "$HERE/summarise_ifn_motifs.R" "${ARGS[@]}" --outdir "$OUT" 2>&1 | grep -E "^\[|time|Error"
  cat "$OUT/checksums.sha256"
)

step_enrichment() (
  # from: 07_sparcs/scripts/enrichment/run_analyse_sparcs.sh
  HERE="$M/scripts/enrichment"
  FIMO="${NCRNA_LTR_FIMO:-$ROOT/01_sources/frozen/ltr/extracted/merged_fimo_parsed_v4.tsv}"
  GTF="${NCRNA_GENCODE_V48_GTF:?set NCRNA_GENCODE_V48_GTF to the SRC02 GTF}"
  RS="${NCRNA_RSCRIPT:-Rscript}"; OFF="unshare -rn"
  ARGS=(--candidates-dir "$M/results/candidates" --ifn-dir "$M/results/interferon_motifs" --gtf "$GTF" --source-status "$ROOT/01_sources/SOURCE_STATUS.tsv"
        --fimo "$FIMO" --archive-members "$ROOT/01_sources/manifests/archive_members.tsv" --msigdb-dir "$ROOT/01_sources/frozen/msigdb"
        --snapshot-checksums "$ROOT/01_sources/manifests/manual_snapshot_checksums.tsv" --params "$ROOT/config/parameters.yaml")
  OUT=$M/results/enrichment
  /usr/bin/time -f "[time] producer elapsed %e s, MaxRSS %M KB, exit %x" $OFF "$RS" "$HERE/analyse_sparcs.R" "${ARGS[@]}" --outdir "$OUT" 2>&1 | grep -E "^\[|time|Error"
  cat "$OUT/checksums.sha256"
)

step_u3ru5_integrity() (
  # from: 07_sparcs/scripts/ltr_structure/run_compare_u3ru5_integrity.sh
  HERE="$M/scripts/ltr_structure"
  U3R="${NCRNA_U3R_CATALOGUE:-$ROOT/01_sources/frozen/ltr/HERV_LTR_U3_R_U5_catalogue.tsv}"
  RS="${NCRNA_RSCRIPT:-Rscript}"; OFF="unshare -rn"
  OUT=$M/results/u3ru5_integrity
  ARGS=(--u3r "$U3R" --source-status "$ROOT/01_sources/SOURCE_STATUS.tsv")
  /usr/bin/time -f "[time] producer elapsed %e s, MaxRSS %M KB, exit %x" $OFF "$RS" "$HERE/compare_u3ru5_integrity.R" --candidates-dir "$M/results/candidates" "${ARGS[@]}" --outdir "$OUT" 2>&1 | grep -E "^\[|time|Error"
  cat "$OUT/checksums.sha256"
)

[[ $# -gt 0 ]] || set -- "${STEPS[@]}"
for s in "$@"; do
  declare -F "step_$s" > /dev/null || { echo "unknown step: $s (steps: ${STEPS[*]})" >&2; exit 1; }
  echo "== 07_sparcs: $s"; "step_$s"
done
