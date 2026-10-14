#!/usr/bin/env bash
# 02_rbp: local step (effective order 29, FMR1 locus screen; see docs/WORKFLOW.md). All other 02_rbp steps run on SLURM.
# The step is the producer stage of the validated driver named in its "from:" line, copied line for line; the
# validation, fixture and historical-comparison stages of that driver are not part of this release.
# Requires: the 02_rbp primary scan results, the 03_transcript_context integrated table, NCRNA_RSCRIPT (optional),
#           python3 with PyYAML.
# Usage: 02_rbp/run_local.sh [step ...]        (no argument: all steps in order)
set -euo pipefail
M="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(cd "$M/.." && pwd)"
STEPS=(fmr1_locus_screen)

step_fmr1_locus_screen() (
  # from: 02_rbp/scripts/fmr1_locus_screen/run_screen_fmr1_loci.sh
  HERE="$M/scripts/fmr1_locus_screen"
  RS="${NCRNA_RSCRIPT:-Rscript}"
  FIMO="$M/results/scanning/primary_fimo_parsed.tsv"
  ATLAS="$ROOT/03_transcript_context/results/integration/transcript_context.with_herv_id.tsv"
  OUT="$M/results/fmr1_locus_screen"; FIG="$M/figure_source_data/figure1b"
  read -r RBP SUB PCUT CMP MINH TOPN < <(python3 -c "import yaml;p=yaml.safe_load(open('$ROOT/config/parameters.yaml'))['rbp']['fmr1_locus_screen'];print(p['rbp'],p['subfamily'],p['pvalue_cut'],p['pvalue_comparison'],p['candidate_min_hits'],p['candidate_top_n'])")
  [[ "$CMP" == "less_equal" ]] || { echo "[error] only the inclusive comparison is implemented (D016)" >&2; exit 1; }
  echo "[params] rbp=$RBP subfamily=$SUB p<=$PCUT min_hits=$MINH top_n=$TOPN"
  /usr/bin/time -f "[time] producer elapsed %e s, MaxRSS %M KB, exit %x" \
    "$RS" "$HERE/screen_fmr1_loci.R" --fimo "$FIMO" --atlas "$ATLAS" --results-dir "$OUT" --figure-dir "$FIG" \
    --rbp "$RBP" --subfamily "$SUB" --p-cut "$PCUT" --min-hits "$MINH" --top-n "$TOPN"
)

[[ $# -gt 0 ]] || set -- "${STEPS[@]}"
for s in "$@"; do
  declare -F "step_$s" > /dev/null || { echo "unknown step: $s (steps: ${STEPS[*]})" >&2; exit 1; }
  echo "== 02_rbp: $s"; "step_$s"
done
