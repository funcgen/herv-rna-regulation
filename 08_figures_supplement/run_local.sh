#!/usr/bin/env bash
# 08_figures_supplement: figure panels, Supplemental Figure 2 and Supplemental Tables 1-8 from the module results
# (effective orders 47-51, 54, 54.1; see docs/WORKFLOW.md). No scientific computation happens here.
# Each step is the producer stage of the validated driver named in its "from:" line, copied line for line; the
# validation, determinism, fixture and historical-comparison stages of those drivers are not part of this release.
# The HERVarium browser panels (Figures 1C, 2E, 3C), the Figure 5A schematic, Supplemental Figure 1 and the final
# multi-panel layouts are manual assets and are not produced here.
# Inputs passed to a driver as an option are read from a variable (default: the module location):
#   NCRNA_DOMAINS_RESULTS   04_domains_terminal_exons results directory (step figure2)
# Requires: the results of modules 02-07, NCRNA_RSCRIPT (optional), unshare (util-linux).
# Usage: 08_figures_supplement/run_local.sh [step ...]        (no argument: all steps in order)
set -euo pipefail
M="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(cd "$M/.." && pwd)"
STEPS=(figure1 figure2 figure3 figure4 figure5 supplement supplementary_figure2)

step_figure1() (
  # from: 08_figures_supplement/scripts/figure1/run_plot_figure1.sh
  RS="${NCRNA_RSCRIPT:-Rscript}"; OFF="unshare -rn"; P="$M/scripts/plot_figure1_panels.R"
  ARGS=(--rbp-dir "$ROOT/02_rbp" --params "$ROOT/config/parameters.yaml")
  /usr/bin/time -f "[time] plotting elapsed %e s, MaxRSS %M KB, exit %x" $OFF "$RS" "$P" "${ARGS[@]}" --sources-out "$M/panel_sources/figure1" --figures-out "$M/figures/figure1" 2>&1 | grep -E "^\[|time|Error"
  cat "$M/panel_sources/figure1/checksums.sha256" "$M/figures/figure1/checksums.sha256"
)

step_figure2() (
  # from: 08_figures_supplement/scripts/figure2/run_plot_figure2.sh
  DOM="${NCRNA_DOMAINS_RESULTS:-$ROOT/04_domains_terminal_exons/results}"
  RS="${NCRNA_RSCRIPT:-Rscript}"; OFF="unshare -rn"; P="$M/scripts/plot_figure2_panels.R"
  ARGS=(--context-dir "$ROOT/03_transcript_context" --domains-results "$DOM" --domains-manifest "$ROOT/04_domains_terminal_exons/RESULTS_MANIFEST.tsv")
  /usr/bin/time -f "[time] plotting elapsed %e s, MaxRSS %M KB, exit %x" $OFF "$RS" "$P" "${ARGS[@]}" --sources-out "$M/panel_sources/figure2" --figures-out "$M/figures/figure2" 2>&1 | grep -E "^\[|time|Error"
  cat "$M/panel_sources/figure2/checksums.sha256" "$M/figures/figure2/checksums.sha256"
)

step_figure3() (
  # from: 08_figures_supplement/scripts/figure3/run_plot_figure3.sh
  RS="${NCRNA_RSCRIPT:-Rscript}"; OFF="unshare -rn"; P="$M/scripts/plot_figure3_panels.R"
  ARGS=(--lncrna-dir "$ROOT/05_lncrna" --params "$ROOT/config/parameters.yaml")
  /usr/bin/time -f "[time] plotting elapsed %e s, MaxRSS %M KB, exit %x" $OFF "$RS" "$P" "${ARGS[@]}" --sources-out "$M/panel_sources/figure3" --figures-out "$M/figures/figure3" 2>&1 | grep -E "^\[|time|Error"
  cat "$M/panel_sources/figure3/checksums.sha256" "$M/figures/figure3/checksums.sha256"
)

step_figure4() (
  # from: 08_figures_supplement/scripts/figure4/run_plot_figure4.sh
  WORK="${TMPDIR:-/tmp}/figure4_work.$$"
  RS="${NCRNA_RSCRIPT:-Rscript}"; OFF="unshare -rn"; P="$M/scripts/plot_figure4.R"
  mkdir -p "$WORK"
  cd "$WORK"
  ARGS=(--function-dir "$ROOT/06_lncrna_function" --lncrna-dir "$ROOT/05_lncrna")
  /usr/bin/time -f "[time] plotting elapsed %e s, MaxRSS %M KB, exit %x" $OFF "$RS" "$P" "${ARGS[@]}" --sources-out "$M/panel_sources/figure4" --figures-out "$M/figures/figure4" 2>&1 | grep -E "^\[|time|^Error"
  cat "$M/panel_sources/figure4/checksums.sha256" "$M/figures/figure4/checksums.sha256"
  cd /; rm -rf "$WORK"
)

step_figure5() (
  # from: 08_figures_supplement/scripts/figure5/run_plot_figure5.sh
  WORK="${TMPDIR:-/tmp}/figure5_work.$$"
  RS="${NCRNA_RSCRIPT:-Rscript}"; OFF="unshare -rn"; P="$M/scripts/plot_figure5_panels.R"
  mkdir -p "$WORK"
  cd "$WORK"
  ARGS=(--sparcs-dir "$ROOT/07_sparcs")
  /usr/bin/time -f "[time] plotting elapsed %e s, MaxRSS %M KB, exit %x" $OFF "$RS" "$P" "${ARGS[@]}" --sources-out "$M/panel_sources/figure5" --figures-out "$M/figures/figure5" 2>&1 | grep -E "^\[|time|^Error"
  cat "$M/panel_sources/figure5/checksums.sha256" "$M/figures/figure5/checksums.sha256"
  cd /; rm -rf "$WORK"
)

step_supplement() (
  # from: 08_figures_supplement/scripts/supplement/run_build_supplement.sh
  HERE="$M/scripts/supplement"
  WORK="${TMPDIR:-/tmp}/supplement_work.$$"
  RS="${NCRNA_RSCRIPT:-Rscript}"; OFF="unshare -rn"; P="$M/scripts/build_supplement.R"
  export TZ=UTC
  rm -rf "$WORK"; mkdir -p "$WORK"
  SECONDS=0; $OFF "$RS" "$P" --release "$ROOT" --out "$WORK/run1" --work "$WORK/tmp1"; echo "[time] run 1: ${SECONDS} s"
  rm -rf "$M/supplement/tables" "$M/supplement/archives" "$M/supplement/checksums.sha256" "$M/panel_sources/supplement"
  mkdir -p "$M/supplement" "$M/panel_sources"
  cp -a "$WORK/run1/supplement/." "$M/supplement/"; cp -a "$WORK/run1/panel_sources/supplement" "$M/panel_sources/"
  cd /; rm -rf "$WORK"
)

step_supplementary_figure2() (
  # from: 08_figures_supplement/scripts/supplementary_figures/run_supplementary_figures.sh
  HERE="$M/scripts/supplementary_figures"
  WORK="${TMPDIR:-/tmp}/supplementary_figures_work.$$"
  RS="${NCRNA_RSCRIPT:-Rscript}"; OFF="unshare -rn"
  export TZ=UTC
  rm -rf "${WORK:?}"; mkdir -p "$WORK"
  for r in 1 2; do $OFF "$RS" "$HERE/plot_supplementary_figure2.R" --mode canonical --release "$ROOT" --figures-out "$WORK/run$r/figures/supplementary" \
    --sources-out "$WORK/run$r/panel_sources/supplementary_figures" 2>&1 | grep -E "^\[OK\]|^Error" | sed "s/^/[run $r] /"; done
  diff -r "$WORK/run1" "$WORK/run2" > /dev/null && echo "[ok] determinism: run 2 byte-identical to run 1 ($(find "$WORK/run1" -type f | wc -l) files)" || { echo "[FAIL] runs differ"; diff -rq "$WORK/run1" "$WORK/run2"; exit 1; }
  for d in figures/supplementary panel_sources/supplementary_figures; do
    ( cd "$WORK/run1/$d" && find . -type f ! -name checksums.sha256 | LC_ALL=C sort | xargs sha256sum > checksums.sha256 )
    rm -rf "${M:?}/${d:?}"; mkdir -p "$M/$d"; cp -a "$WORK/run1/$d/." "$M/$d/"
  done
  cd /; rm -rf "${WORK:?}"
)

[[ $# -gt 0 ]] || set -- "${STEPS[@]}"
for s in "$@"; do
  declare -F "step_$s" > /dev/null || { echo "unknown step: $s (steps: ${STEPS[*]})" >&2; exit 1; }
  echo "== 08_figures_supplement: $s"; "step_$s"
done
