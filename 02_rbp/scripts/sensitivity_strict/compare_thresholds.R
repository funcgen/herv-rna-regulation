#!/usr/bin/env Rscript
# =============================================================================
# 02_rbp / sensitivity_strict : threshold comparison and sensitivity summary
# =============================================================================
# Stage 13. Consolidates the historical scripts 11, 12 and 13:
#
#   11_compare_fimo_thresholds.R              -> sections 2 and 3 below
#   12_diagnose_fimo_1e-5_profile_stability.R -> sections 4 and 5
#   13_make_fimo_threshold_sensitivity_summary.R -> section 6
#
# Those three were sequential, not duplicated: 11 builds the comparison, 12
# diagnoses sparsity from 11's output, 13 repackages 11 and 12 and adds two
# hard-coded narrative tables. They are consolidated here because they form one
# task, and kept as clearly separated sections because they are distinct steps.
#
# READ-ONLY with respect to both analytical branches. Nothing upstream is
# recomputed: this stage consumes the accepted Stage 9/10 primary tables and the
# accepted Stage 12 strict tables, and compares already-adjusted values. In
# particular p_adj_global is NEVER recalculated -- both branches carry the same
# 9,672-test BH universe and that is what makes them comparable at all.
#
# Usage:
#   compare_thresholds.R --primary-subfamily D --strict-subfamily D \
#                        --primary-global D --strict-global D --out-root D
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  hit <- which(args == flag)
  if (length(hit) == 0) stop("missing required argument: ", flag)
  args[[hit + 1L]]
}

PRIMARY_SUBFAMILY <- get_arg("--primary-subfamily")
STRICT_SUBFAMILY  <- get_arg("--strict-subfamily")
PRIMARY_GLOBAL    <- get_arg("--primary-global")
STRICT_GLOBAL     <- get_arg("--strict-global")
OUT_ROOT          <- get_arg("--out-root")

COMPARISON_DIR <- file.path(OUT_ROOT, "threshold_comparison")
SUMMARY_DIR    <- file.path(OUT_ROOT, "sensitivity_summary")
dir.create(COMPARISON_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUMMARY_DIR, recursive = TRUE, showWarnings = FALSE)

families <- c("HERVH", "HERVK")

safe_cor <- function(x, y, method = "pearson") {
  valid <- is.finite(x) & is.finite(y)
  if (sum(valid) < 3L) return(NA_real_)
  cor(x[valid], y[valid], method = method)
}

collapse_rbps <- function(x) {
  x <- sort(unique(x))
  if (length(x) == 0L) return("")
  paste(x, collapse = ", ")
}

# =============================================================================
# 1. Whole-dataset global comparison            (historical script 11, part 1)
# =============================================================================
read_global_null <- function(path, threshold_label) {
  read_tsv(file.path(path, "global_null_summary.tsv"), show_col_types = FALSE) %>%
    mutate(threshold = threshold_label)
}

global_comparison <- bind_rows(
    read_global_null(PRIMARY_GLOBAL, "1e-4"),
    read_global_null(STRICT_GLOBAL, "1e-5")
  ) %>%
  filter(metric %in% c("loci_with_hit", "pct_with_hit", "total_hits",
                       "global_hits_per_kb")) %>%
  select(threshold, metric, observed, null_mean, observed_minus_null,
         observed_to_null_ratio, observed_above_all_nulls) %>%
  arrange(factor(threshold, levels = c("1e-4", "1e-5")),
          match(metric, c("loci_with_hit", "pct_with_hit", "total_hits",
                          "global_hits_per_kb")))

write_tsv(global_comparison,
          file.path(COMPARISON_DIR, "global_threshold_comparison.tsv"))

# =============================================================================
# 2. HERVH/HERVK family burden comparison       (historical script 11, part 2)
# =============================================================================
read_family_summary <- function(path, threshold_label) {
  read_tsv(file.path(path, "rbp_subfamily_length_gc_density_summary.tsv"),
           show_col_types = FALSE) %>%
    filter(subfamily %in% families) %>%
    mutate(threshold = threshold_label)
}

family_global_comparison <- bind_rows(
    read_family_summary(PRIMARY_SUBFAMILY, "1e-4"),
    read_family_summary(STRICT_SUBFAMILY, "1e-5")
  ) %>%
  select(threshold, subfamily, n_loci, total_length_kb, observed_total_hits,
         mean_shuffled_total_hits, observed_hits_per_kb,
         mean_shuffled_hits_per_kb, global_density_ratio,
         median_observed_locus_density, median_shuffled_locus_density,
         median_distinct_rbps_per_locus) %>%
  arrange(subfamily, factor(threshold, levels = c("1e-4", "1e-5")))

write_tsv(family_global_comparison,
          file.path(COMPARISON_DIR, "HERVH_HERVK_global_threshold_comparison.tsv"))

# =============================================================================
# 3. Per-RBP comparison                         (historical script 11, part 3)
# =============================================================================
# full_join, NOT inner_join: every RBP present in either branch is kept, and an
# RBP absent from one becomes status "missing" rather than vanishing.
read_enrichment <- function(path, suffix) {
  x <- read_tsv(file.path(path, "rbp_HERVH_HERVK_composition_controlled_enrichment.tsv"),
                show_col_types = FALSE) %>%
    filter(subfamily %in% families) %>%
    select(subfamily, rbp, observed_total_hits, mean_shuffled_total_hits,
           observed_hits_per_kb, mean_shuffled_hits_per_kb,
           log2_observed_expected, p_value, p_adj_global, status)
  names(x)[-(1:2)] <- paste0(names(x)[-(1:2)], "_", suffix)
  x
}

comparison <- full_join(
    read_enrichment(PRIMARY_SUBFAMILY, "1e4"),
    read_enrichment(STRICT_SUBFAMILY, "1e5"),
    by = c("subfamily", "rbp")
  ) %>%
  mutate(
    status_1e4 = coalesce(status_1e4, "missing"),
    status_1e5 = coalesce(status_1e5, "missing"),
    # DIRECTION is the sign of the effect, independent of significance.
    direction_1e4 = case_when(
      log2_observed_expected_1e4 > 0 ~ "positive",
      log2_observed_expected_1e4 < 0 ~ "negative",
      log2_observed_expected_1e4 == 0 ~ "zero",
      TRUE ~ "missing"),
    direction_1e5 = case_when(
      log2_observed_expected_1e5 > 0 ~ "positive",
      log2_observed_expected_1e5 < 0 ~ "negative",
      log2_observed_expected_1e5 == 0 ~ "zero",
      TRUE ~ "missing"),
    same_effect_direction = direction_1e4 == direction_1e5,
    effect_difference_1e5_minus_1e4 =
      log2_observed_expected_1e5 - log2_observed_expected_1e4
  ) %>%
  arrange(subfamily, rbp)

write_tsv(comparison,
          file.path(COMPARISON_DIR, "HERVH_HERVK_per_RBP_threshold_comparison.tsv"))

effect_concordance <- comparison %>%
  group_by(subfamily) %>%
  summarise(
    n_rbps = n(),
    pearson_effect_correlation = safe_cor(log2_observed_expected_1e4,
                                          log2_observed_expected_1e5, "pearson"),
    spearman_effect_correlation = safe_cor(log2_observed_expected_1e4,
                                           log2_observed_expected_1e5, "spearman"),
    n_same_effect_direction = sum(same_effect_direction, na.rm = TRUE),
    pct_same_effect_direction = 100 * mean(same_effect_direction, na.rm = TRUE),
    median_absolute_effect_change = median(abs(effect_difference_1e5_minus_1e4),
                                           na.rm = TRUE),
    .groups = "drop")

write_tsv(effect_concordance,
          file.path(COMPARISON_DIR, "HERVH_HERVK_effect_concordance.tsv"))

status_counts <- comparison %>%
  select(subfamily, rbp, status_1e4, status_1e5) %>%
  pivot_longer(cols = c(status_1e4, status_1e5),
               names_to = "threshold", values_to = "status") %>%
  mutate(threshold = recode(threshold, status_1e4 = "1e-4", status_1e5 = "1e-5")) %>%
  count(subfamily, threshold, status, name = "n_rbps") %>%
  arrange(subfamily, factor(threshold, levels = c("1e-4", "1e-5")), status)

write_tsv(status_counts,
          file.path(COMPARISON_DIR, "HERVH_HERVK_status_counts_by_threshold.tsv"))

status_transitions <- comparison %>%
  count(subfamily, status_1e4, status_1e5, name = "n_rbps") %>%
  arrange(subfamily, status_1e4, status_1e5)

write_tsv(status_transitions,
          file.path(COMPARISON_DIR, "HERVH_HERVK_status_transitions.tsv"))

set_summary <- comparison %>%
  group_by(subfamily) %>%
  summarise(
    enriched_1e4 = sum(status_1e4 == "enriched"),
    enriched_1e5 = sum(status_1e5 == "enriched"),
    enriched_retained = sum(status_1e4 == "enriched" & status_1e5 == "enriched"),
    enriched_lost = sum(status_1e4 == "enriched" & status_1e5 != "enriched"),
    enriched_gained = sum(status_1e4 != "enriched" & status_1e5 == "enriched"),
    depleted_1e4 = sum(status_1e4 == "depleted"),
    depleted_1e5 = sum(status_1e5 == "depleted"),
    depleted_retained = sum(status_1e4 == "depleted" & status_1e5 == "depleted"),
    depleted_lost = sum(status_1e4 == "depleted" & status_1e5 != "depleted"),
    depleted_gained = sum(status_1e4 != "depleted" & status_1e5 == "depleted"),
    .groups = "drop")

write_tsv(set_summary,
          file.path(COMPARISON_DIR, "HERVH_HERVK_significant_set_summary.tsv"))

rbp_set_changes <- comparison %>%
  group_by(subfamily) %>%
  summarise(
    retained_enriched = collapse_rbps(rbp[status_1e4 == "enriched" & status_1e5 == "enriched"]),
    lost_enriched     = collapse_rbps(rbp[status_1e4 == "enriched" & status_1e5 != "enriched"]),
    gained_enriched   = collapse_rbps(rbp[status_1e4 != "enriched" & status_1e5 == "enriched"]),
    retained_depleted = collapse_rbps(rbp[status_1e4 == "depleted" & status_1e5 == "depleted"]),
    lost_depleted     = collapse_rbps(rbp[status_1e4 == "depleted" & status_1e5 != "depleted"]),
    gained_depleted   = collapse_rbps(rbp[status_1e4 != "depleted" & status_1e5 == "depleted"]),
    .groups = "drop")

write_tsv(rbp_set_changes,
          file.path(COMPARISON_DIR, "HERVH_HERVK_RBP_set_changes.tsv"))

# =============================================================================
# 4. Strict-threshold sparsity diagnostics           (historical script 12)
# =============================================================================
# EVIDENCE is observed strict hits plus the SUMMED strict null across all five
# replicates (5 x the per-replicate mean). It is a count of how much data the
# strict comparison rests on, not a statistic.
# Re-read the per-RBP table from disk rather than reusing the in-memory tibble.
# That intermediate TSV was the documented interface between historical scripts
# 11 and 12, and the text round-trip is part of the values script 12 saw: reusing
# the in-memory doubles changes one p-value by a single ULP
# (HERVK/SRSF2 p_adj_global_1e4, 2.2935436106137897e-8 vs ...94e-8). Scientifically
# nil, but reading the file keeps the historical data flow exact.
x <- read_tsv(file.path(COMPARISON_DIR, "HERVH_HERVK_per_RBP_threshold_comparison.tsv"),
              show_col_types = FALSE) %>%
  mutate(
    strict_total_evidence =
      observed_total_hits_1e5 + 5 * mean_shuffled_total_hits_1e5,
    primary_significant = status_1e4 %in% c("enriched", "depleted"),
    primary_expected_direction = case_when(
      status_1e4 == "enriched" ~ "positive",
      status_1e4 == "depleted" ~ "negative",
      TRUE ~ "not_primary_significant"),
    strict_effect_direction = case_when(
      log2_observed_expected_1e5 > 0 ~ "positive",
      log2_observed_expected_1e5 < 0 ~ "negative",
      log2_observed_expected_1e5 == 0 ~ "zero",
      TRUE ~ "missing"),
    # Direction RETAINED / REVERSED are about the sign of the effect among
    # primary-significant RBPs. Losing significance is NOT a reversal, and a
    # strict effect of exactly zero is counted separately as neither.
    primary_direction_retained =
      primary_significant & primary_expected_direction == strict_effect_direction,
    primary_direction_reversed =
      primary_significant &
      strict_effect_direction %in% c("positive", "negative") &
      primary_expected_direction != strict_effect_direction
  )

sparsity <- x %>%
  group_by(subfamily) %>%
  summarise(
    n_rbps = n(),
    observed_zero = sum(observed_total_hits_1e5 == 0, na.rm = TRUE),
    shuffled_mean_zero = sum(mean_shuffled_total_hits_1e5 == 0, na.rm = TRUE),
    observed_and_shuffled_zero = sum(observed_total_hits_1e5 == 0 &
                                     mean_shuffled_total_hits_1e5 == 0, na.rm = TRUE),
    evidence_at_least_5 = sum(strict_total_evidence >= 5, na.rm = TRUE),
    evidence_at_least_10 = sum(strict_total_evidence >= 10, na.rm = TRUE),
    evidence_at_least_20 = sum(strict_total_evidence >= 20, na.rm = TRUE),
    .groups = "drop")

primary_direction <- x %>%
  filter(primary_significant) %>%
  count(subfamily, status_1e4, strict_effect_direction, name = "n_rbps") %>%
  arrange(subfamily, status_1e4, strict_effect_direction)

primary_direction_summary <- x %>%
  filter(primary_significant) %>%
  group_by(subfamily) %>%
  summarise(
    n_primary_significant = n(),
    direction_retained = sum(primary_direction_retained, na.rm = TRUE),
    direction_reversed = sum(primary_direction_reversed, na.rm = TRUE),
    strict_zero_effect = sum(strict_effect_direction == "zero", na.rm = TRUE),
    pct_direction_retained = 100 * mean(primary_direction_retained, na.rm = TRUE),
    .groups = "drop")

supported_direction <- bind_rows(lapply(c(5, 10, 20), function(minimum_evidence) {
  x %>%
    filter(primary_significant, strict_total_evidence >= minimum_evidence) %>%
    group_by(subfamily) %>%
    summarise(
      minimum_strict_evidence = minimum_evidence,
      n_primary_significant_supported = n(),
      direction_retained = sum(primary_direction_retained, na.rm = TRUE),
      direction_reversed = sum(primary_direction_reversed, na.rm = TRUE),
      pct_direction_retained = 100 * mean(primary_direction_retained, na.rm = TRUE),
      .groups = "drop")
})) %>% arrange(subfamily, minimum_strict_evidence)

correlation_by_support <- bind_rows(lapply(c(0, 5, 10, 20), function(minimum_evidence) {
  x %>%
    filter(strict_total_evidence >= minimum_evidence) %>%
    group_by(subfamily) %>%
    summarise(
      minimum_strict_evidence = minimum_evidence,
      n_rbps = n(),
      pearson = safe_cor(log2_observed_expected_1e4, log2_observed_expected_1e5, "pearson"),
      spearman = safe_cor(log2_observed_expected_1e4, log2_observed_expected_1e5, "spearman"),
      same_direction = sum(direction_1e4 == direction_1e5, na.rm = TRUE),
      pct_same_direction = 100 * mean(direction_1e4 == direction_1e5, na.rm = TRUE),
      .groups = "drop")
})) %>% arrange(subfamily, minimum_strict_evidence)

# A genuine reversal requires significance in BOTH branches with opposite calls.
status_reversals <- x %>%
  filter(status_1e4 %in% c("enriched", "depleted"),
         status_1e5 %in% c("enriched", "depleted"),
         status_1e4 != status_1e5) %>%
  select(subfamily, rbp, status_1e4, status_1e5,
         observed_total_hits_1e4, mean_shuffled_total_hits_1e4,
         log2_observed_expected_1e4, p_adj_global_1e4,
         observed_total_hits_1e5, mean_shuffled_total_hits_1e5,
         strict_total_evidence, log2_observed_expected_1e5, p_adj_global_1e5) %>%
  arrange(subfamily, rbp)

write_tsv(sparsity, file.path(COMPARISON_DIR, "strict_threshold_sparsity.tsv"))
write_tsv(primary_direction,
          file.path(COMPARISON_DIR, "primary_significant_direction_at_1e-5.tsv"))
write_tsv(primary_direction_summary,
          file.path(COMPARISON_DIR, "primary_significant_direction_summary.tsv"))
write_tsv(supported_direction,
          file.path(COMPARISON_DIR, "primary_direction_by_strict_support.tsv"))
write_tsv(correlation_by_support,
          file.path(COMPARISON_DIR, "effect_correlation_by_strict_support.tsv"))
write_tsv(status_reversals,
          file.path(COMPARISON_DIR, "significant_status_reversals.tsv"))

# =============================================================================
# 5. Sensitivity summary                             (historical script 13)
# =============================================================================
# Repackaging plus two narrative tables. The narrative tables are HARD-CODED
# prose in the historical script, not computed; they are reproduced verbatim
# because they carry the project's stated interpretive limits.
global_wide <- global_comparison %>%
  filter(metric %in% c("pct_with_hit", "total_hits", "global_hits_per_kb")) %>%
  select(threshold, metric, observed, null_mean, observed_to_null_ratio)

family_compact <- family_global_comparison %>%
  select(threshold, subfamily, observed_total_hits, mean_shuffled_total_hits,
         observed_hits_per_kb, mean_shuffled_hits_per_kb, global_density_ratio)

functional_summary <- tribble(
  ~threshold, ~analysis, ~result,
  "1e-4", "HERVH enriched-RBP GO ORA", "No significant GO BP term after BH correction",
  "1e-4", "HERVK enriched-RBP GO ORA", "No significant GO BP term after BH correction",
  "1e-4", "HERVH-versus-HERVK ranked GO BP",
  "mRNA splicing, via spliceosome was HERVK-biased; adjusted P = 0.004",
  "1e-5", "HERVH enriched-RBP GO ORA",
  paste0("GO:0051239 regulation of multicellular organismal process; ",
         "adjusted P = 0.025; broad propagated ancestor term supported by ",
         "heterogeneous descendant annotations"),
  "1e-5", "HERVK enriched-RBP GO ORA", "No significant GO BP term after BH correction",
  "1e-5", "HERVH-versus-HERVK ranked GO BP", "No significant GO BP term after BH correction")

interpretation <- tribble(
  ~topic, ~conclusion,
  "Global matched-null result",
  paste0("Robust. Observed motif burden exceeded the dinucleotide-matched ",
         "null at both FIMO thresholds."),
  "HERVH family-level burden",
  paste0("Robust. HERVH showed composition-controlled excess at both ",
         "thresholds."),
  "HERVK family-level burden",
  paste0("Threshold-dependent. HERVK was near the null at 1e-4 but showed ",
         "an excess among the strongest 1e-5 motif calls."),
  "Individual RBP profiles",
  paste0("Threshold-dependent and highly sparse at 1e-5; individual ",
         "enrichment calls should not be presented as threshold-invariant."),
  "Functional categories",
  paste0("Threshold-dependent. No functional claim was reproduced ",
         "consistently across both thresholds."),
  "Primary threshold",
  paste0("Retain 1e-4 as the primary threshold and report 1e-5 as a ",
         "sensitivity analysis."))

write_tsv(global_wide, file.path(SUMMARY_DIR, "fimo_threshold_global_summary.tsv"))
write_tsv(family_compact, file.path(SUMMARY_DIR, "fimo_threshold_HERVH_HERVK_summary.tsv"))
write_tsv(set_summary, file.path(SUMMARY_DIR, "fimo_threshold_significant_RBP_counts.tsv"))
write_tsv(rbp_set_changes, file.path(SUMMARY_DIR, "fimo_threshold_RBP_set_changes.tsv"))
write_tsv(sparsity, file.path(SUMMARY_DIR, "fimo_1e-5_sparsity_summary.tsv"))
write_tsv(functional_summary, file.path(SUMMARY_DIR, "fimo_threshold_functional_summary.tsv"))
write_tsv(interpretation, file.path(SUMMARY_DIR, "fimo_threshold_interpretation.tsv"))

# The prose report is likewise hard-coded narrative in the historical script.
# Reproduced verbatim: it is the project's own statement of what these numbers
# do and do not support.
report <- c(
  "FIMO THRESHOLD SENSITIVITY ANALYSIS",
  "===================================",
  "",
  paste0("At P <= 1e-4, the global observed/null motif-burden ratio was ",
         "1.162. At P <= 1e-5, it remained elevated at 1.120."),
  "",
  paste0("HERVH remained enriched relative to its dinucleotide-preserving ",
         "null at both thresholds: O/E = 1.216 at 1e-4 and 2.171 at 1e-5."),
  "",
  paste0("HERVK was close to its global shuffled expectation at 1e-4 ",
         "(O/E = 1.005), whereas its strongest 1e-5 calls showed an excess ",
         "(O/E = 1.296)."),
  "",
  paste0("The strict analysis was highly sparse: 80 of 93 RBP profiles had ",
         "zero observed and shuffled hits in each family. Consequently, ",
         "individual RBP enrichment calls and functional annotations were ",
         "threshold-dependent."),
  "",
  paste0("The single significant HERVH GO term at 1e-5 was the broad ",
         "ancestor category GO:0051239, whose membership was propagated ",
         "from heterogeneous descendant annotations. The ranked GO analysis ",
         "identified no significant term at 1e-5."),
  "",
  paste0("Conclusion: the global composition-controlled motif excess is ",
         "robust, whereas detailed RBP-specific and pathway-level profiles ",
         "should be interpreted cautiously and in relation to the selected ",
         "FIMO threshold.")
)

writeLines(report, file.path(SUMMARY_DIR, "fimo_threshold_sensitivity_report.txt"))

writeLines(capture.output(sessionInfo()), file.path(OUT_ROOT, "sessionInfo.txt"))

cat("\n=== STRICT-THRESHOLD SPARSITY ===\n"); print(sparsity, n = Inf, width = Inf)
cat("\n=== STATUS TRANSITIONS ===\n"); print(status_transitions, n = Inf)
cat("\n=== SIGNIFICANT SET SUMMARY ===\n"); print(set_summary, n = Inf, width = Inf)
cat("\n=== DIRECTION SUMMARY ===\n"); print(primary_direction_summary, n = Inf, width = Inf)
cat("\nWritten to:\n", COMPARISON_DIR, "\n", SUMMARY_DIR, "\n")
