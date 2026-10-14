#!/usr/bin/env Rscript
# =============================================================================
# 02_rbp / enrichment : global observed-versus-null summary (Stage 8)
# =============================================================================
# Refactor of scripts/rbp/08_summarise_global_dinuc_null.R. The statistics are
# reproduced exactly; only paths, input resolution and validation changed.
#
#   observed FIMO (Stage 3) + five shuffled FIMO scans (Stage 7)
#   + canonical source FASTA (sequence lengths)
#     -> results/enrichment/global_metrics_by_dataset.tsv
#        results/enrichment/global_null_summary.tsv
#        results/enrichment/per_locus_observed_vs_mean_null.tsv
#        results/enrichment/paired_locus_tests.tsv
#
# SEMANTICS PRESERVED FROM THE HISTORICAL SCRIPT -- each of these changes the
# numbers if altered:
#
#   * FULL LOCUS UNIVERSE. A FIMO file lists only loci with >=1 retained call.
#     Every dataset is left_joined onto the complete FASTA locus table and
#     missing counts are replaced with 0, so all 53,176 loci are present in
#     every dataset and the paired test is over the full universe. An inner
#     join would silently shrink the denominator.
#   * LENGTH DENOMINATOR = ACTUAL SEQUENCE LENGTH, summed from the FASTA
#     records themselves -- NOT (end - start) parsed from the sequence name.
#     Those differ by one base per locus (Stage 4), so the choice is material.
#   * PER-LOCUS DENSITY = 1000 * n_hits / length_bp; GLOBAL DENSITY =
#     1000 * sum(n_hits) / sum(length_bp). The global value is a ratio of sums,
#     not a mean of per-locus densities.
#   * NULL AGGREGATION IS PER METRIC, ACROSS DATASETS. null_mean for a metric is
#     the mean of that metric's five per-replicate values. So the null median
#     density is the MEAN OF THE FIVE REPLICATE MEDIANS -- not the median of
#     per-locus means. The per-locus table separately carries mean_null_hits,
#     which is the per-locus mean across replicates and is what the paired test
#     uses. The two are different quantities and both are historical.
#   * PAIRED WILCOXON, two-sided, exact = FALSE, over all 53,176 loci.
#
# Usage:
#   summarise_global_null.R --manifest PATH --fasta PATH --out-dir PATH
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
})

# ------------------------------------------------------------------ arguments
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, required = TRUE) {
  hit <- which(args == flag)
  if (length(hit) == 0) {
    if (required) stop("missing required argument: ", flag)
    return(NULL)
  }
  args[[hit + 1L]]
}

manifest_path <- get_arg("--manifest")
fasta_path    <- get_arg("--fasta")
output_dir    <- get_arg("--out-dir")
# Explicit rather than derived from the manifest path: the manifest sits two
# directories below the module root, and deriving it by repeated dirname() is
# exactly the kind of silent off-by-one that fails late.
module_dir    <- get_arg("--module-dir")

for (p in c(manifest_path, fasta_path)) {
  if (!file.exists(p)) stop("input not found: ", p)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------- datasets from the manifest
# The Stage 7 manifest is the canonical dataset mapping; inputs are never
# rediscovered by globbing.
manifest <- read_tsv(manifest_path, show_col_types = FALSE, progress = FALSE)

required_cols <- c("dataset", "fimo_output", "n_hits", "status")
missing_cols <- setdiff(required_cols, names(manifest))
if (length(missing_cols) > 0) {
  stop("manifest is missing columns: ", paste(missing_cols, collapse = ", "))
}
if (any(manifest$status != "PASS")) {
  stop("manifest contains datasets whose status is not PASS")
}

datasets <- c("observed", sprintf("shuffle%02d", 1:5))
if (!setequal(manifest$dataset, datasets)) {
  stop("manifest datasets are ", paste(manifest$dataset, collapse = ", "),
       "; expected ", paste(datasets, collapse = ", "))
}

if (!dir.exists(module_dir)) stop("module dir not found: ", module_dir)
fimo_paths <- setNames(
  file.path(module_dir, manifest$fimo_output[match(datasets, manifest$dataset)]),
  datasets
)
expected_hits <- setNames(
  as.numeric(manifest$n_hits[match(datasets, manifest$dataset)]),
  datasets
)

missing_files <- fimo_paths[!file.exists(fimo_paths)]
if (length(missing_files) > 0) {
  stop("Missing FIMO files:\n", paste(missing_files, collapse = "\n"))
}

# ------------------------------------------------------------- FASTA lengths
read_fasta_lengths <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]

  header_positions <- which(startsWith(lines, ">"))
  if (length(header_positions) == 0) stop("No FASTA records found in: ", path)

  record_ends <- c(header_positions[-1] - 1L, length(lines))
  ids <- str_replace(lines[header_positions], "^>(\\S+).*$", "\\1")

  lengths <- map2_int(
    header_positions, record_ends,
    function(start_position, end_position) {
      sum(nchar(lines[seq.int(start_position + 1L, end_position)]))
    }
  )

  output <- tibble(locus_id = ids, length_bp = lengths)
  if (anyDuplicated(output$locus_id)) stop("Duplicated FASTA identifiers detected")
  output
}

sequence_lengths <- read_fasta_lengths(fasta_path)
message("FASTA loci: ", nrow(sequence_lengths))
message("Total bases: ", sum(sequence_lengths$length_bp))

# ------------------------------------ one FIMO file -> complete per-locus table
read_fimo_per_locus <- function(path, dataset, sequence_lengths, expect_hits) {
  fimo <- read_tsv(path, comment = "#", show_col_types = FALSE, progress = FALSE)

  required <- c("motif_id", "motif_alt_id", "sequence_name", "p-value")
  missing <- setdiff(required, names(fimo))
  if (length(missing) > 0) {
    stop(dataset, " is missing columns: ", paste(missing, collapse = ", "))
  }

  # Gate on the accepted Stage 3 / Stage 7 hit count before any statistics.
  if (!is.na(expect_hits) && nrow(fimo) != expect_hits) {
    stop(dataset, ": read ", nrow(fimo), " rows, manifest expects ", expect_hits)
  }

  hits <- fimo %>%
    transmute(
      locus_id = str_remove(sequence_name, "__dinuc_shuf\\d+$"),
      rbp = case_when(
        !is.na(motif_alt_id) & motif_alt_id != "" ~ str_extract(motif_alt_id, "^[^|]+"),
        TRUE ~ motif_id
      ),
      p_value = as.numeric(`p-value`)
    )

  unknown_ids <- setdiff(unique(hits$locus_id), sequence_lengths$locus_id)
  if (length(unknown_ids) > 0) {
    stop(dataset, ": FIMO contains identifiers absent from original FASTA. Example: ",
         unknown_ids[[1]])
  }

  hit_summary <- hits %>%
    group_by(locus_id) %>%
    summarise(n_hits = n(), n_unique_rbps = n_distinct(rbp), .groups = "drop")

  # left_join + replace_na restores the full 53,176-locus universe.
  sequence_lengths %>%
    left_join(hit_summary, by = "locus_id") %>%
    mutate(
      n_hits = replace_na(n_hits, 0L),
      n_unique_rbps = replace_na(n_unique_rbps, 0L),
      has_hit = n_hits > 0,
      hits_per_kb = 1000 * n_hits / length_bp,
      unique_rbps_per_kb = 1000 * n_unique_rbps / length_bp,
      dataset = dataset
    )
}

per_locus_list <- imap(
  fimo_paths,
  ~ read_fimo_per_locus(
      path = .x, dataset = .y,
      sequence_lengths = sequence_lengths,
      expect_hits = expected_hits[[.y]]
    )
)
per_locus <- bind_rows(per_locus_list)

# ------------------------------------------------------------- global metrics
global_metrics <- per_locus %>%
  group_by(dataset) %>%
  summarise(
    n_loci = n(),
    loci_with_hit = sum(has_hit),
    pct_with_hit = 100 * mean(has_hit),
    total_hits = sum(n_hits),
    global_hits_per_kb = 1000 * sum(n_hits) / sum(length_bp),
    median_hits_per_locus = median(n_hits),
    median_hits_per_kb = median(hits_per_kb),
    median_unique_rbps = median(n_unique_rbps),
    median_unique_rbps_per_kb = median(unique_rbps_per_kb),
    .groups = "drop"
  ) %>%
  mutate(dataset = factor(dataset, levels = datasets)) %>%
  arrange(dataset) %>%
  mutate(dataset = as.character(dataset))

write_tsv(global_metrics, file.path(output_dir, "global_metrics_by_dataset.tsv"))

# ------------------------------- observed versus the five-replicate null spread
global_long <- global_metrics %>%
  pivot_longer(cols = -dataset, names_to = "metric", values_to = "value")

observed_long <- global_long %>%
  filter(dataset == "observed") %>%
  select(metric, observed = value)

null_summary <- global_long %>%
  filter(dataset != "observed") %>%
  group_by(metric) %>%
  summarise(
    null_mean = mean(value), null_sd = sd(value),
    null_min = min(value), null_max = max(value),
    .groups = "drop"
  )

global_null_summary <- observed_long %>%
  left_join(null_summary, by = "metric") %>%
  mutate(
    observed_minus_null = observed - null_mean,
    observed_to_null_ratio = if_else(null_mean == 0, NA_real_, observed / null_mean),
    observed_above_all_nulls = observed > null_max
  )

write_tsv(global_null_summary, file.path(output_dir, "global_null_summary.tsv"))

# ------------------------------------------- paired locus-level comparison
observed_loci <- per_locus %>%
  filter(dataset == "observed") %>%
  select(
    locus_id, length_bp,
    observed_hits = n_hits,
    observed_hits_per_kb = hits_per_kb,
    observed_unique_rbps = n_unique_rbps,
    observed_has_hit = has_hit
  )

mean_null_loci <- per_locus %>%
  filter(dataset != "observed") %>%
  group_by(locus_id) %>%
  summarise(
    mean_null_hits = mean(n_hits),
    mean_null_hits_per_kb = mean(hits_per_kb),
    mean_null_unique_rbps = mean(n_unique_rbps),
    proportion_null_with_hit = mean(has_hit),
    .groups = "drop"
  )

paired_loci <- observed_loci %>%
  left_join(mean_null_loci, by = "locus_id") %>%
  mutate(
    hits_per_kb_difference = observed_hits_per_kb - mean_null_hits_per_kb,
    unique_rbps_difference = observed_unique_rbps - mean_null_unique_rbps
  )

write_tsv(paired_loci, file.path(output_dir, "per_locus_observed_vs_mean_null.tsv"))

paired_tests <- tibble(
  metric = c("hits_per_locus", "hits_per_kb", "unique_rbps_per_locus"),
  alternative = "two.sided",
  p_value = c(
    wilcox.test(paired_loci$observed_hits, paired_loci$mean_null_hits,
                paired = TRUE, alternative = "two.sided", exact = FALSE)$p.value,
    wilcox.test(paired_loci$observed_hits_per_kb, paired_loci$mean_null_hits_per_kb,
                paired = TRUE, alternative = "two.sided", exact = FALSE)$p.value,
    wilcox.test(paired_loci$observed_unique_rbps, paired_loci$mean_null_unique_rbps,
                paired = TRUE, alternative = "two.sided", exact = FALSE)$p.value
  )
)

write_tsv(paired_tests, file.path(output_dir, "paired_locus_tests.tsv"))

# The V statistic and n are not part of the historical output schema, so they are
# written to a SEPARATE file rather than added as columns. Reporting-only: they do
# not enter any manuscript value.
paired_statistics <- tibble(
  metric = c("hits_per_locus", "hits_per_kb", "unique_rbps_per_locus"),
  test = "wilcox.test(paired = TRUE, alternative = \"two.sided\", exact = FALSE)",
  n_pairs = nrow(paired_loci),
  statistic_V = c(
    wilcox.test(paired_loci$observed_hits, paired_loci$mean_null_hits,
                paired = TRUE, alternative = "two.sided", exact = FALSE)$statistic,
    wilcox.test(paired_loci$observed_hits_per_kb, paired_loci$mean_null_hits_per_kb,
                paired = TRUE, alternative = "two.sided", exact = FALSE)$statistic,
    wilcox.test(paired_loci$observed_unique_rbps, paired_loci$mean_null_unique_rbps,
                paired = TRUE, alternative = "two.sided", exact = FALSE)$statistic
  )
)
write_tsv(paired_statistics, file.path(output_dir, "paired_locus_test_statistics.tsv"))

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

# ------------------------------------------------------------------- console
cat("\n=== GLOBAL METRICS ===\n"); print(global_metrics, n = Inf)
cat("\n=== OBSERVED VERSUS NULL ===\n")
print(global_null_summary %>%
        filter(metric %in% c("pct_with_hit", "total_hits", "global_hits_per_kb",
                             "median_hits_per_kb", "median_unique_rbps")), n = Inf)
cat("\n=== PAIRED TESTS ===\n"); print(paired_tests, n = Inf)
cat("\nResults written to:\n", output_dir, "\n")
