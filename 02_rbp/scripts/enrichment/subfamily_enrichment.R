#!/usr/bin/env Rscript
# =============================================================================
# 02_rbp / enrichment : subfamily composition-controlled enrichment (Stage 9)
# =============================================================================
# Refactor of scripts/rbp/09_subfamily_composition_controlled_enrichment.R.
# Every statistical definition is reproduced exactly; only paths, input
# resolution and input validation changed.
#
#   observed FIMO (Stage 3) + five null FIMO scans (Stage 7)
#   + canonical FASTA (length, GC, subfamily) + Stage 1 motif library (RBP universe)
#     -> results/enrichment/subfamily/*.tsv
#
# THE COMPOSITION CONTROL IS THE SHUFFLE ITSELF. The five null sequence sets
# preserve each locus's length, mononucleotide and dinucleotide composition and
# terminal bases, so no GC / length / k-mer covariate is modelled here. Adding
# one would double-correct.
#
# DEFINITIONS PRESERVED FROM THE HISTORICAL SCRIPT -- each changes published
# numbers if altered:
#
#   * SUBFAMILY = str_remove(locus_id, "-int_.*$").
#   * LENGTH = actual FASTA sequence length, not end - start from the name.
#   * FULL GRID. Each subfamily is expanded to every locus x every RBP in the
#     93-RBP universe, and absent combinations become 0. Zero-hit loci and
#     zero-hit RBPs are therefore present, and the Wilcoxon runs over all loci
#     of the subfamily, not only motif-positive ones.
#   * mean_shuffled_hits is the PER LOCUS x RBP mean across the five replicates,
#     so it is generally fractional and its subfamily sum usually is too.
#   * SUBFAMILY DENSITY = total hits / total subfamily length in kb -- a ratio of
#     sums, not a mean of per-locus densities.
#   * log2 O/E uses a 0.5 PSEUDOCOUNT on both sides, solely so the ratio stays
#     finite when a total is zero. It is applied to the effect size ONLY, never
#     to the counts, densities or the test.
#   * TEST: paired Wilcoxon over loci within a subfamily, observed vs mean
#     shuffled, two-sided, exact = FALSE, correct = TRUE. Returns NA for
#     subfamilies with fewer than MIN_LOCI_FOR_TEST loci or when every
#     difference is zero.
#   * FDR: BH across ALL subfamily x RBP tests globally (p_adj_global) AND, as a
#     separate column, BH within each subfamily (p_adj_within_subfamily).
#     `status` is decided by p_adj_global -- this is the manuscript's criterion.
#
# Usage:
#   subfamily_enrichment.R --module-dir D --manifest F --fasta F --motifs F --out-dir D
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
get_arg <- function(flag) {
  hit <- which(args == flag)
  if (length(hit) == 0) stop("missing required argument: ", flag)
  args[[hit + 1L]]
}

MODULE_DIR <- get_arg("--module-dir")
MANIFEST   <- get_arg("--manifest")
FASTA      <- get_arg("--fasta")
MOTIF_FILE <- get_arg("--motifs")
OUTDIR     <- get_arg("--out-dir")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Historical constants, unchanged.
WRITE_DENSE_LOCUS_RBP_TABLE <- TRUE
MIN_LOCI_FOR_TEST <- 5L
EFFECT_PSEUDOCOUNT <- 0.5
SIGNIFICANCE_ALPHA <- 0.05

# ------------------------------------------------- datasets from the manifest
manifest <- read_tsv(MANIFEST, show_col_types = FALSE, progress = FALSE)
if (any(manifest$status != "PASS")) stop("manifest has datasets whose status is not PASS")

DATASETS <- c("observed", sprintf("shuffle%02d", 1:5))
if (!setequal(manifest$dataset, DATASETS)) {
  stop("manifest datasets are ", paste(manifest$dataset, collapse = ", "))
}

FIMO_PATHS <- setNames(
  file.path(MODULE_DIR, manifest$fimo_output[match(DATASETS, manifest$dataset)]),
  DATASETS
)
EXPECTED_HITS <- setNames(
  as.numeric(manifest$n_hits[match(DATASETS, manifest$dataset)]), DATASETS
)

required_files <- c(FASTA, MOTIF_FILE, unname(FIMO_PATHS))
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

# ------------------------------------------------------------------- helpers
infer_subfamily <- function(locus_id) str_remove(locus_id, "-int_.*$")

safe_wilcox <- function(observed, expected, minimum_loci = 5L) {
  valid <- is.finite(observed) & is.finite(expected)
  observed <- observed[valid]
  expected <- expected[valid]
  if (length(observed) < minimum_loci) return(NA_real_)
  differences <- observed - expected
  if (all(differences == 0)) return(NA_real_)
  suppressWarnings(tryCatch(
    wilcox.test(observed, expected, paired = TRUE,
                alternative = "two.sided", exact = FALSE, correct = TRUE)$p.value,
    error = function(e) NA_real_
  ))
}

format_p_value <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x == 0 ~ "<2.225e-308",
    TRUE ~ format(x, scientific = TRUE, digits = 4)
  )
}

# --------------------------------------- FASTA: locus, subfamily, length, GC
read_fasta_stats <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "ASCII")
  lines <- lines[nzchar(lines)]
  header_positions <- which(startsWith(lines, ">"))
  if (length(header_positions) == 0) stop("No FASTA records found in: ", path)
  record_ends <- c(header_positions[-1] - 1L, length(lines))
  locus_ids <- str_replace(lines[header_positions], "^>(\\S+).*$", "\\1")

  n_records <- length(locus_ids)
  length_bp <- integer(n_records)
  gc_pct <- numeric(n_records)

  for (index in seq_len(n_records)) {
    first_sequence_line <- header_positions[index] + 1L
    last_sequence_line <- record_ends[index]
    if (first_sequence_line > last_sequence_line) {
      sequence <- ""
    } else {
      sequence <- paste0(lines[seq.int(first_sequence_line, last_sequence_line)],
                         collapse = "")
    }
    sequence <- toupper(sequence)
    sequence_length <- nchar(sequence)
    length_bp[index] <- sequence_length
    if (sequence_length == 0) {
      gc_pct[index] <- NA_real_
    } else {
      gc_pct[index] <- 100 * str_count(sequence, "[GC]") / sequence_length
    }
  }

  output <- tibble(locus_id = locus_ids, subfamily = infer_subfamily(locus_ids),
                   length_bp = length_bp, gc_pct = gc_pct)
  if (anyDuplicated(output$locus_id)) stop("Duplicated FASTA identifiers detected")
  output
}

# -------------------------------------------------------------- motif universe
read_motif_library <- function(path) {
  motif_lines <- readLines(path, warn = FALSE)
  motif_lines <- motif_lines[str_detect(motif_lines, "^MOTIF\\s+")]
  if (length(motif_lines) == 0) stop("No MOTIF records found in: ", path)
  split_lines <- str_split(str_squish(motif_lines), "\\s+")
  map_dfr(split_lines, function(fields) {
    motif_id <- fields[2]
    motif_alt_id <- if (length(fields) >= 3) fields[3] else motif_id
    tibble(motif_id = motif_id, motif_alt_id = motif_alt_id,
           rbp = str_extract(motif_alt_id, "^[^|]+"))
  }) %>%
    filter(!is.na(rbp), rbp != "") %>%
    distinct(motif_id, motif_alt_id, rbp)
}

# ------------------------------------------------------------------ FIMO counts
read_fimo_counts <- function(path, dataset, expect_hits) {
  message("Reading ", dataset, ": ", path)
  fimo <- read_tsv(path, comment = "#", show_col_types = FALSE, progress = FALSE,
                   col_types = cols(.default = col_character()))

  required_columns <- c("motif_id", "motif_alt_id", "sequence_name")
  missing_columns <- setdiff(required_columns, names(fimo))
  if (length(missing_columns) > 0) {
    stop(dataset, " is missing columns: ", paste(missing_columns, collapse = ", "))
  }
  # Gate on the accepted Stage 3 / Stage 7 counts before any statistics.
  if (nrow(fimo) != expect_hits) {
    stop(dataset, ": read ", nrow(fimo), " rows, manifest expects ", expect_hits)
  }

  count_column <- paste0("hits_", dataset)

  fimo %>%
    transmute(
      locus_id = str_remove(sequence_name, "__dinuc_shuf\\d+$"),
      rbp = case_when(
        !is.na(motif_alt_id) & motif_alt_id != "" ~ str_extract(motif_alt_id, "^[^|]+"),
        TRUE ~ motif_id
      )
    ) %>%
    filter(!is.na(locus_id), !is.na(rbp), rbp != "") %>%
    count(locus_id, rbp, name = count_column)
}

# ------------------------------------------------------------------ load data
message("Reading FASTA statistics")
fasta_stats <- read_fasta_stats(FASTA)
message("FASTA loci: ", nrow(fasta_stats))
message("Subfamilies: ", n_distinct(fasta_stats$subfamily))

motif_library <- read_motif_library(MOTIF_FILE)
rbp_universe <- sort(unique(motif_library$rbp))

motif_library_summary <- motif_library %>%
  count(rbp, name = "n_motifs") %>%
  arrange(desc(n_motifs), rbp)
write_tsv(motif_library_summary, file.path(OUTDIR, "rbp_motif_library_summary.tsv"))

write_tsv(
  tibble(metric = c("unique_motifs", "unique_rbps"),
         value = c(n_distinct(motif_library$motif_id), length(rbp_universe))),
  file.path(OUTDIR, "motif_library_scope.tsv")
)

count_tables <- imap(FIMO_PATHS,
                     ~ read_fimo_counts(path = .x, dataset = .y,
                                        expect_hits = EXPECTED_HITS[[.y]]))

counts_sparse <- reduce(count_tables, full_join, by = c("locus_id", "rbp"))

hit_columns <- paste0("hits_", DATASETS)
shuffle_columns <- paste0("hits_shuffle", sprintf("%02d", 1:5))

counts_sparse <- counts_sparse %>%
  mutate(across(all_of(hit_columns), ~ replace_na(as.integer(.x), 0L)))

unknown_loci <- setdiff(unique(counts_sparse$locus_id), fasta_stats$locus_id)
if (length(unknown_loci) > 0) {
  stop("FIMO contains locus identifiers absent from the original FASTA. Example: ",
       unknown_loci[1])
}

# Reconciliation gate: the sparse table must sum back to the raw FIMO totals.
for (ds in DATASETS) {
  got <- sum(counts_sparse[[paste0("hits_", ds)]])
  if (got != EXPECTED_HITS[[ds]]) {
    stop(ds, ": sparse counts sum to ", got, ", expected ", EXPECTED_HITS[[ds]])
  }
}

write_tsv(counts_sparse, gzfile(file.path(OUTDIR, "sparse_locus_rbp_counts.tsv.gz")))

# ------------------------------------------------------------ subfamily loop
subfamilies <- sort(unique(fasta_stats$subfamily))
enrichment_results_list <- vector("list", length(subfamilies))
locus_metrics_list <- vector("list", length(subfamilies))

dense_output_path <- file.path(OUTDIR, "locus_rbp_observed_vs_shuffled.tsv.gz")
if (WRITE_DENSE_LOCUS_RBP_TABLE) {
  dense_connection <- gzfile(dense_output_path, open = "wt")
  dense_header_written <- FALSE
}

for (subfamily_index in seq_along(subfamilies)) {
  subfamily_name <- subfamilies[subfamily_index]
  message("[", subfamily_index, "/", length(subfamilies), "] ", subfamily_name)

  loci_subfamily <- fasta_stats %>% filter(subfamily == subfamily_name)
  sparse_subfamily <- counts_sparse %>% filter(locus_id %in% loci_subfamily$locus_id)

  # Complete locus x RBP grid: adds every zero-hit locus/RBP combination.
  dense_subfamily <- expand_grid(locus_id = loci_subfamily$locus_id, rbp = rbp_universe) %>%
    left_join(loci_subfamily, by = "locus_id") %>%
    left_join(sparse_subfamily, by = c("locus_id", "rbp")) %>%
    mutate(across(all_of(hit_columns), ~ replace_na(as.numeric(.x), 0)))

  dense_subfamily$mean_shuffled_hits <-
    rowMeans(as.matrix(dense_subfamily[shuffle_columns]))

  dense_subfamily <- dense_subfamily %>%
    mutate(
      observed_hits = hits_observed,
      observed_hits_per_kb = 1000 * observed_hits / length_bp,
      mean_shuffled_hits_per_kb = 1000 * mean_shuffled_hits / length_bp
    ) %>%
    select(locus_id, subfamily, length_bp, gc_pct, rbp, observed_hits,
           all_of(shuffle_columns), mean_shuffled_hits,
           observed_hits_per_kb, mean_shuffled_hits_per_kb)

  if (WRITE_DENSE_LOCUS_RBP_TABLE) {
    write.table(dense_subfamily, file = dense_connection, sep = "\t",
                row.names = FALSE, col.names = !dense_header_written, quote = FALSE)
    dense_header_written <- TRUE
  }

  total_length_kb <- sum(loci_subfamily$length_bp) / 1000

  enrichment_subfamily <- dense_subfamily %>%
    group_by(rbp) %>%
    summarise(
      subfamily = first(subfamily),
      n_loci = n(),
      observed_total_hits = sum(observed_hits),
      mean_shuffled_total_hits = sum(mean_shuffled_hits),
      observed_hits_per_kb = observed_total_hits / total_length_kb,
      mean_shuffled_hits_per_kb = mean_shuffled_total_hits / total_length_kb,
      density_difference_per_kb = observed_hits_per_kb - mean_shuffled_hits_per_kb,
      log2_observed_expected = log2(
        (observed_total_hits + EFFECT_PSEUDOCOUNT) /
        (mean_shuffled_total_hits + EFFECT_PSEUDOCOUNT)
      ),
      loci_with_observed_hit = sum(observed_hits > 0),
      median_observed_hits_per_locus = median(observed_hits),
      median_shuffled_hits_per_locus = median(mean_shuffled_hits),
      p_value = safe_wilcox(observed_hits, mean_shuffled_hits,
                            minimum_loci = MIN_LOCI_FOR_TEST),
      .groups = "drop"
    ) %>%
    select(subfamily, rbp, everything())

  enrichment_results_list[[subfamily_index]] <- enrichment_subfamily

  locus_metrics <- dense_subfamily %>%
    group_by(locus_id, subfamily, length_bp, gc_pct) %>%
    summarise(
      observed_total_hits = sum(observed_hits),
      observed_unique_rbps = sum(observed_hits > 0),
      mean_shuffled_total_hits = sum(mean_shuffled_hits),
      observed_hits_per_kb = 1000 * observed_total_hits / first(length_bp),
      mean_shuffled_hits_per_kb = 1000 * mean_shuffled_total_hits / first(length_bp),
      .groups = "drop"
    )

  locus_metrics_list[[subfamily_index]] <- locus_metrics

  rm(dense_subfamily, sparse_subfamily, loci_subfamily)
  invisible(gc())
}

if (WRITE_DENSE_LOCUS_RBP_TABLE) close(dense_connection)

# -------------------------------------- multiple testing and classification
enrichment_results <- bind_rows(enrichment_results_list) %>%
  mutate(p_adj_global = p.adjust(p_value, method = "BH")) %>%
  group_by(subfamily) %>%
  mutate(p_adj_within_subfamily = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    # status is decided by the GLOBAL adjustment, not the within-subfamily one.
    status = case_when(
      !is.na(p_adj_global) & p_adj_global < SIGNIFICANCE_ALPHA &
        log2_observed_expected > 0 ~ "enriched",
      !is.na(p_adj_global) & p_adj_global < SIGNIFICANCE_ALPHA &
        log2_observed_expected < 0 ~ "depleted",
      TRUE ~ "not_significant"
    ),
    p_value_display = format_p_value(p_value),
    p_adj_global_display = format_p_value(p_adj_global)
  ) %>%
  arrange(subfamily, p_adj_global, desc(abs(log2_observed_expected)))

write_tsv(enrichment_results,
          file.path(OUTDIR, "rbp_subfamily_composition_controlled_enrichment.tsv"))

# ----------------------------------------------------- subfamily descriptives
locus_metrics_all <- bind_rows(locus_metrics_list) %>%
  # Rename the per-locus density columns BEFORE summarising. Without this,
  # dplyr::summarise() resolves median(observed_hits_per_kb) against the
  # subfamily-level column of the same name created earlier in the same call,
  # so the "median" silently becomes the global density. This rename is present
  # in the historical script and is required for correctness.
  rename(
    locus_observed_hits_per_kb = observed_hits_per_kb,
    locus_mean_shuffled_hits_per_kb = mean_shuffled_hits_per_kb
  )

subfamily_summary <- locus_metrics_all %>%
  group_by(subfamily) %>%
  summarise(
    n_loci = n(),
    total_length_kb = sum(length_bp) / 1000,
    median_length_bp = median(length_bp),
    length_q1_bp = quantile(length_bp, 0.25, names = FALSE),
    length_q3_bp = quantile(length_bp, 0.75, names = FALSE),
    median_gc_pct = median(gc_pct, na.rm = TRUE),
    gc_q1_pct = quantile(gc_pct, 0.25, na.rm = TRUE, names = FALSE),
    gc_q3_pct = quantile(gc_pct, 0.75, na.rm = TRUE, names = FALSE),
    observed_total_hits = sum(observed_total_hits),
    mean_shuffled_total_hits = sum(mean_shuffled_total_hits),
    observed_hits_per_kb = observed_total_hits / total_length_kb,
    mean_shuffled_hits_per_kb = mean_shuffled_total_hits / total_length_kb,
    global_density_ratio = observed_hits_per_kb / mean_shuffled_hits_per_kb,
    median_observed_locus_density = median(locus_observed_hits_per_kb, na.rm = TRUE),
    median_shuffled_locus_density = median(locus_mean_shuffled_hits_per_kb, na.rm = TRUE),
    median_distinct_rbps_per_locus = median(observed_unique_rbps),
    .groups = "drop"
  ) %>%
  arrange(desc(n_loci))

write_tsv(subfamily_summary,
          file.path(OUTDIR, "rbp_subfamily_length_gc_density_summary.tsv"))

# ------------------------------------------------------ focused HERVH/HERVK
hervh_hervk <- enrichment_results %>%
  filter(subfamily %in% c("HERVH", "HERVK")) %>%
  arrange(subfamily, p_adj_global, desc(abs(log2_observed_expected)))

write_tsv(hervh_hervk,
          file.path(OUTDIR, "rbp_HERVH_HERVK_composition_controlled_enrichment.tsv"))

significant_hervh_hervk <- hervh_hervk %>% filter(p_adj_global < SIGNIFICANCE_ALPHA)
write_tsv(significant_hervh_hervk, file.path(OUTDIR, "rbp_HERVH_HERVK_significant.tsv"))

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))

# ------------------------------------------------------------------- console
cat("\n=== MOTIF LIBRARY ===\n")
cat("Unique motifs:", n_distinct(motif_library$motif_id), "\n")
cat("Unique RBPs:", length(rbp_universe), "\n")

cat("\n=== SUBFAMILY SUMMARY: HERVH AND HERVK ===\n")
print(subfamily_summary %>% filter(subfamily %in% c("HERVH", "HERVK")), n = Inf, width = Inf)

cat("\n=== SIGNIFICANT PAIRS ===\n")
print(enrichment_results %>% count(subfamily, status) %>%
        filter(status != "not_significant") %>% arrange(subfamily, status), n = Inf)

cat("\nResults written to:\n", OUTDIR, "\n")
