#!/usr/bin/env Rscript
# ==============================================================================
# screen_fmr1_loci.R
# 02_rbp stage 14 (frozen global order 28, effective order 29), LOCAL
# ==============================================================================
# FMR1 motif-guided HERVH locus prioritisation: rank HERVH internal regions by the
# FMR1 RNA-binding-protein motif occurrences they carry at FIMO P <= 1e-5. (FMR1 is
# the RBP motif, not the FMR1 gene; no host-gene field is read.)
#
# Extracted from the historical
#   C_rbp_and_rna_regulation/06_rbp_subfamily_enrichment_and_atlas_FMR1_panel_B.R
# FMR1 branch ONLY. The surrounding raw-count Fisher enrichment, heatmaps and xlsx
# export are superseded (D008/D009) and not migrated; atlas summaries and plots are
# dropped (plotting belongs to 08_figures_supplement).
#
# Rules (historical, unchanged; the verbatim branch reproduces the saved outputs)
#   motif hits : the accepted stage-4 PRIMARY parsed table (results/scanning/
#                primary_fimo_parsed.tsv); never the raw fimo.tsv, never the stage-11
#                strict datasets. The P cut is applied HERE (D016): p-value <= p_cut.
#   RBP        : first '|' token of motif_alt_id == rbp (FMR1: ATtRACT_78, _91, _M016_0.6)
#   subfamily  : '-int_.*$' removed -- from sequence_name for hits, from HERV_id for loci
#   universe   : every locus of the subfamily in the integrated table (zero-hit loci kept)
#   per locus  : hits = n; unique motifs = n_distinct(motif_alt_id); best p / q = min;
#                best score = max; density = 1000 * hits / length; burden and density
#                ranks = min_rank(desc()); density percentile = 100 * percent_rank()
#   candidates : hits >= min_hits and density not NA, ordered by density desc, hits
#                desc, best p asc, HERV_id; candidate_rank = row number; first top_n
#
# Canonical change: the historical loader required February-atlas fields
# (mirna_potential, dsrna_potential, dsrna_structural, dsrna_bidirectional, ...) and the
# exports carried the whole February table through. This script reads only HERV_id,
# locid and length from the accepted 53,176 x 42 table, and writes the screen's own
# columns. The locus chosen for Figure 1C (MC01) is a manual curated decision and is
# NOT computed or flagged here.
#
# Usage:
#   Rscript screen_fmr1_loci.R --fimo <primary_fimo_parsed.tsv> --atlas <transcript_context.with_herv_id.tsv>
#     --results-dir <dir> --figure-dir <dir> --rbp FMR1 --subfamily HERVH --p-cut 1e-5
#     --min-hits 2 --top-n 50 [--expected-fimo-sha256 X] [--expected-atlas-sha256 X]
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

ACCEPTED_FIMO_SHA256  <- "1aeace7924318f2bc14bf5abe567ae1ff67d8e422b23ef85b8976c3f24217cef"  # stage 4
ACCEPTED_ATLAS_SHA256 <- "cbaead1273c916fbcfe4337c7602a9d12c1d5b8bdacfd2a0dd05c0b51f662995"  # effective order 18
PARSED_HEADER <- c("motif_id", "motif_alt_id", "sequence_name", "start", "stop", "strand", "score",
                   "p-value", "q-value", "matched_sequence", "chrom", "internal_start", "internal_end",
                   "internal_strand", "match_strand")

# ------------------------------------------------------------------------------
# CLI and input guards
# ------------------------------------------------------------------------------
parse_args <- function(argv) {
  a <- list(expected_fimo_sha256 = ACCEPTED_FIMO_SHA256, expected_atlas_sha256 = ACCEPTED_ATLAS_SHA256)
  i <- 1
  while (i <= length(argv)) {
    k <- sub("^--", "", argv[[i]]); a[[gsub("-", "_", k)]] <- argv[[i + 1]]; i <- i + 2
  }
  need <- c("fimo", "atlas", "results_dir", "figure_dir", "rbp", "subfamily", "p_cut", "min_hits", "top_n")
  miss <- setdiff(need, names(a))
  if (length(miss)) stop("missing arguments: ", paste(miss, collapse = ", "))
  a$p_cut <- as.numeric(a$p_cut); a$min_hits <- as.integer(a$min_hits); a$top_n <- as.integer(a$top_n)
  a
}

sha256_file <- function(path) substr(system2("sha256sum", shQuote(path), stdout = TRUE)[[1]], 1, 64)

guard_inputs <- function(a) {
  if (grepl("sensitivity_strict|fimo_1e-5|strict", a$fimo))
    stop("stage 14 must not read the stage-11/12 strict datasets (D016); use the primary parsed hits")
  if (grepl("internal_integrated", a$atlas))
    stop("the February atlas is not a production input; use transcript_context.with_herv_id.tsv")
  hdr <- strsplit(readLines(a$fimo, n = 1), "\t")[[1]]
  if (!identical(hdr, PARSED_HEADER))
    stop("--fimo is not the stage-4 primary parsed table (header ", paste(head(hdr, 4), collapse = ","), ", ...)")
  for (x in list(c(a$fimo, a$expected_fimo_sha256), c(a$atlas, a$expected_atlas_sha256))) {
    got <- sha256_file(x[1])
    if (got != x[2]) stop(x[1], ": sha256 ", got, " != expected ", x[2])
    message("[ok] ", basename(x[1]), " sha256 ", got)
  }
}

# ------------------------------------------------------------------------------
# Helpers (historical definitions)
# ------------------------------------------------------------------------------
as_num <- function(x) suppressWarnings(as.numeric(x))
safe_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  suppressWarnings(min(x, na.rm = TRUE))
}
safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  suppressWarnings(max(x, na.rm = TRUE))
}
infer_subfamily <- function(x) str_remove(x, "-int_.*$")

# Parsed primary hits with the historical derived fields (column types guessed by
# read_tsv exactly as historically; p/q/score coerced to numeric).
load_hits <- function(path) {
  read_tsv(path, show_col_types = FALSE) %>%
    mutate(
      subfamily = infer_subfamily(sequence_name),
      rbp       = str_extract(motif_alt_id, "^[^|]+"),
      `p-value` = as_num(`p-value`),
      `q-value` = as_num(`q-value`),
      score     = as_num(score)
    )
}

# Minimal locus contract: the three fields the screen consumes.
load_loci <- function(path) {
  raw <- read_tsv(path, col_types = cols(.default = col_character()), progress = FALSE)
  miss <- setdiff(c("HERV_id", "locid", "length"), names(raw))
  if (length(miss)) stop("locus table lacks: ", paste(miss, collapse = ", "))
  raw %>%
    select(HERV_id, locid, length) %>%
    mutate(subfamily = infer_subfamily(HERV_id), length_num = as_num(length))
}

# ------------------------------------------------------------------------------
# The screen
# ------------------------------------------------------------------------------
screen <- function(hits, loci, a) {
  focus_hits <- hits %>%
    filter(subfamily == a$subfamily, rbp == a$rbp, !is.na(`p-value`), `p-value` <= a$p_cut)

  by_locus <- focus_hits %>%
    group_by(sequence_name) %>%
    summarise(
      focus_rbp_hits = n(),
      focus_rbp_unique_motifs = n_distinct(motif_alt_id),
      focus_rbp_best_p = safe_min(`p-value`),
      focus_rbp_best_q = safe_min(`q-value`),
      focus_rbp_best_score = safe_max(score),
      .groups = "drop"
    )

  all_loci <- loci %>%
    filter(subfamily == a$subfamily) %>%
    left_join(by_locus, by = c("locid" = "sequence_name")) %>%
    mutate(
      focus_rbp_hits = coalesce(focus_rbp_hits, 0L),
      focus_rbp_unique_motifs = coalesce(focus_rbp_unique_motifs, 0L),
      focus_rbp_density_kb = if_else(!is.na(length_num) & length_num > 0,
                                     1000 * focus_rbp_hits / length_num, NA_real_),
      focus_rbp_burden_rank = min_rank(desc(focus_rbp_hits)),
      focus_rbp_density_rank = min_rank(desc(focus_rbp_density_kb)),
      focus_rbp_density_percentile = 100 * percent_rank(focus_rbp_density_kb),
      has_focus_rbp_motif = focus_rbp_hits > 0
    )

  candidates <- all_loci %>%
    filter(focus_rbp_hits >= a$min_hits, !is.na(focus_rbp_density_kb)) %>%
    arrange(desc(focus_rbp_density_kb), desc(focus_rbp_hits), focus_rbp_best_p, HERV_id) %>%
    mutate(candidate_rank = row_number())

  list(focus_hits = focus_hits, by_locus = by_locus, all_loci = all_loci, candidates = candidates)
}

METRIC_COLS <- c("HERV_id", "locid", "subfamily", "length_num", "focus_rbp_hits", "focus_rbp_unique_motifs",
                 "focus_rbp_density_kb", "focus_rbp_best_p", "focus_rbp_best_q", "focus_rbp_best_score",
                 "focus_rbp_burden_rank", "focus_rbp_density_rank", "focus_rbp_density_percentile",
                 "has_focus_rbp_motif")

write_atomic <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".part")
  write_tsv(df, tmp)
  if (!file.rename(tmp, path)) stop("could not write ", path)
  message("[ok] ", path, ": ", nrow(df), " rows x ", ncol(df), " cols")
}

main <- function() {
  a <- parse_args(commandArgs(trailingOnly = TRUE))
  guard_inputs(a)
  hits <- load_hits(a$fimo)
  loci <- load_loci(a$atlas)
  if (anyDuplicated(loci$locid) || anyDuplicated(loci$HERV_id)) stop("locus identifiers not unique")
  s <- screen(hits, loci, a)

  unmatched <- setdiff(s$by_locus$sequence_name, loci$locid)
  if (length(unmatched)) stop(length(unmatched), " qualifying FIMO sequences have no locus: ", paste(head(unmatched), collapse = ", "))
  if (sum(s$all_loci$focus_rbp_hits) != nrow(s$focus_hits)) stop("qualifying hits lost in the locus join")
  message("[qc] ", a$rbp, " rows: primary ", sum(hits$rbp == a$rbp, na.rm = TRUE),
          "; p <= ", a$p_cut, " ", sum(hits$rbp == a$rbp & !is.na(hits$`p-value`) & hits$`p-value` <= a$p_cut, na.rm = TRUE),
          "; ", a$subfamily, " ", nrow(s$focus_hits), " in ", nrow(s$by_locus), " loci of ", nrow(s$all_loci),
          "; candidates (>= ", a$min_hits, " hits) ", nrow(s$candidates))

  prefix <- paste(a$subfamily, a$rbp, sep = "_")
  all_out <- s$all_loci %>%
    select(all_of(METRIC_COLS)) %>%
    arrange(desc(focus_rbp_density_kb), desc(focus_rbp_hits), HERV_id)
  cand_out <- s$candidates %>%
    select(candidate_rank, all_of(METRIC_COLS)) %>%
    slice_head(n = a$top_n)
  hits_out <- s$focus_hits %>%
    left_join(distinct(loci, HERV_id, locid), by = c("sequence_name" = "locid")) %>%
    select(HERV_id, everything())

  write_atomic(all_out,  file.path(a$results_dir, paste0(prefix, "_all_loci.tsv")))
  write_atomic(cand_out, file.path(a$results_dir, paste0(prefix, "_top_candidates.tsv")))
  write_atomic(hits_out, file.path(a$results_dir, paste0(prefix, "_FIMO_hits.tsv")))

  # Figure 1B source data: every candidate, in the historical panel order (density
  # ascending, bottom to top) with its label. The panel-C selection is NOT here: it is
  # the manual curated decision MC01, joined only at plotting time (08_figures_supplement).
  panel_b <- s$candidates %>%
    arrange(focus_rbp_density_kb) %>%
    mutate(panel_order = row_number(), motif_count_label = paste0(focus_rbp_hits, " motifs")) %>%
    select(panel_order, candidate_rank, HERV_id, locid, length_num, focus_rbp_hits,
           focus_rbp_density_kb, motif_count_label)
  write_atomic(panel_b, file.path(a$figure_dir, paste0(prefix, "_panel_B_multiple_motif_loci.tsv")))
  writeLines(capture.output(sessionInfo()), file.path(a$results_dir, "sessionInfo.txt"))
  message("[OK] stage 14 complete")
}

main()
