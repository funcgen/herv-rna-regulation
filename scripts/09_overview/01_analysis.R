#!/usr/bin/env Rscript
# ==============================================================================
# HERVarium — ncRNA layer (Genome-wide atlas)  [v2 input adaptation]
# ==============================================================================
# Author: Tomàs Montserrat-Ayuso
#
# Purpose
#   Generate genome-wide “atlas” summaries for the HERVarium ncRNA layer:
#     (1) Global prevalence of ncRNA features across all HERV internal regions
#     (2) Subfamily-level prevalence and context statistics
#     (3) Selection of “highlight” subfamilies for main figures
#
# Input (v2)
#   ../data/internal_integrated.with_hervs_id.v2.tsv
#
# Outputs (TSV)
#   ../results/atlas_summary/atlas_global_summary.tsv
#   ../results/atlas_summary/atlas_subfamily_summary.tsv
#   ../results/atlas_summary/atlas_subfamilies_to_highlight.tsv
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(tidyr)
})

# -----------------------------
# I/O
# -----------------------------
in_tsv  <- "../../../data/internal_integrated.with_hervs_id.v2.tsv"
out_dir <- "../../../results/atlas_summary"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

DEBUG <- FALSE

# -----------------------------
# Helpers
# -----------------------------
infer_subfamily <- function(x) {
  # "HERVK14-int_00001" -> "HERVK14"
  # "ERVL-E-int_00015"  -> "ERVL-E"
  str_replace(x, "-int_[0-9]+$", "")
}

infer_class <- function(subfamily) {
  case_when(
    str_starts(subfamily, "HERVK")  ~ "HERVK",
    str_starts(subfamily, "HERVH")  ~ "HERVH",
    str_starts(subfamily, "HERVW")  ~ "HERVW",
    str_starts(subfamily, "HERV9")  ~ "HERV9",
    str_starts(subfamily, "HERV17") ~ "HERV17",
    TRUE ~ "other"
  )
}

prev <- function(x) mean(x, na.rm = TRUE)

as_num <- function(x) suppressWarnings(as.numeric(x))
as_int <- function(x) suppressWarnings(as.integer(x))

# -----------------------------
# Load data (robust to "." / blanks)
# -----------------------------
df <- read_tsv(
  in_tsv,
  col_types = cols(.default = col_character()),
  progress = FALSE
)

stopifnot("HERV_id" %in% names(df))

if (DEBUG) {
  message("Checking HERV_id format...")
  message("Fraction matching '-int_<digits>$': ",
          mean(str_detect(df$HERV_id, "-int_[0-9]+$")))
}

# -----------------------------
# Atlas-ready derived columns (v2)
# -----------------------------
# v2 has:
# - ltr5_length / ltr3_length (numeric-ish, may be ".")
# - mirna_potential, dsrna_potential, dsrna_structural, dsrna_bidirectional (0/1-ish, may be ".")
# - domain_count, rbp_burden, rbp_unique, length (numeric-ish)
df_atlas <- df %>%
  mutate(
    subfamily = infer_subfamily(HERV_id),
    class = infer_class(subfamily),
    
    # numeric coercions used downstream
    length_num       = as_num(length),
    domain_count_num = as_num(domain_count),
    rbp_burden_num   = as_num(rbp_burden),
    rbp_unique_num   = as_num(rbp_unique),
    ltr5_len_num     = as_num(ltr5_length),
    ltr3_len_num     = as_num(ltr3_length),
    
    mirna_pot_num    = as_int(mirna_potential),
    dsrna_pot_num    = as_int(dsrna_potential),
    dsrna_struct_num = as_int(dsrna_structural),
    dsrna_bidir_num  = as_int(dsrna_bidirectional),
    
    # Binary feature flags
    has_domain = !is.na(domain_count_num) & domain_count_num > 0,
    
    # v2: infer LTR presence from length > 0 (status is optional/noisy)
    has_ltr5 = !is.na(ltr5_len_num) & ltr5_len_num > 0,
    has_ltr3 = !is.na(ltr3_len_num) & ltr3_len_num > 0,
    has_ltr_any = has_ltr5 | has_ltr3,
    
    has_rbp   = !is.na(rbp_burden_num) & rbp_burden_num > 0,
    has_mirna = (mirna_pot_num == 1),
    
    has_dsrna_any    = (dsrna_pot_num == 1),
    has_dsrna_struct = (dsrna_struct_num == 1),
    has_dsrna_bidir  = (dsrna_bidir_num == 1)
  )

# -----------------------------
# (1) Global prevalence summary
# -----------------------------
global_summary <- df_atlas %>%
  summarise(
    n_total = n(),
    
    n_rbp   = sum(has_rbp, na.rm = TRUE),
    n_mirna = sum(has_mirna, na.rm = TRUE),
    
    n_dsrna_any    = sum(has_dsrna_any, na.rm = TRUE),
    n_dsrna_struct = sum(has_dsrna_struct, na.rm = TRUE),
    n_dsrna_bidir  = sum(has_dsrna_bidir, na.rm = TRUE),
    
    n_domain = sum(has_domain, na.rm = TRUE),
    n_ltr_any = sum(has_ltr_any, na.rm = TRUE),
    
    n_intron     = sum(feature_overlap == "intron", na.rm = TRUE),
    n_exon       = sum(feature_overlap == "exon", na.rm = TRUE),
    n_intergenic = sum(feature_overlap == "intergenic", na.rm = TRUE)
  ) %>%
  mutate(
    across(starts_with("n_"), as.integer),
    
    p_rbp   = n_rbp / n_total,
    p_mirna = n_mirna / n_total,
    
    p_dsrna_any    = n_dsrna_any / n_total,
    p_dsrna_struct = n_dsrna_struct / n_total,
    p_dsrna_bidir  = n_dsrna_bidir / n_total,
    
    p_domain = n_domain / n_total,
    p_ltr_any = n_ltr_any / n_total
  )

write_tsv(global_summary, file.path(out_dir, "atlas_global_summary.tsv"))

baseline_p_dsrna <- global_summary %>% pull(p_dsrna_any) %>% as.numeric()
baseline_p_mirna <- global_summary %>% pull(p_mirna) %>% as.numeric()
stopifnot(length(baseline_p_dsrna) == 1, is.finite(baseline_p_dsrna))
stopifnot(length(baseline_p_mirna) == 1, is.finite(baseline_p_mirna))

if (DEBUG) {
  message("baseline_p_dsrna = ", baseline_p_dsrna)
  message("baseline_p_mirna = ", baseline_p_mirna)
}

# -----------------------------
# (2) Subfamily summary table
# -----------------------------
subfam_summary <- df_atlas %>%
  group_by(subfamily) %>%
  summarise(
    n = n(),
    
    p_rbp   = prev(has_rbp),
    p_mirna = prev(has_mirna),
    
    p_dsrna_any    = prev(has_dsrna_any),
    p_dsrna_struct = prev(has_dsrna_struct),
    p_dsrna_bidir  = prev(has_dsrna_bidir),
    
    p_domain  = prev(has_domain),
    p_ltr_any = prev(has_ltr_any),
    
    p_intron     = prev(feature_overlap == "intron"),
    p_intergenic = prev(feature_overlap == "intergenic"),
    
    med_len        = median(length_num, na.rm = TRUE),
    med_rbp_burden = median(rbp_burden_num, na.rm = TRUE),
    med_rbp_unique = median(rbp_unique_num, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(across(starts_with("p_"), as.numeric)) %>%
  arrange(desc(n))

write_tsv(subfam_summary, file.path(out_dir, "atlas_subfamily_summary.tsv"))

# -----------------------------
# (3) Subfamilies to highlight
# -----------------------------
min_n <- 50

pool <- subfam_summary %>%
  filter(n >= min_n) %>%
  mutate(
    enr_dsrna = round(p_dsrna_any / baseline_p_dsrna, 3),
    enr_mirna = round(p_mirna / baseline_p_mirna, 3)
  )

classic_prefixes <- c("HERVK", "HERVH", "HERVW", "HERV9", "HERV17")
classic_subfams <- pool %>%
  filter(str_detect(subfamily, paste0("^(", paste(classic_prefixes, collapse = "|"), ")"))) %>%
  mutate(reason = "classic")

dsrna_highlights <- pool %>%
  filter(p_dsrna_any >= 0.05 | enr_dsrna >= 3) %>%
  arrange(desc(p_dsrna_any), desc(enr_dsrna), desc(n)) %>%
  slice_head(n = 12) %>%
  mutate(reason = "high_dsRNA")

mirna_highlights <- pool %>%
  filter(p_mirna >= 0.03 | enr_mirna >= 3) %>%
  arrange(desc(p_mirna), desc(enr_mirna), desc(n)) %>%
  slice_head(n = 12) %>%
  mutate(reason = "high_miRNA")

rbp_highlights <- pool %>%
  arrange(desc(med_rbp_burden), desc(med_rbp_unique), desc(n)) %>%
  slice_head(n = 10) %>%
  mutate(reason = "high_RBP")

subfams_to_highlight <- bind_rows(
  classic_subfams,
  dsrna_highlights,
  mirna_highlights,
  rbp_highlights
) %>%
  distinct(subfamily, .keep_all = TRUE) %>%
  arrange(desc(n))

write_tsv(subfams_to_highlight, file.path(out_dir, "atlas_subfamilies_to_highlight.tsv"))

# -----------------------------
# Console report
# -----------------------------
message("== HERVarium ncRNA atlas (v2) ==")
message("Total loci: ", global_summary$n_total[[1]])
message("RBP+: ", global_summary$n_rbp[[1]], " (", round(100 * global_summary$p_rbp[[1]], 2), "%)")
message("miRNA+: ", global_summary$n_mirna[[1]], " (", round(100 * global_summary$p_mirna[[1]], 2), "%)")
message("dsRNA(any)+: ", global_summary$n_dsrna_any[[1]], " (", round(100 * global_summary$p_dsrna_any[[1]], 2), "%)")
message("Outputs written to: ", out_dir)

# -----------------------------
# Reproducibility: session info
# -----------------------------
sess_file <- file.path(out_dir, "sessionInfo.txt")
sess_lines <- capture.output(sessionInfo())

pkg_versions <- c(
  paste0("dplyr: ", as.character(packageVersion("dplyr"))),
  paste0("stringr: ", as.character(packageVersion("stringr"))),
  paste0("readr: ", as.character(packageVersion("readr"))),
  paste0("tidyr: ", as.character(packageVersion("tidyr")))
)

writeLines(
  c(
    "HERVarium ncRNA atlas — reproducibility record",
    paste0("Timestamp: ", format(Sys.time(), tz = "UTC"), " UTC"),
    paste0("Input: ", normalizePath(in_tsv)),
    paste0("Output dir: ", normalizePath(out_dir)),
    "",
    "Package versions (key deps):",
    pkg_versions,
    "",
    "sessionInfo():",
    sess_lines
  ),
  con = sess_file
)

message("Session info written to: ", sess_file)
