#!/usr/bin/env Rscript

# ==============================================================================
# 01_build_herv_gene_overlap_tables.R
# ==============================================================================
# Purpose
#   Build the core transcript-context overlap tables for HERV internal regions.
#
# Main questions addressed
#   1. Which HERV internal regions overlap genes on the same strand?
#   2. Among those overlaps, how many involve lncRNAs vs protein-coding genes?
#   3. Which genes are overlapped by HERV internal regions?
#   4. Which lncRNAs show exonic HERV overlaps?
#   5. What fraction of overlaps occur in exons vs introns?
#
# Inputs
#   - ../../../data/internal_integrated.with_hervs_id.v2.tsv
#
# Outputs
#   Written to:
#     ../../../results/lncRNA/
#
#   Main output tables:
#     - lncRNA_overlapping_HERVs.herv_level.tsv
#     - protein_coding_overlapping_HERVs.herv_level.tsv
#     - lncRNA_overlapping_HERVs.gene_level.tsv
#     - protein_coding_overlapping_HERVs.gene_level.tsv
#     - lncRNA_overlapping_HERVs.exonic.herv_level.tsv
#     - lncRNA_overlapping_HERVs.exonic.gene_level.tsv
#     - protein_coding_feature_overlap_summary.tsv
#     - lncRNA_feature_overlap_summary.tsv
#
# Dependencies
#   - utils_transcript_context.R
#
# Notes
#   - This is the first script to run in the transcript-context module.
#   - It depends on the integrated atlas and produces the base overlap tables
#     consumed by downstream annotation / plotting steps.
#   - All heavy lifting is delegated to helper functions defined in
#     utils_transcript_context.R.
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

# ------------------------------------------------------------------------------
# Load shared helper functions
# ------------------------------------------------------------------------------
# This file provides:
#   - atlas loading helpers
#   - overlap parsing helpers
#   - HERV-level / gene-level table builders
#   - small summary helpers
source("utils_transcript_context.R")
# This script assumes it is run from its own directory so that source() works.

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
# atlas_in : integrated atlas used as the source table
# out_dir  : directory where all transcript-context overlap tables are written
CFG <- list(
  atlas_in = "../../../data/transcript_context.with_herv_id.tsv",
  out_dir  = "../../../results/lncRNA/"
)

# Ensure output directory exists
dir.create(CFG$out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1) Load and prepare atlas
# ------------------------------------------------------------------------------
# load_atlas() returns:
#   - raw     : original table
#   - derived : same table with parsed metadata and convenient feature flags
atlas <- load_atlas(CFG$atlas_in)

# Keep only rows with trusted same-strand overlap annotation and add Boolean flags:
#   - overlap_lncrna
#   - overlap_protein_coding
atlas_ov_same_strand <- add_overlap_flags(atlas$derived)

# ------------------------------------------------------------------------------
# 2) Define main subsets used downstream
# ------------------------------------------------------------------------------
# HERVs overlapping lncRNAs on the same strand
atlas_only_lncrna <- atlas_ov_same_strand %>% filter(overlap_lncrna)

# HERVs overlapping protein-coding genes on the same strand
atlas_only_protein_coding <- atlas_ov_same_strand %>% filter(overlap_protein_coding)

# HERVs overlapping lncRNAs specifically at exonic positions
atlas_lncrna_exonic <- atlas_only_lncrna %>%
  filter(has_feature(feature_overlap, "exon"))

# ------------------------------------------------------------------------------
# 3) Print high-level overlap prevalence
# ------------------------------------------------------------------------------
# Percentages are computed relative to all internal regions in the atlas.
message(
  pct(sum(atlas_ov_same_strand$overlap_lncrna, na.rm = TRUE), nrow(atlas$derived)),
  "% of HERV internal regions overlap any lncRNA on the same strand."
)

message(
  pct(sum(atlas_ov_same_strand$overlap_protein_coding, na.rm = TRUE), nrow(atlas$derived)),
  "% of HERV internal regions overlap any protein-coding gene on the same strand."
)

# ------------------------------------------------------------------------------
# 4) Build HERV-level overlap tables
# ------------------------------------------------------------------------------
# One row per HERV internal region, with overlapping genes collapsed into
# semicolon-separated fields.
lnc_herv_overlap_table <- build_herv_level_overlap_table(
  atlas_only_lncrna,
  target = "lncRNA"
)

pc_herv_overlap_table <- build_herv_level_overlap_table(
  atlas_only_protein_coding,
  target = "protein_coding"
)

# Exonic lncRNA-specific HERV-level table
lnc_herv_exonic_overlap_table <- build_herv_level_overlap_table(
  atlas_lncrna_exonic,
  target = "lncRNA",
  include_feature_overlap_per_gene = TRUE
)

# ------------------------------------------------------------------------------
# 5) Build gene-level overlap tables
# ------------------------------------------------------------------------------
# One row per unique gene, summarizing all overlapping HERV internal regions.
lnc_gene_overlap_table <- build_gene_level_overlap_table(
  atlas_only_lncrna,
  target = "lncRNA"
)

pc_gene_overlap_table <- build_gene_level_overlap_table(
  atlas_only_protein_coding,
  target = "protein_coding"
)

# Exonic lncRNA-specific gene-level table
lnc_gene_exonic_overlap_table <- build_gene_level_overlap_table(
  atlas_lncrna_exonic,
  target = "lncRNA",
  exonic_only = TRUE
)

# ------------------------------------------------------------------------------
# 6) Write main overlap tables
# ------------------------------------------------------------------------------
write_tsv(
  lnc_herv_overlap_table,
  file.path(CFG$out_dir, "lncRNA_overlapping_HERVs.herv_level.tsv")
)

write_tsv(
  pc_herv_overlap_table,
  file.path(CFG$out_dir, "protein_coding_overlapping_HERVs.herv_level.tsv")
)

write_tsv(
  lnc_gene_overlap_table,
  file.path(CFG$out_dir, "lncRNA_overlapping_HERVs.gene_level.tsv")
)

write_tsv(
  pc_gene_overlap_table,
  file.path(CFG$out_dir, "protein_coding_overlapping_HERVs.gene_level.tsv")
)

write_tsv(
  lnc_herv_exonic_overlap_table,
  file.path(CFG$out_dir, "lncRNA_overlapping_HERVs.exonic.herv_level.tsv")
)

write_tsv(
  lnc_gene_exonic_overlap_table,
  file.path(CFG$out_dir, "lncRNA_overlapping_HERVs.exonic.gene_level.tsv")
)

# ------------------------------------------------------------------------------
# 7) Summarize exon/intron overlap composition
# ------------------------------------------------------------------------------
# These small summary tables are useful both for QC and for downstream figures.
pc_feature_summary  <- summarise_feature_overlap(atlas_only_protein_coding)
lnc_feature_summary <- summarise_feature_overlap(atlas_only_lncrna)

write_tsv(
  pc_feature_summary,
  file.path(CFG$out_dir, "protein_coding_feature_overlap_summary.tsv")
)

write_tsv(
  lnc_feature_summary,
  file.path(CFG$out_dir, "lncRNA_feature_overlap_summary.tsv")
)

# ------------------------------------------------------------------------------
# 8) Report output locations
# ------------------------------------------------------------------------------
message(
  "Written lncRNA HERV-level table: ",
  file.path(CFG$out_dir, "lncRNA_overlapping_HERVs.herv_level.tsv")
)

message(
  "Written protein-coding HERV-level table: ",
  file.path(CFG$out_dir, "protein_coding_overlapping_HERVs.herv_level.tsv")
)

message(
  "Written lncRNA gene-level table: ",
  file.path(CFG$out_dir, "lncRNA_overlapping_HERVs.gene_level.tsv")
)

message(
  "Written protein-coding gene-level table: ",
  file.path(CFG$out_dir, "protein_coding_overlapping_HERVs.gene_level.tsv")
)

message(
  "Written exonic lncRNA HERV-level table: ",
  file.path(CFG$out_dir, "lncRNA_overlapping_HERVs.exonic.herv_level.tsv")
)

message(
  "Written exonic lncRNA gene-level table: ",
  file.path(CFG$out_dir, "lncRNA_overlapping_HERVs.exonic.gene_level.tsv")
)

message(
  "Written protein-coding feature summary: ",
  file.path(CFG$out_dir, "protein_coding_feature_overlap_summary.tsv")
)

message(
  "Written lncRNA feature summary: ",
  file.path(CFG$out_dir, "lncRNA_feature_overlap_summary.tsv")
)

# ------------------------------------------------------------------------------
# 9) Report simple counts for quick inspection
# ------------------------------------------------------------------------------
message("Number of HERV rows overlapping lncRNA exons: ", nrow(lnc_herv_exonic_overlap_table))
message("Number of unique lncRNA genes with exonic HERV overlaps: ", nrow(lnc_gene_exonic_overlap_table))
message("Number of HERV rows overlapping lncRNAs: ", nrow(lnc_herv_overlap_table))
message("Number of HERV rows overlapping protein-coding genes: ", nrow(pc_herv_overlap_table))
message("Number of unique lncRNA genes overlapped: ", nrow(lnc_gene_overlap_table))
message("Number of unique protein-coding genes overlapped: ", nrow(pc_gene_overlap_table))

# ------------------------------------------------------------------------------
# 10) Optional console summaries for QC
# ------------------------------------------------------------------------------
# Print raw feature_overlap distributions
print(table(atlas_only_protein_coding$feature_overlap, useNA = "ifany"))
print(pc_feature_summary)

print(table(atlas_only_lncrna$feature_overlap, useNA = "ifany"))
print(lnc_feature_summary)
