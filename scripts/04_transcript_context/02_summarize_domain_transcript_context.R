#!/usr/bin/env Rscript

# ==============================================================================
# 02_summarize_domain_transcript_context.R
# ==============================================================================
# Purpose
#   Summarize the transcript context of conserved HERV domain overlaps.
#
# Main questions addressed
#   1. Are conserved exon-overlapping HERV domains enriched in terminal exons?
#   2. In protein-coding genes, are terminal-exon conserved domains located
#      preferentially in 3′ UTRs, CDS, or 5′ UTRs?
#   3. Which retroviral domain classes are retained in terminal exons of lncRNA
#      and protein-coding transcripts?
#
# Inputs
#   - ../../../data/domain_tx_feature_overlaps.wide.tsv
#
# Outputs
#   Written to:
#     ../../../results/domain_transcript_context/
#
#   Summary tables:
#     - panelB_terminal_vs_nonterminal_by_gene_type.tsv
#     - panelC_terminal_exon_context_in_protein_coding.tsv
#     - panelD_domain_composition_terminal_exons.tsv
#     - panelC_mixed_feature_diagnostic.tsv
#
#   Figures (if save_plots = TRUE):
#     - panelB_terminal_vs_nonterminal_by_gene_type.pdf
#     - panelC_terminal_exon_context_in_protein_coding.pdf
#     - panelD_domain_composition_terminal_exons.pdf
#
# Dependencies
#   - utils_transcript_context.R
#
# Notes
#   - This is the second script to run in the transcript-context module.
#   - It operates on the domain/transcript overlap table and generates the
#     summary tables later used for interpretation and figure assembly.
#   - Plotting code is kept here because these summaries are tightly linked to
#     the corresponding visual outputs.
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(forcats)
})

# ------------------------------------------------------------------------------
# Load shared helper functions
# ------------------------------------------------------------------------------
# Provides:
#   - read_domain_overlap_table()
#   - summarise_panel_b()
#   - summarise_panel_c()
#   - summarise_panel_d()
#   - theme_pub()
#   - domain_label_map()
#   - fill_gene_class
source("utils_transcript_context.R")

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
# domain_in   : domain/transcript overlap table
# out_dir     : output directory for summary tables and figures
# save_plots  : whether to write PDF versions of the plots
CFG <- list(
  domain_in   = "../../../data/domain_tx_feature_overlaps.wide.tsv",
  out_dir     = "../../../results/domain_transcript_context/",
  save_plots  = TRUE
)

# Ensure output directory exists
dir.create(CFG$out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1) Load input table and define working subsets
# ------------------------------------------------------------------------------
# read_domain_overlap_table() standardizes gene_type labels by adding
# gene_type_label (e.g. "Protein-coding" instead of "protein_coding").
df <- read_domain_overlap_table(CFG$domain_in)

# Split by broad transcript class
df_lncrna <- df %>% filter(gene_type == "lncRNA")
df_pc     <- df %>% filter(gene_type == "protein_coding")

# Keep only conserved domains overlapping terminal exons
df_lncrna_last_exon <- df_lncrna %>% filter(any_terminal_exon_overlap)
df_pc_last_exon     <- df_pc %>% filter(any_terminal_exon_overlap)

# ------------------------------------------------------------------------------
# 2) Build summary tables
# ------------------------------------------------------------------------------
# Panel B:
#   Fraction of conserved exon-overlapping domains in terminal vs non-terminal
#   exons, split by transcript class.
panelB_df <- summarise_panel_b(df)

# Panel C:
#   For protein-coding genes only, summarize where terminal-exon conserved
#   domains are embedded (3′ UTR / CDS / 5′ UTR / undefined).
panelC_summary <- summarise_panel_c(df_pc_last_exon)

# Panel D:
#   Domain-class composition of terminal-exon conserved domains, split by
#   transcript class.
panelD_counts <- summarise_panel_d(df_lncrna_last_exon, df_pc_last_exon)

# Diagnostic table:
#   Check whether protein-coding terminal-exon overlaps simultaneously hit
#   multiple feature classes (e.g. CDS + 3′ UTR).
mixed_feature_diagnostic <- df_pc_last_exon %>%
  mutate(
    has_cds  = CDS != 0,
    has_3utr = three_prime_UTR != 0,
    has_5utr = five_prime_UTR != 0
  ) %>%
  count(has_cds, has_3utr, has_5utr)

# ------------------------------------------------------------------------------
# 3) Write summary tables
# ------------------------------------------------------------------------------
write_tsv(
  panelB_df,
  file.path(CFG$out_dir, "panelB_terminal_vs_nonterminal_by_gene_type.tsv")
)

write_tsv(
  panelC_summary,
  file.path(CFG$out_dir, "panelC_terminal_exon_context_in_protein_coding.tsv")
)

write_tsv(
  panelD_counts,
  file.path(CFG$out_dir, "panelD_domain_composition_terminal_exons.tsv")
)

write_tsv(
  mixed_feature_diagnostic,
  file.path(CFG$out_dir, "panelC_mixed_feature_diagnostic.tsv")
)

# ------------------------------------------------------------------------------
# 4) Build plots
# ------------------------------------------------------------------------------
# Panel B plot:
#   Terminal vs non-terminal exon overlap by transcript class
pB_bytype <- ggplot(panelB_df, aes(x = gene_type_label, y = percent, fill = exon_position)) +
  geom_col(width = 0.65, color = NA) +
  geom_text(
    aes(label = ifelse(percent > 0.05, percent(percent, accuracy = 0.1), "")),
    position = position_stack(vjust = 0.5),
    size = 3.8,
    color = "white"
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "Non-terminal exon" = "grey75",
      "Terminal exon" = "#4C5C78"
    )
  ) +
  labs(
    x = NULL,
    y = "Fraction of conserved exon-overlapping domains",
    fill = NULL,
    title = "Terminal exon enrichment is observed in both lncRNA and protein-coding transcripts"
  ) +
  theme_pub() +
  theme(legend.position = "top")

# Panel C plot:
#   Context of terminal-exon conserved domains in protein-coding genes
pC <- ggplot(panelC_summary, aes(x = percent, y = context)) +
  geom_segment(
    aes(x = 0, xend = percent, y = context, yend = context),
    color = "grey80",
    linewidth = 1
  ) +
  geom_point(
    aes(fill = context),
    shape = 21,
    size = 4.8,
    color = "black",
    stroke = 0.3
  ) +
  geom_text(
    aes(label = paste0(n, " (", percent(percent, accuracy = 0.1), ")")),
    nudge_x = 0.02,
    hjust = 0,
    size = 3.8
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, max(panelC_summary$percent) + 0.2),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = c(
      "3′ UTR embedded"     = "#4C5C78",
      "CDS embedded"        = "#A2A9B0",
      "5′ UTR embedded"     = "#D7DCE2",
      "CDS/UTR not defined" = "#C9CDD3"
    )
  ) +
  labs(
    x = "Fraction of terminal-exon conserved domains",
    y = NULL,
    title = "Within terminal exons, conserved HERV domains are most often embedded in 3′ UTRs"
  ) +
  theme_pub() +
  theme(
    legend.position = "none",
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# Panel D plot:
#   Domain composition in terminal exons, split by transcript class
pD_bar <- ggplot(panelD_counts, aes(x = domain_class, y = n, fill = gene_class)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68, color = NA) +
  geom_text(
    aes(label = n),
    position = position_dodge(width = 0.75),
    vjust = -0.3,
    size = 3.4
  ) +
  scale_fill_manual(values = fill_gene_class) +
  labs(
    x = NULL,
    y = "Number of unique conserved domain hits",
    fill = NULL,
    title = "Multiple retroviral protein modules are retained in terminal exons"
  ) +
  theme_pub() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

# ------------------------------------------------------------------------------
# 5) Save plots
# ------------------------------------------------------------------------------
if (isTRUE(CFG$save_plots)) {
  ggsave(
    file.path(CFG$out_dir, "panelB_terminal_vs_nonterminal_by_gene_type.pdf"),
    pB_bytype,
    width = 6,
    height = 4
  )
  
  ggsave(
    file.path(CFG$out_dir, "panelC_terminal_exon_context_in_protein_coding.pdf"),
    pC,
    width = 7,
    height = 4.5
  )
  
  ggsave(
    file.path(CFG$out_dir, "panelD_domain_composition_terminal_exons.pdf"),
    pD_bar,
    width = 7,
    height = 4.5
  )
}

# ------------------------------------------------------------------------------
# 6) Print summaries to console for quick inspection
# ------------------------------------------------------------------------------
print(panelB_df)
print(panelC_summary)
print(panelD_counts)
print(mixed_feature_diagnostic)