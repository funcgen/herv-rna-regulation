#!/usr/bin/env Rscript

# ==============================================================================
# 04_plot_herv_transcript_context_figure.R
# ==============================================================================
# Purpose
#   Assemble the main multi-panel figure summarizing the transcript context of
#   HERV internal regions and conserved HERV domains.
#
# Modes
#   - use_precomputed = TRUE  → read summary tables (recommended for paper)
#   - use_precomputed = FALSE → recompute from raw data
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(forcats)
  library(patchwork)
  library(cowplot)
})

source("utils_transcript_context.R")

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------
CFG <- list(
  use_precomputed = TRUE,
  
  # Raw inputs (used only if recomputing)
  atlas_in  = "../../../data/internal_integrated.with_hervs_id.v2.tsv",
  domain_in = "../../../data/domain_tx_feature_overlaps.wide.tsv",
  
  # Precomputed inputs (from previous scripts)
  panelA_in = "../../../results/figure_herv_transcript_overlap/panelA_global_overlap_context.tsv",
  panelB_in = "../../../results/domain_transcript_context/panelB_terminal_vs_nonterminal_by_gene_type.tsv",
  panelC_in = "../../../results/domain_transcript_context/panelC_terminal_exon_context_in_protein_coding.tsv",
  panelD_in = "../../../results/domain_transcript_context/panelD_domain_composition_terminal_exons.tsv",
  
  # Output
  out_dir    = "../../../results/figure_herv_transcript_overlap/",
  save_plots = TRUE
)

dir.create(CFG$out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Load or compute panel data
# ------------------------------------------------------------------------------

if (isTRUE(CFG$use_precomputed)) {
  
  message("Using precomputed panel tables...")
  
  panelA_clean  <- read_tsv(CFG$panelA_in, show_col_types = FALSE)
  panelB_df     <- read_tsv(CFG$panelB_in, show_col_types = FALSE)
  panelC_summary <- read_tsv(CFG$panelC_in, show_col_types = FALSE)
  panelD_counts <- read_tsv(CFG$panelD_in, show_col_types = FALSE)
  
} else {
  
  message("Recomputing panel tables from raw data...")
  
  # Panel A
  atlas <- load_atlas(CFG$atlas_in)
  atlas_ov_same_strand <- add_overlap_flags(atlas$derived)
  
  atlas_only_lncrna <- atlas_ov_same_strand %>% filter(overlap_lncrna)
  atlas_only_protein_coding <- atlas_ov_same_strand %>% filter(overlap_protein_coding)
  
  panelA_clean <- tibble(
    category = c(
      "lncRNA intronic overlap",
      "lncRNA exonic overlap",
      "Protein-coding intronic overlap",
      "Protein-coding exonic overlap"
    ),
    n = c(
      sum(panelA_df$is_intronic_lncrna & panelA_df$overlap_lncrna, na.rm = TRUE),
      sum(panelA_df$is_exonic_lncrna   & panelA_df$overlap_lncrna, na.rm = TRUE),
      sum(panelA_df$is_intronic_pc     & panelA_df$overlap_protein_coding, na.rm = TRUE),
      sum(panelA_df$is_exonic_pc       & panelA_df$overlap_protein_coding, na.rm = TRUE)
    )
  ) %>%
    mutate(percent = n / nrow(atlas$derived))
  
  # Panels B–D
  df <- read_domain_overlap_table(CFG$domain_in)
  
  df_lncrna <- df %>% filter(gene_type == "lncRNA")
  df_pc     <- df %>% filter(gene_type == "protein_coding")
  
  df_lncrna_last_exon <- df_lncrna %>% filter(any_terminal_exon_overlap)
  df_pc_last_exon     <- df_pc %>% filter(any_terminal_exon_overlap)
  
  panelB_df     <- summarise_panel_b(df)
  panelC_summary <- summarise_panel_c(df_pc_last_exon)
  panelD_counts <- summarise_panel_d(df_lncrna_last_exon, df_pc_last_exon)
}

# ------------------------------------------------------------------------------
# Build plots
# ------------------------------------------------------------------------------

# Panel A
pA <- ggplot(panelA_clean, aes(x = percent, y = fct_rev(category), fill = category)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = percent(percent, accuracy = 0.1)), hjust = -0.15, size = 3.8) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, max(panelA_clean$percent) + 0.05)) +
  scale_fill_manual(values = c(
    "lncRNA intron" = "#9EC5C4",
    "lncRNA exon"   = "#3E7C80",
    "PC intron"     = "#D8C7B8",
    "PC exon"       = "#8C6D5A"
  )) +
  labs(x = "Fraction of all HERV internal regions", y = NULL) +
  theme_pub(base_size = 14) +
  theme(legend.position = "none")

# Panel B
pB <- ggplot(panelB_df, aes(x = gene_type_label, y = percent, fill = exon_position)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = ifelse(percent > 0.05, percent(percent, accuracy = 0.1), "")),
            position = position_stack(vjust = 0.5), color = "white") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_manual(values = c("Non-terminal exon" = "grey75", "Terminal exon" = "#4C5C78")) +
  labs(x = NULL, y = "Fraction of conserved exon domains") +
  theme_pub(base_size = 14)

# Panel C
panelC_summary <- panelC_summary %>%
  mutate(
    context = as.character(context),
    context = recode(
      context,
      "CDS/UTR not defined" = "Not defined"
    )
  ) %>%
  complete(
    context = c(
      "3′ UTR embedded",
      "CDS embedded",
      "5′ UTR embedded",
      "Not defined"
    ),
    fill = list(n = 0, percent = 0)
  ) %>%
  mutate(
    context = factor(
      context,
      levels = rev(c(
        "3′ UTR embedded",
        "CDS embedded",
        "5′ UTR embedded",
        "Not defined"
      ))
    )
  )
pC <- ggplot(panelC_summary, aes(x = percent, y = fct_rev(context))) +
  geom_segment(
    aes(x = 0, xend = percent, y = context, yend = context),
    color = "grey80",
    linewidth = 1
  ) +
  geom_point(aes(fill = context), shape = 21, size = 4.5) +
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
  scale_fill_manual(values = c(
    "3′ UTR embedded" = "#4C5C78",
    "CDS embedded"    = "#A2A9B0",
    "5′ UTR embedded" = "#D7DCE2",
    "Not defined"     = "#C9CDD3"
  )) +
  labs(
    x = "Fraction of terminal-exon conserved domains",
    y = NULL
  ) +
  theme_pub(base_size = 14) +
  theme(
    legend.position = "none",
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(margin = margin(t = 1))
  )




# Panel D
panelD_counts <- panelD_counts %>%
  mutate(
    domain_class = factor(
      domain_class,
      levels = c(
        "Gag",
        "dUTPase",
        "Protease",
        "Reverse transcriptase",
        "RNase H",
        "Integrase",
        "ORFX",
        "Env"
      )
    )
  )
pD <- ggplot(panelD_counts, aes(x = domain_class, y = n, fill = gene_class)) +
  geom_col(position = position_dodge(width = 0.75)) +
  geom_text(aes(label = n), position = position_dodge(width = 0.75), vjust = -0.3) +
  scale_fill_manual(values = fill_gene_class) +
  labs(x = NULL, y = "Number of conserved domain hits") +
  theme_pub(base_size = 14) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

# ------------------------------------------------------------------------------
# Combine figure
# ------------------------------------------------------------------------------
panel_margin <- theme(
  plot.margin = margin(10, 10, 10, 10)
)

pA <- pA + panel_margin
pB <- pB + panel_margin
pC <- pC + panel_margin
pD <- pD + panel_margin

p_multi <- plot_grid(
  pA, pB, pC, pD,
  labels = c("A", "B", "C", "D"),
  ncol = 2,
  label_size = 14,
  label_fontface = "bold"
)


# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------
if (isTRUE(CFG$save_plots)) {
  ggsave(file.path(CFG$out_dir, "FigureX_herv_transcript_overlap.pdf"),
         p_multi, width = 16, height = 12)
  ggsave(file.path(CFG$out_dir, "FigureX_herv_transcript_overlap.png"),
         p_multi, width = 16, height = 12)
}

# ------------------------------------------------------------------------------
# Save tables (useful for debugging/reproducibility)
# ------------------------------------------------------------------------------
# write_tsv(panelA_clean, file.path(CFG$out_dir, "panelA_global_overlap_context.tsv"))
write_tsv(panelB_df, file.path(CFG$out_dir, "panelB_terminal_vs_nonterminal_by_gene_type.tsv"))
write_tsv(panelC_summary, file.path(CFG$out_dir, "panelC_terminal_exon_context_in_protein_coding.tsv"))
write_tsv(panelD_counts, file.path(CFG$out_dir, "panelD_domain_composition_terminal_exons.tsv"))

print(p_multi)
