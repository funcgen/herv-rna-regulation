#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(patchwork)
})

# ============================================================
# Load data
# ============================================================

res_cov <- readRDS("../../../results/lncRNA/res_cov.RDS")

plot_data <- res_cov$full_annotated %>%
  filter(!is.na(gene_body_herv_fraction))

plot_tail <- plot_data %>%
  filter(gene_body_herv_fraction > 0.5)

# ============================================================
# Parameters
# ============================================================

bw_full <- 0.02
bw_tail <- 0.02

# ============================================================
# Panel A: Full distribution
# ============================================================

p_full <- ggplot(plot_data, aes(x = gene_body_herv_fraction)) +
  geom_histogram(
    binwidth = bw_full,
    boundary = 0,
    fill = "#4C78A8",
    color = "white",
    alpha = 0.7
  ) +
  geom_vline(
    xintercept = c(0.1, 0.5, 0.8),
    linetype = "dashed",
    color = "black",
    linewidth = 0.2
  ) +
  annotate("text", x = 0.05, y = Inf, label = "<10%", vjust = 1.8, size = 5.5) +
  annotate("text", x = 0.30, y = Inf, label = "10–50%", vjust = 1.8, size = 5.5) +
  annotate("text", x = 0.65, y = Inf, label = "50–80%", vjust = 1.8, size = 5.5) +
  annotate("text", x = 0.90, y = Inf, label = ">80%", vjust = 1.8, size = 5.5) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = label_number(accuracy = 0.1),
    expand = c(0, 0)
  ) +
  labs(
    x = "Fraction of lncRNA gene body overlapping HERV sequence",
    y = "Number of lncRNAs"
  ) +
  theme_classic(base_size = 16)

# ============================================================
# Panel B: Zoom (>0.5)
# ============================================================

p_zoom <- ggplot(plot_tail, aes(x = gene_body_herv_fraction)) +
  geom_histogram(
    binwidth = bw_tail,
    boundary = 0.5,
    fill = "#4C78A8",
    color = "white",
    alpha = 0.7
  ) +
  geom_density(
    aes(y = after_stat(count * bw_tail)),
    color = "#C21F1F",
    linewidth = 1,
    adjust = 1.2
  ) +
  geom_vline(
    xintercept = 0.8,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  annotate("text", x = 0.9, y = Inf, label = "HERV-dominated", vjust = 1.8, size = 5.5) +
  scale_x_continuous(
    limits = c(0.5, 1),
    breaks = seq(0.5, 1, 0.1),
    labels = label_number(accuracy = 0.1),
    expand = c(0, 0)
  ) +
  labs(
    x = "Fraction of lncRNA gene body overlapping HERV sequence",
    y = "Number of lncRNAs"
  ) +
  theme_classic(base_size = 16)

# ============================================================
# Combine and save
# ============================================================

plotAB <- (p_full | p_zoom) +
  plot_annotation(tag_levels = "A")

ggsave(
  path = "../../../results/figure_herv_lncrna_overlap",
  filename = "histograms_full_and_zoom.svg",
  plot = plotAB,
  width = 15,
  height = 9
)

ggsave(
  path = "../../../results/figure_herv_lncrna_overlap",
  filename = "histograms_full_and_zoom.png",
  plot = plotAB,
  width = 15,
  height = 9
)