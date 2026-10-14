#!/usr/bin/env Rscript
# ==============================================================================
# 08_figures_supplement step 5 (frozen global order 51, effective order 51), LOCAL, offline: Figure 5 panels B and C.
# NO SCIENCE HERE: Panel B reads the accepted 07_sparcs order-45 figure5B_source_data.tsv (unique LTRs and unique gene_id per
# branch, D033) and Panel C the accepted order-45 GO_ORA_BP_sparcs_pc_vs_proteinCodingUniverse.tsv (clusterProfiler BP ORA,
# D033). No candidate definition, unique counting, enrichment, multiple-testing correction or GO background exists here.
#
# REFACTOR_SCRIPT of the plotting blocks of analysis/scripts/scripts_categories/E_sparcs_and_dsrna_and_antisense_ltr/
# 01_analysis_sparcs.R (sha256 a017bced..., mtime 2026-04-27 16:08):
#   Panel B  lines 876-1012: two bar charts (patchwork), 8.5 x 4.8 in, PDF + PNG (300 dpi)  -> panel_global_counts_composition.*
#   Panel C  lines 418-451:  GO dot plot of every row of the GO table, 12 x 6 in, SVG     -> go_plot_pc_with_sparks.svg
# drawing code verbatim. The frozen plan asks for PNG and SVG of B and C, so B is also written as SVG and C as PNG.
# Panel A (the SPARCS-like schematic) is hand-drawn on analysis/figures/FigureX3_in_construction.odg and the A-C assembly is
# manual: neither is produced here.
# Usage: plot_figure5_panels.R --sparcs-dir 07_sparcs --sources-out <dir> --figures-out <dir>
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(digest)
})
for (ns in c("svglite", "ragg")) if (!requireNamespace(ns, quietly = TRUE)) stop("required package missing: ", ns, " (ggsave's SVG / PNG device)")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("missing argument ", flag)
  args[[i + 1]]
}
SPA <- get_arg("--sparcs-dir"); SRC_OUT <- get_arg("--sources-out"); FIG_OUT <- get_arg("--figures-out")
PDF_FIXED_DATE <- "D:20000101000000"    # PDF CreationDate / ModDate normalised for byte-level determinism

# ------------------------------------------------------------------------------
# Inputs: accepted order-45 contracts, verified against the results/enrichment manifest
# ------------------------------------------------------------------------------
verified <- function(rel) {
  p <- file.path(SPA, "results", rel)
  man <- read.table(file.path(dirname(p), "checksums.sha256"), col.names = c("sha", "file"), colClasses = "character")
  want <- man$sha[man$file %in% c(basename(p), paste0("./", basename(p)))]
  if (length(want) != 1 || digest(file = p, algo = "sha256") != want) stop(rel, " does not match the accepted 07_sparcs manifest")
  message("[ok] 07_sparcs/results/", rel, " sha256 ", want)
  list(path = p, sha = want, rel = paste0("07_sparcs/results/", rel))
}
IN <- list(
  panelB = verified("enrichment/figure5B_source_data.tsv"),
  panelC = verified("enrichment/GO_ORA_BP_sparcs_pc_vs_proteinCodingUniverse.tsv")
)
b_txt <- read_tsv(IN$panelB$path, col_types = cols(.default = col_character()))
c_txt <- read_tsv(IN$panelC$path, col_types = cols(.default = col_character()))

# guards: the tables must be what the panels show (nothing is recounted)
BRANCHES <- c("coding_3UTR", "lncRNA_last_exon")
if (!identical(b_txt$branch, BRANCHES)) stop("figure5B_source_data.tsv must hold the branches ", paste(BRANCHES, collapse = ", "), " in that order")
for (k in c("n_insertions", "n_affected")) {
  v <- as.integer(b_txt[[k]])
  frac <- as.numeric(b_txt[[sub("^n_(insertion|affected)s?$", "\\1_pct", k)]])
  if (anyNA(v) || any(abs(frac - v / sum(v)) > 1e-12)) stop("figure5B ", k, ": the percentage columns are not the branch shares")
}
if (anyDuplicated(c_txt$ID) || nrow(c_txt) == 0) stop("the GO table must hold one row per GO ID")
ratio_parts <- do.call(rbind, strsplit(c_txt$GeneRatio, "/", fixed = TRUE))
if (ncol(ratio_parts) != 2 || any(ratio_parts[, 1] != c_txt$Count) || any(ratio_parts[, 2] != ratio_parts[1, 2]))
  stop("GeneRatio must be Count / one foreground size")
if (any(lengths(strsplit(c_txt$geneID, "/", fixed = TRUE)) != as.integer(c_txt$Count))) stop("geneID lists do not match Count")

# ============================================================
# Panel B: Global counts / composition of SPARCS-like candidates (historical lines 876-1012)
# ============================================================
BRANCH_LABEL <- setNames(sub(" (", "\n(", b_txt$class, fixed = TRUE), b_txt$branch)    # the historical two-line axis labels

plot_df <- data.table(
  class = factor(
    unname(BRANCH_LABEL[BRANCHES]),
    levels = unname(BRANCH_LABEL[BRANCHES])
  ),
  n_insertions = as.integer(b_txt$n_insertions),
  n_affected = as.integer(b_txt$n_affected)
)

# -----------------------
# Colors
# -----------------------
class_cols <- setNames(c("#4C78A8", "#54A24B"), unname(BRANCH_LABEL[BRANCHES]))

# -----------------------
# Left: unique insertions
# -----------------------
p_insertions <- ggplot(
  plot_df,
  aes(x = class, y = n_insertions, fill = class)
) +
  geom_col(width = 0.68, color = "black", linewidth = 0.25, show.legend = FALSE) +
  geom_text(
    aes(label = comma(n_insertions)),
    vjust = -0.45,
    size = 4.4
  ) +
  scale_fill_manual(values = class_cols) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Unique LTR SPARCS-like insertions",
    x = NULL,
    y = "Count"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 11)
  )

# -----------------------
# Right: affected genes / loci
# -----------------------
p_affected <- ggplot(
  plot_df,
  aes(x = class, y = n_affected, fill = class)
) +
  geom_col(width = 0.68, color = "black", linewidth = 0.25, show.legend = FALSE) +
  geom_text(
    aes(label = comma(n_affected)),
    vjust = -0.45,
    size = 4.4
  ) +
  scale_fill_manual(values = class_cols) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Unique loci with SPARCS-like insertions",
    x = NULL,
    y = "Count"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 11)
  )

# -----------------------
# Combine
# -----------------------
p_counts_panel <- p_insertions + p_affected +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
    )
  )

# ============================================================
# Panel C: GO enrichment plot (historical lines 418-451); every row of the accepted GO table
# ============================================================

# Load data
go <- read_tsv(IN$panelC$path, show_col_types = FALSE)   # historical read (readr types); the contract keeps the text

# Clean and prepare
go <- go %>%
  separate(GeneRatio, into = c("num", "den"), sep = "/", convert = TRUE) %>%
  mutate(GeneRatio_num = num / den)

# Plot
go_plot <- ggplot(go, aes(x = GeneRatio_num,
                   y = reorder(Description, GeneRatio_num))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_continuous(trans = "reverse") +
  labs(
    x = "Gene Ratio",
    y = "",
    color = "Adjusted p-value",
    size = "Gene count"
  ) +
  theme_minimal(base_size = 17)

# ============================================================
# Output (staged, then published)
# ============================================================
stage_dir <- function(final) {
  parent <- dirname(normalizePath(final, mustWork = FALSE)); dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  s <- tempfile(pattern = ".stage_figure5.", tmpdir = parent); dir.create(s); s
}
ST_SRC <- stage_dir(SRC_OUT); ST_FIG <- stage_dir(FIG_OUT)
options(error = function() { unlink(c(ST_SRC, ST_FIG), recursive = TRUE); quit(save = "no", status = 1) })
OUT_PANEL_COUNTS_PDF <- file.path(ST_FIG, "panel_global_counts_composition.pdf")
OUT_PANEL_COUNTS_PNG <- file.path(ST_FIG, "panel_global_counts_composition.png")

# save (historical)
ggsave(OUT_PANEL_COUNTS_PDF, p_counts_panel, width = 8.5, height = 4.8, bg = "white")
ggsave(OUT_PANEL_COUNTS_PNG, p_counts_panel, width = 8.5, height = 4.8, dpi = 300, bg = "white")
ggsave(
  plot = go_plot,
  path = ST_FIG,
  filename = "go_plot_pc_with_sparks.svg",
  width = 12,
  height = 6,
  bg = "white"
)
# formats the frozen plan adds (same plots, same sizes)
ggsave(file.path(ST_FIG, "panel_global_counts_composition.svg"), p_counts_panel, width = 8.5, height = 4.8, bg = "white")
ggsave(file.path(ST_FIG, "go_plot_pc_with_sparks.png"), go_plot, width = 12, height = 6, dpi = 300, bg = "white")

raw <- readBin(OUT_PANEL_COUNTS_PDF, "raw", file.info(OUT_PANEL_COUNTS_PDF)$size)
for (tag in c("/CreationDate (D:", "/ModDate (D:")) {                  # byte-level, same length: xref offsets unchanged
  pos <- grepRaw(tag, raw, fixed = TRUE, all = TRUE)
  if (length(pos) != 1) stop("PDF metadata tag not found once: ", tag)
  at <- pos + nchar(tag) - 2L
  raw[at:(at + 15L)] <- charToRaw(PDF_FIXED_DATE)
}
writeBin(raw, OUT_PANEL_COUNTS_PDF)

# ------------------------------------------------------------------------------
# Source contracts
# ------------------------------------------------------------------------------
ord <- levels(reorder(go$Description, go$GeneRatio_num))                        # bottom-to-top display order
src_b <- b_txt %>% mutate(panel_order = seq_len(n()), axis_label = gsub("\n", " | ", unname(BRANCH_LABEL[branch])),
                          label_insertions = comma(as.integer(n_insertions)), label_affected = comma(as.integer(n_affected)),
                          fill_colour = unname(class_cols[BRANCH_LABEL[branch]]), unit_left = "unique LTR (a_name)",
                          unit_right = "unique gene_id") %>%
  select(panel_order, branch, class, axis_label, n_insertions, n_affected, insertion_pct, affected_pct, label_insertions, label_affected,
         fill_colour, unit_left, unit_right)
src_c <- c_txt %>% mutate(display_rank_top_down = match(Description, rev(ord)), gene_ratio_plotted = as.character(go$GeneRatio_num[match(ID, go$ID)])) %>%
  arrange(display_rank_top_down) %>%
  select(display_rank_top_down, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count, gene_ratio_plotted, geneID)
contract <- tibble(
  panel = c("B", "C"),
  upstream_module = "07_sparcs",
  effective_order = c("45 (D033)", "45 (D033)"),
  upstream_file = c(IN$panelB$rel, IN$panelC$rel),
  sha256 = c(IN$panelB$sha, IN$panelC$sha),
  row_unit = c("branch (coding_3UTR, lncRNA_last_exon); left = unique LTRs (a_name), right = unique gene_id; branches shown separately (no sum or union)",
               "GO BP term (clusterProfiler ORA, BH; foreground = coding_3UTR protein-coding genes; universe = GENCODE v48 protein-coding)"),
  rows_used = c(paste0(nrow(b_txt), " of ", nrow(b_txt)), paste0(nrow(c_txt), " of ", nrow(c_txt), " (every significant term; no subset)")),
  display_transformation = c("bar height = count; label comma(count); two-line axis labels",
                             "x = GeneRatio (Count / foreground); y ordered by GeneRatio; size = Count; colour = p.adjust (reversed continuous scale)")
)
write_tsv(src_b, file.path(ST_SRC, "figure5B_source_data.tsv"))
write_tsv(src_c, file.path(ST_SRC, "figure5C_source_data.tsv"))
write_tsv(contract, file.path(ST_SRC, "figure5_source_contracts.tsv"))
for (d in c(ST_SRC, ST_FIG)) {
  fs <- sort(list.files(d))
  writeLines(paste0(vapply(fs, function(f) digest(file = file.path(d, f), algo = "sha256"), ""), "  ", fs), file.path(d, "checksums.sha256"))
  Sys.chmod(list.files(d, full.names = TRUE), "0644")
}
for (pair in list(c(ST_SRC, SRC_OUT), c(ST_FIG, FIG_OUT))) {
  if (dir.exists(pair[2])) unlink(pair[2], recursive = TRUE)
  if (!file.rename(pair[1], pair[2])) stop("could not publish ", pair[2])
  Sys.chmod(pair[2], "0755")
}
message("[OK] Figure 5B: ", paste(b_txt$branch, b_txt$n_insertions, b_txt$n_affected, sep = " ", collapse = "; "), "; 5C: ", nrow(c_txt),
        " GO terms (", paste(src_c$Description, collapse = " | "), ")")
message("[OK] wrote ", SRC_OUT, " and ", FIG_OUT)
