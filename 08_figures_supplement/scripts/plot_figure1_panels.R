#!/usr/bin/env Rscript
# ==============================================================================
# 08_figures_supplement step 1 (frozen global order 47, effective order 47), LOCAL, offline: Figure 1 panels A and B.
# NO SCIENCE HERE: every plotted value is read from an accepted 02_rbp contract; this script only reshapes and draws.
#
# REFACTOR_SCRIPT of the plotting code of
#   A  ncRNA_layer/scripts/rbp/10_hervh_hervk_composition_profiles_go_gsea.R section 9 (pheatmap; sha256 cadf0078...)
#   B  C_rbp_and_rna_regulation/06_rbp_subfamily_enrichment_and_atlas_FMR1_panel_B.R "Panel B" block (ggplot2; ac70222f...)
# drawing code verbatim; data and parameters from:
#   A  02_rbp results/enrichment/hervh_hervk_followup/Figure1A_heatmap_values.tsv (66 RBPs significant in >= 1 family at
#      global BH < 0.05; PRIMARY FIMO p <= 1e-4 composition-controlled log2 O/E, capped at +-2) and Figure1A_render_parameters.tsv;
#      the matrix row input order (which fixes hclust's leaf order) = first appearance in the accepted Stage-9 table
#      results/enrichment/subfamily/rbp_HERVH_HERVK_composition_controlled_enrichment.tsv, as historically
#   B  02_rbp figure_source_data/figure1b/HERVH_FMR1_panel_B_multiple_motif_loci.tsv (STRICT FMR1 p <= 1e-5 screen, 4 loci,
#      historical panel order) + MC01_fmr1_locus_selection.tsv (the manual choice, joined here; never derived from rank);
#      subtitle / axis values from config/parameters.yaml rbp.fmr1_locus_screen
# Panel C (HERVarium browser view) is a manual snapshot (effective order 52) and the multi-panel ODG assembly is manual;
# neither is produced here.
# Usage: plot_figure1_panels.R --rbp-dir 02_rbp --params config/parameters.yaml --sources-out <dir> --figures-out <dir>
# ==============================================================================
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
  library(digest)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("missing argument ", flag)
  args[[i + 1]]
}
RBP <- get_arg("--rbp-dir"); PARAMS <- get_arg("--params"); SRC_OUT <- get_arg("--sources-out"); FIG_OUT <- get_arg("--figures-out")

# presentation-only constants of the historical panel B (06 lines 842-918); no scientific meaning
PANEL_B_STYLE <- list(segment_linewidth = 0.7, segment_alpha = 0.5, shape = c(`FALSE` = 16, `TRUE` = 18), size = c(`FALSE` = 3, `TRUE` = 4.5),
                      count_label_hjust = -0.25, count_label_size = 3.5, selected_label = "Selected for panel C", selected_hjust = -0.15,
                      selected_vjust = -1.2, selected_size = 3.5, x_expand = c(0.02, 0.38), width_in = 7, height_in = 7, png_dpi = 300)
PDF_FIXED_DATE <- "D:20000101000000"    # PDF CreationDate / ModDate normalised for byte-level determinism

# ------------------------------------------------------------------------------
# Inputs: accepted 02_rbp contracts, each verified against the 02_rbp module manifest
# ------------------------------------------------------------------------------
man <- read.table(file.path(RBP, "checksums.sha256"), col.names = c("sha", "file"), colClasses = "character")
verified <- function(rel) {
  p <- file.path(RBP, rel)
  want <- man$sha[man$file %in% c(rel, paste0("./", rel))]
  if (length(want) != 1 || digest(file = p, algo = "sha256") != want) stop(rel, " does not match the accepted 02_rbp manifest")
  message("[ok] 02_rbp/", rel, " sha256 ", want)
  p
}
F_VALUES <- verified("results/enrichment/hervh_hervk_followup/Figure1A_heatmap_values.tsv")
F_RENDER <- verified("results/enrichment/hervh_hervk_followup/Figure1A_render_parameters.tsv")
F_STAGE9 <- verified("results/enrichment/subfamily/rbp_HERVH_HERVK_composition_controlled_enrichment.tsv")
F_PANELB <- verified("figure_source_data/figure1b/HERVH_FMR1_panel_B_multiple_motif_loci.tsv")
F_MC01   <- verified("figure_source_data/figure1b/MC01_fmr1_locus_selection.tsv")
if (grepl("sensitivity_strict", F_STAGE9)) stop("Figure 1A uses the PRIMARY Stage-9 results, not the strict sensitivity analysis")

yaml_block <- function(path, block) {
  txt <- readLines(path, warn = FALSE)
  s <- grep(paste0("^\\s+", block, ":\\s*$"), txt)
  if (length(s) != 1) stop("parameters block ", block, " not found")
  ind <- nchar(sub("\\S.*$", "", txt[s]))
  out <- list()
  for (l in txt[(s + 1):length(txt)]) {
    if (!nzchar(trimws(l))) next
    if (nchar(sub("\\S.*$", "", l)) <= ind) break
    l <- sub("\\s+#.*$", "", l)                                        # drop trailing comments
    key <- trimws(sub(":.*$", "", l)); val <- trimws(sub("^[^:]*:", "", l))
    out[[key]] <- gsub('^"|"$', "", val)
  }
  out
}
FMR1P <- yaml_block(PARAMS, "fmr1_locus_screen")
if (is.null(FMR1P$pvalue_cut) || is.na(suppressWarnings(as.numeric(FMR1P$pvalue_cut))) || is.null(FMR1P$candidate_min_hits) || FMR1P$rbp != "FMR1")
  stop("rbp.fmr1_locus_screen parameters (rbp, pvalue_cut, candidate_min_hits) not readable")
rp <- read_tsv(F_RENDER, show_col_types = FALSE, col_types = cols(.default = col_character()))
R_ <- setNames(rp$value, rp$parameter)

# ------------------------------------------------------------------------------
# Panel A source contract
# ------------------------------------------------------------------------------
vals <- read_tsv(F_VALUES, col_types = cols(.default = col_character()))     # passthrough as text; numbers via as.numeric
stage9 <- read_tsv(F_STAGE9, show_col_types = FALSE) %>% filter(subfamily %in% c("HERVH", "HERVK"))
if (nrow(vals) != 2 * length(unique(vals$rbp)) || !setequal(vals$subfamily, c("HERVH", "HERVK")))
  stop("Figure1A_heatmap_values.tsv must hold one HERVH and one HERVK row per RBP")
row_input_order <- intersect(unique(stage9$rbp), unique(vals$rbp))   # historical pivot_wider order
if (length(row_input_order) != length(unique(vals$rbp))) stop("an RBP of Figure 1A is absent from the Stage-9 table")
heatmap_matrix <- vals %>%
  transmute(rbp, subfamily, heatmap_value_capped = as.numeric(heatmap_value_capped)) %>%
  pivot_wider(names_from = subfamily, values_from = heatmap_value_capped) %>%
  mutate(rbp = factor(rbp, levels = row_input_order)) %>% arrange(rbp) %>% mutate(rbp = as.character(rbp)) %>%
  column_to_rownames("rbp") %>%
  as.matrix()
heatmap_matrix <- heatmap_matrix[, c("HERVH", "HERVK"), drop = FALSE]
cap <- as.numeric(R_[["cap_abs"]])
heatmap_colors <- colorRampPalette(c(R_[["colour_low"]], R_[["colour_mid"]], R_[["colour_high"]]))(as.integer(R_[["breaks_n"]]) - 1L)
heatmap_breaks <- seq(-cap, cap, length.out = as.integer(R_[["breaks_n"]]))

plot_heatmap_hervh_hervk <- function(heatmap_matrix, out_file = NULL, width = as.numeric(R_[["width_in"]]), height = as.numeric(R_[["height_in"]])) {
  if (!is.null(out_file)) {
    ext <- tolower(tools::file_ext(out_file))
    if (ext == "pdf") {
      grDevices::pdf(out_file, width = width, height = height)
    } else if (ext == "png") {
      grDevices::png(out_file, width = width, height = height, units = "in", res = as.numeric(R_[["png_dpi"]]))
    } else if (ext == "svg") {
      grDevices::svg(out_file, width = width, height = height)
    } else {
      stop("Unsupported format: use .pdf, .png, or .svg")
    }
  }
  ph <- pheatmap::pheatmap(
    heatmap_matrix,
    scale = "none",
    color = heatmap_colors,
    breaks = heatmap_breaks,
    na_col = R_[["na_colour"]],
    cluster_rows = as.logical(R_[["cluster_rows"]]),
    cluster_cols = as.logical(R_[["cluster_cols"]]),
    border_color = NA,
    fontsize = as.numeric(R_[["fontsize"]]),
    fontsize_row = as.numeric(R_[["fontsize_row"]]),
    fontsize_col = as.numeric(R_[["fontsize_col"]])
  )
  if (!is.null(out_file)) grDevices::dev.off()
  invisible(ph)
}

# ------------------------------------------------------------------------------
# Panel B source contract (MC01 joined; nothing re-ranked or re-selected)
# ------------------------------------------------------------------------------
pb <- read_tsv(F_PANELB, col_types = cols(.default = col_character())) %>% mutate(panel_order = as.integer(panel_order))
mc01 <- read_tsv(F_MC01, show_col_types = FALSE)
if (!identical(sort(pb$HERV_id), sort(mc01$HERV_id)) || sum(mc01$selected_for_panel_c) != 1) stop("MC01 must cover the panel-B loci and select exactly one")
panel_b_df <- pb %>%
  arrange(panel_order) %>%
  left_join(dplyr::select(mc01, HERV_id, selected_for_panel_c, mc01_decision_id = decision_id), by = "HERV_id") %>%
  mutate(HERV_id_plot = factor(HERV_id, levels = HERV_id), focus_rbp_density_plot = as.numeric(focus_rbp_density_kb))
focus_p_label <- format(as.numeric(FMR1P$pvalue_cut), scientific = TRUE)
subtitle_b <- paste0("Loci with at least ", as.integer(FMR1P$candidate_min_hits), " motifs at FIMO P <= ", focus_p_label)
xlab_b <- paste0(FMR1P$rbp, " motif occurrences per kb")

p_panel_b <- ggplot(panel_b_df, aes(y = HERV_id_plot, x = focus_rbp_density_plot)) +
  geom_segment(aes(x = 0, xend = focus_rbp_density_plot, yend = HERV_id_plot), linewidth = PANEL_B_STYLE$segment_linewidth, alpha = PANEL_B_STYLE$segment_alpha) +
  geom_point(aes(shape = selected_for_panel_c, size = selected_for_panel_c)) +
  geom_text(aes(label = motif_count_label), hjust = PANEL_B_STYLE$count_label_hjust, size = PANEL_B_STYLE$count_label_size) +
  geom_text(data = filter(panel_b_df, selected_for_panel_c), aes(label = PANEL_B_STYLE$selected_label),
            hjust = PANEL_B_STYLE$selected_hjust, vjust = PANEL_B_STYLE$selected_vjust, size = PANEL_B_STYLE$selected_size, fontface = "italic") +
  scale_shape_manual(values = PANEL_B_STYLE$shape, guide = "none") +
  scale_size_manual(values = PANEL_B_STYLE$size, guide = "none") +
  scale_x_continuous(expand = expansion(mult = PANEL_B_STYLE$x_expand)) +
  labs(x = xlab_b, y = NULL, subtitle = subtitle_b) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 13),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(size = 14),
    plot.subtitle = element_text(size = 13),
    plot.margin = margin(8, 20, 8, 8)
  )

# ------------------------------------------------------------------------------
# Output (both directories staged, then published)
# ------------------------------------------------------------------------------
stage_dir <- function(final) {
  parent <- dirname(normalizePath(final, mustWork = FALSE)); dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  s <- tempfile(pattern = ".stage_figure1.", tmpdir = parent); dir.create(s); s
}
ST_SRC <- stage_dir(SRC_OUT); ST_FIG <- stage_dir(FIG_OUT)
options(error = function() { unlink(c(ST_SRC, ST_FIG), recursive = TRUE); quit(save = "no", status = 1) })
fa <- "Figure1A_HERVH_HERVK_composition_controlled_heatmap"; fb <- "HERVH_FMR1_panel_B_multiple_motif_loci"
ph <- plot_heatmap_hervh_hervk(heatmap_matrix, file.path(ST_FIG, paste0(fa, ".pdf")))
invisible(plot_heatmap_hervh_hervk(heatmap_matrix, file.path(ST_FIG, paste0(fa, ".png"))))
invisible(plot_heatmap_hervh_hervk(heatmap_matrix, file.path(ST_FIG, paste0(fa, ".svg"))))
pdf_path <- file.path(ST_FIG, paste0(fa, ".pdf"))
raw <- readBin(pdf_path, "raw", file.info(pdf_path)$size)
for (tag in c("/CreationDate (D:", "/ModDate (D:")) {                  # byte-level, same length: xref offsets unchanged
  pos <- grepRaw(tag, raw, fixed = TRUE, all = TRUE)
  if (length(pos) != 1) stop("PDF metadata tag not found once: ", tag)
  at <- pos + nchar(tag) - 2L                                          # start of "D:" + 14 digits
  raw[at:(at + 15L)] <- charToRaw(PDF_FIXED_DATE)
}
writeBin(raw, pdf_path)
ggsave(file.path(ST_FIG, paste0(fb, ".png")), p_panel_b, width = PANEL_B_STYLE$width_in, height = PANEL_B_STYLE$height_in, dpi = PANEL_B_STYLE$png_dpi)
ggsave(file.path(ST_FIG, paste0(fb, ".svg")), p_panel_b, width = PANEL_B_STYLE$width_in, height = PANEL_B_STYLE$height_in)

display <- rownames(heatmap_matrix)[ph$tree_row$order]
src_a <- vals %>%
  mutate(row_input_order = match(rbp, row_input_order), display_row_order = match(rbp, display),
         column_order = match(subfamily, c("HERVH", "HERVK")), plotted_value = heatmap_value_capped,   # text as upstream
         fimo_threshold = "primary p <= 1e-4 (composition-controlled Stage 9)", colour_scale = paste0("[-", cap, ", +", cap, "] ", R_[["colour_low"]], "/", R_[["colour_mid"]], "/", R_[["colour_high"]])) %>%
  arrange(display_row_order, column_order)
src_b <- panel_b_df %>%
  transmute(panel_order, candidate_rank, HERV_id, locid, length_num, focus_rbp_hits, focus_rbp_density_kb, motif_count_label,
            selected_for_panel_c, mc01_decision_id, fimo_threshold = paste0("strict FMR1 post-filter p <= ", focus_p_label),
            subtitle = subtitle_b, x_label = xlab_b)
write_tsv(src_a, file.path(ST_SRC, "figure1A_source_data.tsv"))
write_tsv(src_b, file.path(ST_SRC, "figure1B_source_data.tsv"))
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
message("[OK] Figure 1A: ", nrow(heatmap_matrix), " RBPs x 2; first displayed rows ", paste(head(display, 3), collapse = ", "),
        "; Figure 1B: ", nrow(panel_b_df), " loci; selected (MC01) ", panel_b_df$HERV_id[panel_b_df$selected_for_panel_c])
message("[OK] wrote ", SRC_OUT, " and ", FIG_OUT)
