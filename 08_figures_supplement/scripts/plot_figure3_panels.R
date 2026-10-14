#!/usr/bin/env Rscript
# ==============================================================================
# 08_figures_supplement step 3 (frozen global order 49, effective order 49), LOCAL, offline: Figure 3 panels A and B.
# NO SCIENCE HERE: coverage comes from the accepted 05_lncrna order-35 table (D027) and category membership from the
# order-37 table (D029, the single source of categories). This script classifies nothing: it never compares a coverage
# value with a threshold. Panel B membership is a set of order-37 categories; the dashed lines and their labels come from
# config/parameters.yaml lncrna.coverage_bins and the order-37 category labels.
#
# REFACTOR_SCRIPT of analysis/scripts/scripts_categories/B_transcript_context/06_plot_herv_lncrna_overlap.R
# (sha256 cd34c679..., mtime 2026-04-22 07:52): drawing code verbatim (two ggplot2 histograms, bin width 0.02, right-closed
# bins; panel-B density with adjust 1.2 scaled to counts; theme_classic(base_size = 16); patchwork A | B; 15 x 9 in; ggsave
# SVG (svglite) and PNG (ragg, 300 dpi)). Canonical changes, both required by D027 / D029:
#   x = the exact ratio herv_bp_in_gene_body / gene_body_bp (historically the 4-decimal fraction), so every bar lies in one
#       category and NRIR (0.100038) is drawn above the 10% line;
#   the panel-A labels are the order-37 labels (<=10%, >10-50%, >50-80%, >80%; historically <10%, 10-50%, 50-80%, >80%).
# Panel B: historically fraction > 0.5; canonically the order-37 categories dominant + predominant (exact ratio > 1/2).
# Panel C (HERVarium view of ENSG00000286892, MC05) and the multi-panel ODG assembly (MC07) are manual: not produced here.
# Usage: plot_figure3_panels.R --lncrna-dir 05_lncrna --params config/parameters.yaml --sources-out <dir> --figures-out <dir>
# ==============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(patchwork)
  library(readr)
  library(digest)
})
for (ns in c("svglite", "ragg")) if (!requireNamespace(ns, quietly = TRUE)) stop("required package missing: ", ns, " (ggsave's SVG / PNG device)")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("missing argument ", flag)
  args[[i + 1]]
}
LNC <- get_arg("--lncrna-dir"); PARAMS <- get_arg("--params"); SRC_OUT <- get_arg("--sources-out"); FIG_OUT <- get_arg("--figures-out")

# ------------------------------------------------------------------------------
# Inputs: accepted 05_lncrna contracts, each verified against its results-directory manifest
# ------------------------------------------------------------------------------
verified <- function(rel) {
  p <- file.path(LNC, "results", rel)
  man <- read.table(file.path(dirname(p), "checksums.sha256"), col.names = c("sha", "file"), colClasses = "character")
  want <- man$sha[man$file %in% c(basename(p), paste0("./", basename(p)))]
  if (length(want) != 1 || digest(file = p, algo = "sha256") != want) stop(rel, " does not match the accepted 05_lncrna manifest")
  message("[ok] 05_lncrna/results/", rel, " sha256 ", want)
  list(path = p, sha = want, rel = paste0("results/", rel))
}
IN <- list(
  coverage   = verified("herv_coverage/lncRNA_herv_coverage.tsv"),                          # order 35 (D027)
  categories = verified("coverage_categories/lncRNA_herv_coverage_categories.tsv"),          # order 37 (D029)
  counts     = verified("coverage_categories/lncRNA_herv_coverage_category_counts.tsv")      # order 37 (D029)
)
txt <- function(k) read_tsv(IN[[k]]$path, col_types = cols(.default = col_character()), na = character())

yaml_block <- function(path, block) {
  lines <- readLines(path, warn = FALSE)
  s <- grep(paste0("^\\s+", block, ":\\s*$"), lines)
  if (length(s) != 1) stop("parameters block ", block, " not found")
  ind <- nchar(sub("\\S.*$", "", lines[s]))
  out <- list()
  for (l in lines[(s + 1):length(lines)]) {
    if (!nzchar(trimws(l)) || grepl("^\\s*#", l)) next
    if (nchar(sub("\\S.*$", "", l)) <= ind) break
    l <- sub("\\s+#.*$", "", l)
    out[[trimws(sub(":.*$", "", l))]] <- gsub('^"|"$', "", trimws(sub("^[^:]*:", "", l)))
  }
  out
}
BINS <- yaml_block(PARAMS, "coverage_bins")
sibling <- function(key) {                                   # keys beside coverage_bins in the lncrna block (D029)
  lines <- readLines(PARAMS, warn = FALSE)
  top <- grep("^lncrna:\\s*$", lines); nxt <- grep("^\\S", lines); end <- min(c(nxt[nxt > top], length(lines) + 1)) - 1
  hit <- grep(paste0("^  ", key, ":"), lines[top:end], value = TRUE)
  if (length(hit) != 1) stop("lncrna.", key, " not found once")
  gsub('^"|"$', "", trimws(sub("\\s+#.*$", "", sub("^[^:]*:", "", hit))))
}
BINS$boundary_rule <- sibling("boundary_rule"); BINS$classification_value <- sibling("classification_value")
if (!identical(BINS$boundary_rule, "right_closed") || !identical(BINS$classification_value, "exact_ratio"))
  stop("lncrna.coverage_bins must be right_closed on the exact_ratio (D029)")
CUT <- c(minor = BINS$minor_max, substantial = BINS$substantial_max, dominant = BINS$dominant_max)   # category upper bounds, text

# ------------------------------------------------------------------------------
# The Figure 3 universe: order-37 categories joined to order-35 coverage (consistency guards only)
# ------------------------------------------------------------------------------
cov <- txt("coverage")
cats <- txt("categories")
cnt <- txt("counts")
CATEGORY_ORDER <- c("minor", "substantial", "dominant", "predominant")
tot <- cnt$n[cnt$herv_coverage_category == "total"]
if (anyDuplicated(cats$gene_stable_id) || anyDuplicated(cov$gene_stable_id)) stop("duplicate gene_stable_id in the order-35 / order-37 tables")
if (!identical(cats$gene_stable_id, cov$gene_stable_id)) stop("the order-37 categories and the order-35 coverage do not hold the same genes in the same order")
if (length(tot) != 1 || nrow(cats) != as.integer(tot)) stop("the gene universe (", nrow(cats), ") differs from the order-37 total (", tot, ")")
if (!identical(cats$herv_bp_in_gene_body, cov$herv_bp_in_gene_body) || !identical(cats$gene_body_bp, cov$gene_body_bp))
  stop("the order-37 integers differ from the order-35 coverage")
if (!identical(cats$exact_ratio, paste0(cats$herv_bp_in_gene_body, "/", cats$gene_body_bp))) stop("exact_ratio is not herv_bp_in_gene_body/gene_body_bp")
if (any(!nzchar(cats$herv_bp_in_gene_body) | !nzchar(cats$gene_body_bp) | !nzchar(cats$herv_coverage_category))) stop("missing coverage or category")
if (!all(cats$herv_coverage_category %in% CATEGORY_ORDER)) stop("unknown herv_coverage_category")
lab <- setNames(cnt$herv_coverage_category_label, cnt$herv_coverage_category)[CATEGORY_ORDER]
n_cat <- setNames(as.integer(cnt$n), cnt$herv_coverage_category)[CATEGORY_ORDER]
if (anyNA(lab) || !identical(unname(n_cat), as.integer(table(factor(cats$herv_coverage_category, CATEGORY_ORDER)))))
  stop("category counts of the order-37 table differ from the order-37 count table")
if (!all(cats$herv_coverage_category_label == lab[cats$herv_coverage_category])) stop("category labels inconsistent within order 37")
bounds_txt <- sub("^.*<= ?", "", cnt$interval_on_exact_ratio[match(names(CUT), cnt$herv_coverage_category)])
if (!isTRUE(all(abs(as.numeric(bounds_txt) - as.numeric(CUT)) == 0))) stop("parameters.yaml coverage_bins disagree with the order-37 intervals")

plot_data <- cats %>%
  mutate(
    herv_bp = as.numeric(herv_bp_in_gene_body), body_bp = as.numeric(gene_body_bp),
    exact_ratio_value = herv_bp / body_bp,                     # the exact integer ratio as a double, for drawing only
    display_fraction = gene_body_herv_fraction_display
  )
if (any(abs(plot_data$exact_ratio_value - as.numeric(plot_data$display_fraction)) > 5e-5 + 1e-12)) stop("display fraction is not the 4-decimal rounding of the exact ratio")

# Panel B = the order-37 categories above the substantial_max boundary (historically fraction > 0.5; D029 right-closed)
PANEL_B_CATEGORIES <- c("dominant", "predominant")
plot_tail <- plot_data %>% filter(herv_coverage_category %in% PANEL_B_CATEGORIES)

# ============================================================
# Parameters (historical presentation constants)
# ============================================================

bw_full <- 0.02
bw_tail <- 0.02
X_MAX <- 1                                                   # upper end of the fraction axis
V_LINES <- as.numeric(CUT)                                   # dashed category boundaries (configured thresholds)
B_LOW <- as.numeric(CUT[["substantial"]])                    # panel-B lower edge = the substantial / dominant boundary
B_LINE <- as.numeric(CUT[["dominant"]])                      # panel-B dashed line = the dominant / predominant boundary
LABEL_X <- (c(0, V_LINES) + c(V_LINES, X_MAX)) / 2           # interval midpoints (historical 0.05 / 0.30 / 0.65 / 0.90)
LABELS <- unname(lab)                                        # order-37 labels
B_LABEL_X <- (B_LINE + X_MAX) / 2
B_LABEL <- "HERV-dominated"                                  # historical panel-B annotation

# ============================================================
# Panel A: Full distribution
# ============================================================

p_full <- ggplot(plot_data, aes(x = exact_ratio_value)) +
  geom_histogram(
    binwidth = bw_full,
    boundary = 0,
    fill = "#4C78A8",
    color = "white",
    alpha = 0.7
  ) +
  geom_vline(
    xintercept = V_LINES,
    linetype = "dashed",
    color = "black",
    linewidth = 0.2
  ) +
  annotate("text", x = LABEL_X[1], y = Inf, label = LABELS[1], vjust = 1.8, size = 5.5) +
  annotate("text", x = LABEL_X[2], y = Inf, label = LABELS[2], vjust = 1.8, size = 5.5) +
  annotate("text", x = LABEL_X[3], y = Inf, label = LABELS[3], vjust = 1.8, size = 5.5) +
  annotate("text", x = LABEL_X[4], y = Inf, label = LABELS[4], vjust = 1.8, size = 5.5) +
  scale_x_continuous(
    limits = c(0, X_MAX),
    breaks = seq(0, X_MAX, 0.2),
    labels = label_number(accuracy = 0.1),
    expand = c(0, 0)
  ) +
  labs(
    x = "Fraction of lncRNA gene body overlapping HERV sequence",
    y = "Number of lncRNAs"
  ) +
  theme_classic(base_size = 16)

# ============================================================
# Panel B: Zoom (dominant + predominant)
# ============================================================

p_zoom <- ggplot(plot_tail, aes(x = exact_ratio_value)) +
  geom_histogram(
    binwidth = bw_tail,
    boundary = B_LOW,
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
    xintercept = B_LINE,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  annotate("text", x = B_LABEL_X, y = Inf, label = B_LABEL, vjust = 1.8, size = 5.5) +
  scale_x_continuous(
    limits = c(B_LOW, X_MAX),
    breaks = seq(B_LOW, X_MAX, 0.1),
    labels = label_number(accuracy = 0.1),
    expand = c(0, 0)
  ) +
  labs(
    x = "Fraction of lncRNA gene body overlapping HERV sequence",
    y = "Number of lncRNAs"
  ) +
  theme_classic(base_size = 16)

# ============================================================
# Combine and save (staged, then published)
# ============================================================

plotAB <- (p_full | p_zoom) +
  plot_annotation(tag_levels = "A")

stage_dir <- function(final) {
  parent <- dirname(normalizePath(final, mustWork = FALSE)); dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  s <- tempfile(pattern = ".stage_figure3.", tmpdir = parent); dir.create(s); s
}
ST_SRC <- stage_dir(SRC_OUT); ST_FIG <- stage_dir(FIG_OUT)
options(error = function() { unlink(c(ST_SRC, ST_FIG), recursive = TRUE); quit(save = "no", status = 1) })

ggsave(
  path = ST_FIG,
  filename = "histograms_full_and_zoom.svg",
  plot = plotAB,
  width = 15,
  height = 9
)

ggsave(
  path = ST_FIG,
  filename = "histograms_full_and_zoom.png",
  plot = plotAB,
  width = 15,
  height = 9
)

# ------------------------------------------------------------------------------
# Source contracts (what is drawn, from ggplot2's own layer data; nothing recomputed)
# ------------------------------------------------------------------------------
la <- layer_data(p_full, 1L); lb <- layer_data(p_zoom, 1L); ld <- layer_data(p_zoom, 2L)
bins <- bind_rows(
  tibble(panel = "A", bin = seq_len(nrow(la)), xmin = la$xmin, xmax = la$xmax, count = as.integer(round(la$count))),
  tibble(panel = "B", bin = seq_len(nrow(lb)), xmin = lb$xmin, xmax = lb$xmax, count = as.integer(round(lb$count)))
) %>% mutate(closed = "right", xmin = round(xmin, 10), xmax = round(xmax, 10))
if (sum(bins$count[bins$panel == "A"]) != nrow(plot_data) || sum(bins$count[bins$panel == "B"]) != nrow(plot_tail))
  stop("a plotted gene fell outside the histogram range")
dens <- tibble(point = seq_len(nrow(ld)), x = ld$x, density = ld$density, y_drawn = ld$y, n_observations = ld$n)
gene_cols <- function(d, panel) d %>%
  transmute(gene_stable_id, gene_name, herv_bp_in_gene_body, gene_body_bp, exact_ratio, plotted_value = exact_ratio_value,
            gene_body_herv_fraction_display = display_fraction, gene_body_herv_percent_display, herv_coverage_category,
            herv_coverage_category_label, in_panel_B = herv_coverage_category %in% PANEL_B_CATEGORIES)
src_a <- gene_cols(plot_data)
src_b <- gene_cols(plot_tail) %>% select(-in_panel_B)
summary_tab <- tibble(herv_coverage_category = CATEGORY_ORDER, herv_coverage_category_label = LABELS,
                      upper_bound_parameter = c(CUT, "none"), n_order37 = unname(n_cat),
                      n_panel_A = as.integer(table(factor(plot_data$herv_coverage_category, CATEGORY_ORDER))),
                      n_panel_B = as.integer(table(factor(plot_tail$herv_coverage_category, CATEGORY_ORDER))),
                      label_x_panel_A = LABEL_X)
contract <- tibble(
  panel = c("A,B", "A,B", "A,B", "A,B"),
  upstream_module = c("05_lncrna", "05_lncrna", "05_lncrna", "config"),
  effective_order = c("35 (D027)", "37 (D029)", "37 (D029)", "37 (D029)"),
  upstream_file = c(IN$coverage$rel, IN$categories$rel, IN$counts$rel, "config/parameters.yaml lncrna.coverage_bins"),
  sha256 = c(IN$coverage$sha, IN$categories$sha, IN$counts$sha, digest(file = PARAMS, algo = "sha256")),
  row_unit = c("lncRNA gene (gene body)", "lncRNA gene (gene body)", "coverage category", "parameter"),
  rows_used = c(paste0(nrow(cov), " (consistency: integers, gene set and order)"),
                paste0(nrow(plot_data), " in A; ", nrow(plot_tail), " in B (categories ", paste(PANEL_B_CATEGORIES, collapse = " + "), ")"),
                "4 categories + total (labels, intervals, n)", "minor_max / substantial_max / dominant_max, right_closed, exact_ratio"),
  display_transformation = c("not drawn; herv_bp_in_gene_body and gene_body_bp must equal order 37",
                             paste0("x = herv_bp_in_gene_body / gene_body_bp (exact ratio as a double); right-closed histogram bins of ", bw_full,
                                    " from 0 (A) and from the substantial_max boundary (B); B density: gaussian, bw nrd0 x 1.2, scaled count x ", bw_tail),
                             "panel-A labels = herv_coverage_category_label, placed at interval midpoints",
                             "dashed lines: A at the three upper bounds; B at dominant_max; B axis from substantial_max")
)
write_tsv(src_a, file.path(ST_SRC, "figure3A_source_data.tsv"))
write_tsv(src_b, file.path(ST_SRC, "figure3B_source_data.tsv"))
write_tsv(bins, file.path(ST_SRC, "figure3_histogram_bins.tsv"))
write_tsv(dens, file.path(ST_SRC, "figure3B_density_curve.tsv"))
write_tsv(summary_tab, file.path(ST_SRC, "figure3_category_summary.tsv"))
write_tsv(contract, file.path(ST_SRC, "figure3_source_contracts.tsv"))
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
message("[OK] Figure 3A: ", nrow(plot_data), " lncRNA genes (", paste(CATEGORY_ORDER, summary_tab$n_panel_A, sep = " ", collapse = ", "),
        "); 3B: ", nrow(plot_tail), " genes (", paste(PANEL_B_CATEGORIES, collapse = " + "), ")")
message("[OK] wrote ", SRC_OUT, " and ", FIG_OUT)
