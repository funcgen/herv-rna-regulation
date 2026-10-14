#!/usr/bin/env Rscript
# ==============================================================================
# 08_figures_supplement step 4 (frozen global order 50, effective order 50), LOCAL, offline: Figure 4 (panels A-C, one composite).
# NO SCIENCE HERE: every plotted value is read from an accepted 06_lncrna_function order-42 figure-source table; membership
# (evidence class, functional categories, ceRNA, coverage category) is never re-derived. No regex, matching, NLP or coverage
# classification exists in this script; it only counts genes per accepted membership, as a consistency guard.
#
# REFACTOR_SCRIPT of the plotting blocks of analysis/scripts/scripts_categories/G_literature_integration/exploration_lit_mining.R
# (sha256 393780db..., mtime 2026-04-21 15:58): Panel A lines 492-560, Panel B 564-658 (plot 628-658), Panel C 663-731,
# assembly 736-771 (cowplot plot_grid 1 x 3, labels A-C size 32, ggsave 20 x 9 in, PNG (ragg, 300 dpi) and SVG (svglite)).
# The historical script output IS the published Figure 4 (Figure4.png = figure_herv_lncrna_annotations.png, byte-identical);
# drawing code verbatim. Canonical change: Panel C uses the order-37 D029 categories (via order 42) and prints the order-37
# labels (<=10%, >10-50%, >50-80%, >80%; historically <10%, 10-50%, 50-80%, >80%).
# Usage: plot_figure4.R --function-dir 06_lncrna_function --lncrna-dir 05_lncrna --sources-out <dir> --figures-out <dir>
# ==============================================================================
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(scales)
  library(cowplot)
  library(digest)
})
for (ns in c("svglite", "ragg")) if (!requireNamespace(ns, quietly = TRUE)) stop("required package missing: ", ns, " (ggsave's SVG / PNG device)")
# cowplot::plot_grid builds the grobs on the current device, as the historical run did on Rscript's default PDF device
# (Rplots.pdf; its metrics shape the layout). The same PDF metrics are used here without writing a file.
grDevices::pdf(NULL)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("missing argument ", flag)
  args[[i + 1]]
}
FUN <- get_arg("--function-dir"); LNC <- get_arg("--lncrna-dir"); SRC_OUT <- get_arg("--sources-out"); FIG_OUT <- get_arg("--figures-out")

# ------------------------------------------------------------------------------
# Inputs: accepted contracts, each verified against its results-directory manifest
# ------------------------------------------------------------------------------
verified <- function(root, module, rel) {
  p <- file.path(root, "results", rel)
  man <- read.table(file.path(dirname(p), "checksums.sha256"), col.names = c("sha", "file"), colClasses = "character")
  want <- man$sha[man$file %in% c(basename(p), paste0("./", basename(p)))]
  if (length(want) != 1 || digest(file = p, algo = "sha256") != want) stop(rel, " does not match the accepted ", module, " manifest")
  message("[ok] ", module, "/results/", rel, " sha256 ", want)
  list(path = p, sha = want, rel = paste0(module, "/results/", rel))
}
IA <- "integrated_annotation/"
IN <- list(
  membership = verified(FUN, "06_lncrna_function", paste0(IA, "annotation_membership.tsv")),
  panelA     = verified(FUN, "06_lncrna_function", paste0(IA, "figure4_panelA_annotation_sources.tsv")),
  panelB     = verified(FUN, "06_lncrna_function", paste0(IA, "figure4_panelB_functional_categories.tsv")),
  panelB_gc  = verified(FUN, "06_lncrna_function", paste0(IA, "figure4_panelB_gene_categories.tsv")),
  panelC     = verified(FUN, "06_lncrna_function", paste0(IA, "figure4_panelC_annotation_by_coverage_category.tsv")),
  cat_counts = verified(LNC, "05_lncrna", "coverage_categories/lncRNA_herv_coverage_category_counts.tsv")
)
txt <- function(k) read_tsv(IN[[k]]$path, col_types = cols(.default = col_character()), na = character())
mem <- txt("membership"); a_txt <- txt("panelA"); b_txt <- txt("panelB"); gc <- txt("panelB_gc"); c_txt <- txt("panelC"); cc <- txt("cat_counts")

# ------------------------------------------------------------------------------
# Consistency guards: the figure tables must be the gene counts of the accepted memberships (nothing is re-derived)
# ------------------------------------------------------------------------------
CATEGORY_ORDER <- c("minor", "substantial", "dominant", "predominant")
if (anyDuplicated(mem$gene_stable_id)) stop("annotation_membership.tsv must hold one row per gene")
if (nrow(mem) != as.integer(cc$n[cc$herv_coverage_category == "total"])) stop("the membership universe differs from the order-37 total")
nA <- table(factor(mem$evidence_class, a_txt$evidence_class))
if (any(is.na(match(mem$evidence_class, a_txt$evidence_class))) || !identical(unname(as.integer(nA)), as.integer(a_txt$n)))
  stop("Panel A counts differ from the order-42 memberships")
annotated <- mem$gene_stable_id[mem$has_any_ann == "1"]
if (!identical(sort(annotated), sort(mem$gene_stable_id[mem$evidence_class != "No annotation"]))) stop("has_any_ann disagrees with the Panel A evidence class")
if (anyDuplicated(gc[c("gene_stable_id", "category")]) || !all(gc$gene_stable_id %in% annotated)) stop("Panel B gene categories must be unique gene x category rows of annotated genes")
nB <- table(factor(gc$category, b_txt$category))
if (any(is.na(match(gc$category, b_txt$category))) || !identical(unname(as.integer(nB)), as.integer(b_txt$n_genes)))
  stop("Panel B counts differ from the order-42 gene x category memberships")
if (!identical(sort(gc$gene_stable_id[gc$category == "ceRNA / miRNA sponge"]), sort(mem$gene_stable_id[mem$cerna_related == "1"])))
  stop("the ceRNA genes of Panel B differ from annotation_membership cerna_related")
if (!identical(c_txt$herv_contribution, CATEGORY_ORDER) || !all(mem$herv_coverage_category %in% CATEGORY_ORDER)) stop("unexpected coverage categories")
nC_total <- as.integer(table(factor(mem$herv_coverage_category, CATEGORY_ORDER)))
nC_ann <- as.integer(table(factor(mem$herv_coverage_category[mem$has_any_ann == "1"], CATEGORY_ORDER)))
if (!identical(nC_total, as.integer(c_txt$n_total)) || !identical(nC_ann, as.integer(c_txt$n_annotated)))
  stop("Panel C counts differ from the order-42 memberships")
if (!identical(nC_total, as.integer(cc$n[match(CATEGORY_ORDER, cc$herv_coverage_category)]))) stop("Panel C denominators differ from the order-37 category counts")
C_LABELS <- setNames(cc$herv_coverage_category_label[match(CATEGORY_ORDER, cc$herv_coverage_category)], CATEGORY_ORDER)

# ============================================================
# Panel A — Functional annotation coverage of HERV-associated lncRNAs
# ============================================================

panelA_df <- a_txt %>%
  mutate(
    n = as.integer(n),
    evidence_class = factor(evidence_class, levels = a_txt$evidence_class),   # order-42 order (historical factor levels)
    fraction = as.numeric(fraction),
    label = paste0(comma(n), "\n(", percent(fraction, accuracy = 0.1), ")")
  )

pA <- ggplot(panelA_df, aes(x = evidence_class, y = n, fill = evidence_class)) +
  geom_col(width = 0.92, color = NA) +
  geom_text(
    aes(label = label),
    vjust = -0.25,
    size = 5.8,
    lineheight = 1
  ) +
  scale_fill_manual(values = c(
    "No annotation"    = "#D9D9D9",
    "GO only"          = "#BFD3E6",
    "LncRNAWiki only"  = "#9ECAE1",
    "Literature only"  = "#6BAED6",
    "Multiple sources" = "#3182BD"
  )) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.10)),
    labels = comma
  ) +
  labs(
    x = NULL,
    y = "HERV-associated lncRNAs"
  ) +
  theme_classic(base_size = 19) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 19),
    axis.text.x = element_text(size = 17, angle = 20, hjust = 1),
    axis.text.y = element_text(size = 17),
    plot.margin = margin(8, 8, 8, 8)
  )

# ============================================================
# Panel B: broad functional categories among annotated HERV-lncRNAs
# ============================================================

LABEL_INSIDE_MIN <- 50   # historical label-placement rule (inside the bar when it holds > 50 genes); presentation only

panelB_df <- b_txt %>%
  mutate(
    n_genes = as.integer(n_genes),
    fraction = as.numeric(fraction),
    label = paste0(percent(fraction, accuracy = 0.1))
  )

pB <- ggplot(panelB_df, aes(x = n_genes, y = fct_reorder(category, n_genes))) +
  geom_col(width = 0.92, fill = "#4C78A8") +
  geom_text(
    aes(
      label = label,
      hjust = ifelse(n_genes > LABEL_INSIDE_MIN, 1.05, -0.1)  # inside if large, outside if small
    ),
    color = ifelse(panelB_df$n_genes > LABEL_INSIDE_MIN, "white", "black"),
    size = 6.8
  ) +
  scale_x_continuous(
    limits = c(0, max(panelB_df$n_genes) * 1.1),
    expand = c(0, 0),
    labels = comma
  ) +
  labs(
    x = "HERV-associated lncRNAs",
    y = NULL
  ) +
  theme_classic(base_size = 19) +
  theme(
    axis.title = element_text(size = 19),
    axis.text = element_text(size = 17)
  )

# ============================================================
# Panel C: annotation by HERV coverage category (order 37 D029 categories, via order 42)
# ============================================================

panelC_df <- c_txt %>%
  mutate(
    n_total = as.integer(n_total),
    n_annotated = as.integer(n_annotated),
    frac_annotated = as.numeric(frac_annotated),
    herv_contribution = factor(
      herv_contribution,
      levels = CATEGORY_ORDER
    ),
    label = paste0(
      percent(frac_annotated, accuracy = 0.1),
      "\n(",
      n_annotated, "/", n_total,
      ")"
    )
  )

pC <- ggplot(panelC_df, aes(x = herv_contribution, y = frac_annotated, group = 1)) +
  geom_line(linewidth = 1.4, color = "black") +
  geom_point(size = 5.5, color = "#4C78A8") +
  geom_point(
    data = panelC_df %>% filter(herv_contribution == "predominant"),
    size = 6.5,
    color = "#C21F1F"
  ) +
  geom_label(
    aes(label = label),
    vjust = -0.5,
    size = 5.5,
    fill = scales::alpha("white", 0.6),
    label.size = 0
  ) +
  scale_x_discrete(
    labels = C_LABELS
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, max(panelC_df$frac_annotated) * 1.10),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "HERV contribution to lncRNA gene body",
    y = "Annotated lncRNAs"
  ) +
  theme_classic(base_size = 19) +
  theme(
    axis.title = element_text(size = 19),
    axis.text = element_text(size = 17),
    plot.margin = margin(8, 8, 8, 8)
  )

# ============================================================
# Assembly (historical cowplot composite) and output (staged, then published)
# ============================================================

plot_margin <- theme(plot.margin = margin(10, 10, 10, 10)) # top, right, bottom, left (in pts)
p_multi <- plot_grid(
  pA + plot_margin,
  pB + plot_margin,
  pC + plot_margin,
  labels = c("A", "B", "C"),
  ncol = 3,
  label_size = 32,
  label_fontface = "plain"
)

stage_dir <- function(final) {
  parent <- dirname(normalizePath(final, mustWork = FALSE)); dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  s <- tempfile(pattern = ".stage_figure4.", tmpdir = parent); dir.create(s); s
}
ST_SRC <- stage_dir(SRC_OUT); ST_FIG <- stage_dir(FIG_OUT)
options(error = function() { unlink(c(ST_SRC, ST_FIG), recursive = TRUE); quit(save = "no", status = 1) })

ggsave(
  filename = file.path(ST_FIG, "figure_herv_lncrna_annotations.png"),
  plot = p_multi,
  width = 20,
  height = 9
)

ggsave(
  filename = file.path(ST_FIG, "figure_herv_lncrna_annotations.svg"),
  plot = p_multi,
  width = 20,
  height = 9
)

# ------------------------------------------------------------------------------
# Source contracts
# ------------------------------------------------------------------------------
b_order <- levels(fct_reorder(panelB_df$category, panelB_df$n_genes))          # bottom-to-top display order
src_a <- a_txt %>% mutate(panel_order = seq_len(n()), label = gsub("\n", " ", panelA_df$label), denominator = as.character(nrow(mem)),
                          fill_colour = c("No annotation" = "#D9D9D9", "GO only" = "#BFD3E6", "LncRNAWiki only" = "#9ECAE1",
                                          "Literature only" = "#6BAED6", "Multiple sources" = "#3182BD")[evidence_class]) %>%
  select(panel_order, evidence_class, n, fraction, denominator, label, fill_colour)
src_b <- b_txt %>% mutate(display_rank_top_down = match(category, rev(b_order)), label = panelB_df$label,
                          denominator = as.character(length(annotated)),
                          label_inside = panelB_df$n_genes > LABEL_INSIDE_MIN) %>%
  select(display_rank_top_down, category, n_genes, fraction, denominator, label, label_inside)
src_c <- c_txt %>% mutate(x_label = unname(C_LABELS[herv_contribution]), label = gsub("\n", " ", panelC_df$label),
                          highlighted = herv_contribution == "predominant") %>%
  select(herv_contribution, x_label, n_total, n_annotated, frac_annotated, label, highlighted)
genes <- mem %>% select(gene_stable_id, gene_name, has_lncRNAwiki, has_lit, has_go, has_any_ann, evidence_class, cerna_related, herv_coverage_category)
contract <- tibble(
  panel = c("A,B,C", "A", "B", "B", "C", "C"),
  upstream_module = c("06_lncrna_function", "06_lncrna_function", "06_lncrna_function", "06_lncrna_function", "06_lncrna_function", "05_lncrna"),
  effective_order = c("42", "42", "42", "42", "42 (D029 categories)", "37 (D029)"),
  upstream_file = vapply(IN[c("membership", "panelA", "panelB", "panelB_gc", "panelC", "cat_counts")], `[[`, "", "rel"),
  sha256 = vapply(IN[c("membership", "panelA", "panelB", "panelB_gc", "panelC", "cat_counts")], `[[`, "", "sha"),
  row_unit = c("lncRNA gene (4,074)", "evidence class (exclusive; 'Multiple sources' = any 2 or 3)", "functional category (unique annotated genes, non-exclusive)",
               "gene x category", "coverage category (gene-body HERV coverage)", "coverage category"),
  rows_used = c(paste0(nrow(mem), " (consistency: A/B/C gene counts, ceRNA set)"), paste0(nrow(a_txt), " of ", nrow(a_txt)),
                paste0(nrow(b_txt), " of ", nrow(b_txt)), paste0(nrow(gc), " (consistency only)"), paste0(nrow(c_txt), " of ", nrow(c_txt)),
                "4 labels + denominators (consistency)"),
  display_transformation = c("not drawn; gene-level inspection table", "bar height = n; label comma(n) and percent(fraction of all genes, 0.1)",
                             "bar length = n_genes, ordered by n_genes; label percent(fraction of annotated genes, 0.1)", "not drawn",
                             "point = frac_annotated (n_annotated / n_total); label percent(0.1) and n_annotated/n_total; predominant in red",
                             "x-axis labels = herv_coverage_category_label")
)
write_tsv(src_a, file.path(ST_SRC, "figure4A_source_data.tsv"))
write_tsv(src_b, file.path(ST_SRC, "figure4B_source_data.tsv"))
write_tsv(src_c, file.path(ST_SRC, "figure4C_source_data.tsv"))
write_tsv(genes, file.path(ST_SRC, "figure4_gene_membership.tsv"))
write_tsv(gc, file.path(ST_SRC, "figure4B_gene_categories.tsv"))
write_tsv(contract, file.path(ST_SRC, "figure4_source_contracts.tsv"))
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
message("[OK] Figure 4A: ", paste(a_txt$evidence_class, a_txt$n, sep = " ", collapse = "; "), "; 4B: ", nrow(b_txt), " categories over ",
        length(annotated), " annotated genes (ceRNA ", b_txt$n_genes[b_txt$category == "ceRNA / miRNA sponge"], "); 4C: ",
        paste(c_txt$n_annotated, c_txt$n_total, sep = "/", collapse = ", "))
message("[OK] wrote ", SRC_OUT, " and ", FIG_OUT)
