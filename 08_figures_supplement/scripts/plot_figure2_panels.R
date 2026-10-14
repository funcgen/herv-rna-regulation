#!/usr/bin/env Rscript
# ==============================================================================
# 08_figures_supplement step 2 (frozen global order 48, effective order 48), LOCAL, offline: Figure 2 panels A to D.
# NO SCIENCE HERE: every plotted value is read from an accepted 03_transcript_context / 04_domains_terminal_exons contract;
# this script only selects the rows the figure shows, reshapes and draws. No test, null or median is computed here; the only
# derived numbers are display quantities (percent scaling, interval / n, the n of the plotted rows, which must equal the upstream n).
#
# REFACTOR_SCRIPT of analysis/scripts/scripts_categories/B_transcript_context/
#   05_plot_revised_herv_transcript_context_figure_v10_no_panelE_log_spaced.R (sha256 5662b8f4..., mtime 2026-08-05 12:50);
# drawing code verbatim (ggplot2 layers, theme, colours, patchwork layout, 12 x 9 in, svglite + ragg 300 dpi). v1-v9 do not migrate.
# Panel map (manuscript snapshot 2026-09-22 and the frozen panel traceability):
#   A  same-strand overlap of HERV internal regions with lncRNA / protein-coding introns and exons
#      03 results/gene_overlap/figure2a_source/panelA_global_overlap_context.tsv (effective order 19; D022 / D024 values);
#      denominator = the 53,176 loci of panelA_global_overlap_context_df.tsv; the four categories are not mutually exclusive
#   B  observed terminal localisation of conserved domains vs the transcript-aware null (exact 95% null range)
#      04 results/positional_null/domain_centered_terminal_exon_null_summary_v2.tsv (effective order 31), primary scope;
#      cross-checked against results/domain_transcript_table/conserved/domain_centered_counts.tsv (order 30).
#      Erratum D039 (2026-10-09): the panel shows the order-31 rows whose terminal rule matches the observed count (>= 80% of the
#      domain in the terminal exon, "stringent"); --panel-b-rule any selects the historical any-overlap rows and is used only by
#      the historical reproduction (validation/figure2/run_historical_reproduction.sh)
#   C  terminal localisation, conserved (HMM coverage >= 0.40) vs lower-coverage predictions, all domains and one per locus
#      04 results/conservation_comparison/ (order 32) and results/locus_representative/ (order 33), contingency + Fisher
#   D  terminal-exon length of multi-exon transcripts: background / non-terminal-domain / terminal-domain
#      04 results/terminal_exon_architecture/transcript_architecture_{groups,summary}.tsv (order 34)
#   E  HERVarium browser view (MC03) and the multi-panel SVG assembly (MC04) are manual: not produced here.
# Usage: plot_figure2_panels.R --context-dir 03_transcript_context --domains-results <04 results dir>
#          --domains-manifest 04_domains_terminal_exons/RESULTS_MANIFEST.tsv --sources-out <dir> --figures-out <dir>
# ==============================================================================
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(digest)
})
for (ns in c("svglite", "ragg")) if (!requireNamespace(ns, quietly = TRUE)) stop("required package missing: ", ns)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("missing argument ", flag)
  args[[i + 1]]
}
CTX <- get_arg("--context-dir"); DOM <- get_arg("--domains-results"); DOM_MAN <- get_arg("--domains-manifest")
SRC_OUT <- get_arg("--sources-out"); FIG_OUT <- get_arg("--figures-out")
PANEL_B_RULE <- if ("--panel-b-rule" %in% args) get_arg("--panel-b-rule") else "stringent"  # D039; "any" = historical display
if (!PANEL_B_RULE %in% c("stringent", "any")) stop("--panel-b-rule must be 'stringent' or 'any'")
B_RULE_TAG <- if (PANEL_B_RULE == "stringent") "stringent_" else ""

# ------------------------------------------------------------------------------
# Inputs: accepted contracts, each verified against its module manifest
# ------------------------------------------------------------------------------
ctx_man <- read.table(file.path(CTX, "results/gene_overlap/checksums.sha256"), col.names = c("sha", "file"), colClasses = "character")
verified_ctx <- function(rel) {
  p <- file.path(CTX, "results/gene_overlap", rel)
  want <- ctx_man$sha[ctx_man$file %in% c(rel, paste0("./", rel))]
  if (length(want) != 1 || digest(file = p, algo = "sha256") != want) stop(rel, " does not match the accepted 03_transcript_context manifest")
  message("[ok] 03_transcript_context/results/gene_overlap/", rel, " sha256 ", want)
  list(path = p, sha = want, rel = paste0("results/gene_overlap/", rel))
}
dom_man <- read_tsv(DOM_MAN, comment = "#", show_col_types = FALSE, col_types = cols(.default = col_character()))
verified_dom <- function(rel) {                         # RESULTS_MANIFEST paths are relative to results/domain_transcript_table/
  key <- if (startsWith(rel, "domain_transcript_table/")) paste0("./", sub("^domain_transcript_table/", "", rel)) else paste0("../", rel)
  want <- dom_man$sha256[dom_man$path == key & dom_man$status == "ACCEPTED"]
  p <- file.path(DOM, rel)
  if (length(want) != 1 || digest(file = p, algo = "sha256") != want) stop(rel, " does not match the accepted 04_domains_terminal_exons manifest")
  message("[ok] 04_domains_terminal_exons/results/", rel, " sha256 ", want)
  list(path = p, sha = want, rel = paste0("results/", rel))
}
IN <- list(
  A_counts   = verified_ctx("figure2a_source/panelA_global_overlap_context.tsv"),
  A_loci     = verified_ctx("figure2a_source/panelA_global_overlap_context_df.tsv"),
  B_null     = verified_dom("positional_null/domain_centered_terminal_exon_null_summary_v2.tsv"),
  B_counts   = verified_dom("domain_transcript_table/conserved/domain_centered_counts.tsv"),
  C_cont     = verified_dom("conservation_comparison/conservation_terminal_exon_contingency.tsv"),
  C_fisher   = verified_dom("conservation_comparison/conservation_terminal_exon_fisher.tsv"),
  C_lcont    = verified_dom("locus_representative/locus_conservation_terminal_exon_contingency.tsv"),
  C_lfisher  = verified_dom("locus_representative/locus_conservation_terminal_exon_fisher.tsv"),
  D_groups   = verified_dom("terminal_exon_architecture/transcript_architecture_groups.tsv"),
  D_summary  = verified_dom("terminal_exon_architecture/transcript_architecture_summary.tsv")
)
txt <- function(k) read_tsv(IN[[k]]$path, col_types = cols(.default = col_character()))   # passthrough as text; numbers via as.numeric

# Which upstream rows the figure shows (the historical v10 selection except panel B, which follows D039; a display choice, recorded in
# the source contracts)
SELECT <- list(
  A_order = c("lncRNA intron", "lncRNA exon", "PC intron", "PC exon"),
  A_flag  = c("lncRNA intron" = "is_intronic_lncrna", "lncRNA exon" = "is_exonic_lncrna", "PC intron" = "is_intronic_pc", "PC exon" = "is_exonic_pc"),
  A_class = c("lncRNA intron" = "overlap_lncrna", "lncRNA exon" = "overlap_lncrna", "PC intron" = "overlap_protein_coding", "PC exon" = "overlap_protein_coding"),
  B_scope = "observed_qualifying_transcripts", B_role = "primary",
  B_analyses = setNames(c("terminal_in_any_transcript", "terminal_in_any_multi_exon_transcript"),   # observed counts (order 30)
                        paste0(c("all_domains_", "multi_exon_domains_"), B_RULE_TAG, "terminal_in_any_transcript")),
  B_n = setNames(c("unique_qualifying_domains", "domains_with_multi_exon_transcript"),
                 paste0(c("all_domains_", "multi_exon_domains_"), B_RULE_TAG, "terminal_in_any_transcript")),
  C_domain = "all_qualifying_domains_terminal_in_any_transcript",
  C_locus  = "one_representative_domain_per_locus_terminal_in_any_transcript",
  D_groups = c("background_gene_without_conserved_domain", "nonterminal_conserved_domain", "terminal_conserved_domain"),
  D_summary = c(level = "transcript", scope = "all_gene_types", metric = "terminal_exon_length_bp")
)

# presentation-only constants of the historical v10 script (lines 65-176); no scientific meaning
CFG <- list(width_in = 12, height_in = 9, dpi = 300, row_gap = 0.08,
            output_stem = "Figure2_herv_terminal_exon_architecture_v10_no_panelE_log_spaced")

format_p <- function(p_value, digits = 2) {
  if (length(p_value) != 1L || is.na(p_value)) {
    return("NA")
  }
  if (p_value < 0.001) {
    return(format(p_value, scientific = TRUE, digits = digits))
  }
  formatC(p_value, format = "f", digits = 3)
}

COL <- list(
  lncRNA_intron = "#9EC5C4",
  lncRNA_exon   = "#3E7C80",
  pc_intron     = "#D8C7B8",
  pc_exon       = "#8C6D5A",

  observed      = "#4C5C78",
  expected      = "#4A4A4A",
  null_interval = "#C9CDD3",

  lower_coverage = "#A2A9B0",
  conserved      = "#4C5C78",

  background     = "#D7DCE2",
  nonterminal    = "#A2A9B0",
  terminal       = "#4C5C78"
)

theme_figure <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans", colour = "#222222"),
      axis.text = element_text(colour = "#222222"),
      axis.title = element_text(colour = "#222222", size = base_size + 0.8),
      plot.margin = margin(10, 18, 10, 10),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_blank(),
      legend.key.width = grid::unit(16, "pt"),
      legend.text = element_text(size = base_size - 0.4),
      panel.grid = element_blank()
    )
}
in_limits <- function(x, lim, what) if (any(is.na(x) | x < lim[1] | x > lim[2])) stop(what, " lies outside the fixed axis limits [", lim[1], ", ", lim[2], "]; ggplot2 would drop it silently")

# ==============================================================================
# PANEL A
# ==============================================================================
a_txt <- txt("A_counts")
if (!identical(sort(a_txt$category), sort(SELECT$A_order)) || nrow(a_txt) != 4) stop("panel A needs exactly the four categories ", paste(SELECT$A_order, collapse = ", "))
loci <- read_tsv(IN$A_loci$path, col_select = c(HERV_id, overlap_lncrna, overlap_protein_coding, is_exonic_lncrna, is_intronic_lncrna, is_exonic_pc, is_intronic_pc),
                 col_types = cols(.default = col_character()))
if (anyDuplicated(loci$HERV_id)) stop("panelA_global_overlap_context_df.tsv must hold one row per HERV internal region")
N_LOCI <- nrow(loci)
for (k in seq_len(nrow(a_txt))) {                       # consistency only (order-19 definition sum(is_X & overlap_Y)); the counts table is drawn
  flag_n <- sum(loci[[SELECT$A_flag[[a_txt$category[k]]]]] == "TRUE" & loci[[SELECT$A_class[[a_txt$category[k]]]]] == "TRUE")
  if (flag_n != as.integer(a_txt$n[k])) stop("panel A ", a_txt$category[k], ": n ", a_txt$n[k], " != ", flag_n, " flagged loci")
  if (abs(as.numeric(a_txt$percent[k]) - flag_n / N_LOCI) > 1e-12) stop("panel A ", a_txt$category[k], ": percent is not n / ", N_LOCI)
}
panelA_data <- a_txt %>% mutate(n = as.numeric(n), percent = as.numeric(percent))

panelA_order <- SELECT$A_order

panelA_labels <- c(
  "lncRNA intron" = "lncRNA intron",
  "lncRNA exon"   = "lncRNA exon",
  "PC intron"     = "Protein-coding intron",
  "PC exon"       = "Protein-coding exon"
)

panelA_colours <- c(
  "lncRNA intron" = COL$lncRNA_intron,
  "lncRNA exon"   = COL$lncRNA_exon,
  "PC intron"     = COL$pc_intron,
  "PC exon"       = COL$pc_exon
)

panelA_plot_data <- panelA_data %>%
  mutate(
    category = factor(category, levels = rev(panelA_order)),
    category_label = factor(
      unname(panelA_labels[as.character(category)]),
      levels = rev(unname(panelA_labels[panelA_order]))
    ),
    value_label = sprintf("%.1f%%  (n = %s)", 100 * percent, comma(n))
  )

panelA_xmax <- max(panelA_plot_data$percent, na.rm = TRUE) * 1.28

pA <- ggplot(
  panelA_plot_data,
  aes(x = percent, y = category_label, fill = category)
) +
  geom_col(width = 0.64) +
  geom_text(
    aes(label = value_label),
    hjust = -0.06,
    size = 3.8,
    colour = "#222222"
  ) +
  scale_x_continuous(
    labels = label_percent(accuracy = 1),
    limits = c(0, panelA_xmax),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_manual(values = panelA_colours, drop = FALSE) +
  labs(
    x = "HERV internal regions (%)",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_figure() +
  theme(
    legend.position = "none",
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 12.2)
  )

# ==============================================================================
# PANEL B
# ==============================================================================
b_txt <- txt("B_null") %>% filter(scope == SELECT$B_scope, scope_role == SELECT$B_role, analysis %in% names(SELECT$B_analyses))
if (nrow(b_txt) != 2 || anyDuplicated(b_txt$analysis)) stop("panel B needs one primary null-summary row per analysis ", paste(names(SELECT$B_analyses), collapse = ", "))
cnt <- txt("B_counts"); cnt_n <- setNames(cnt$n, cnt$metric)
for (k in seq_len(nrow(b_txt))) {                       # the null summary and the order-30 domain counts must agree
  a <- b_txt$analysis[k]
  if (b_txt$n_domains[k] != cnt_n[[SELECT$B_n[[a]]]] || b_txt$observed_terminal_domains[k] != cnt_n[[SELECT$B_analyses[[a]]]])
    stop("panel B ", a, ": the null summary (", b_txt$observed_terminal_domains[k], "/", b_txt$n_domains[k], ") disagrees with domain_centered_counts.tsv")
}
null_summary <- b_txt %>% mutate(across(c(n_domains, observed_terminal_domains, exact_95pct_lower, exact_95pct_upper), as.integer),
                                 across(c(observed_fraction, expected_fraction, exact_p_enrichment), as.numeric))

panelB_analyses <- names(SELECT$B_analyses)

panelB_data <- null_summary %>%
  filter(scope_role == "primary", analysis %in% panelB_analyses) %>%
  mutate(
    y_base = c(2, 1)[match(analysis, panelB_analyses)],   # all conserved domains on top, multi-exon below
    observed_percent = 100 * observed_fraction,
    expected_percent = 100 * expected_fraction,
    lower_percent = 100 * exact_95pct_lower / n_domains,
    upper_percent = 100 * exact_95pct_upper / n_domains,
    observed_label = sprintf("%.1f%% (%d/%d)", observed_percent, observed_terminal_domains, n_domains),
    expected_label = sprintf("%.1f%%", expected_percent),
    p_label = paste0("P = ", vapply(exact_p_enrichment, format_p, character(1), digits = 3))
  ) %>%
  arrange(y_base)
in_limits(unlist(panelB_data[c("observed_percent", "expected_percent", "lower_percent", "upper_percent")]), c(75, 100), "a panel-B value")

pB <- ggplot(panelB_data) +
  geom_segment(
    aes(
      x = lower_percent, xend = upper_percent,
      y = y_base - 0.12, yend = y_base - 0.12
    ),
    linewidth = 6,
    lineend = "round",
    colour = COL$null_interval,
    show.legend = FALSE
  ) +
  geom_point(
    aes(
      x = expected_percent,
      y = y_base - 0.12,
      shape = "Expected"
    ),
    size = 3.6,
    stroke = 1.0,
    fill = "white",
    colour = COL$expected,
    show.legend = TRUE
  ) +
  geom_point(
    aes(
      x = observed_percent,
      y = y_base + 0.12,
      shape = "Observed"
    ),
    size = 3.6,
    colour = COL$observed,
    show.legend = TRUE
  ) +
  geom_text(
    aes(x = observed_percent, y = y_base + 0.12, label = observed_label),
    nudge_y = 0.15,
    vjust = 0,
    size = 3.6,
    colour = COL$observed
  ) +
  geom_text(
    aes(x = expected_percent, y = y_base - 0.12, label = expected_label),
    nudge_y = -0.15,
    vjust = 1,
    size = 3.6,
    colour = COL$expected
  ) +
  geom_text(
    aes(x = 99.4, y = y_base, label = p_label),
    hjust = 1,
    size = 3.8
  ) +
  scale_shape_manual(values = c("Observed" = 16, "Expected" = 21)) +
  guides(
    shape = guide_legend(
      order = 1,
      override.aes = list(
        fill = c("white", COL$observed),
        colour = c(COL$expected, COL$observed)
      )
    )
  ) +
  scale_x_continuous(
    limits = c(75, 100),
    breaks = seq(75, 100, by = 5),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0.55, 2.45),
    breaks = c(1, 2),
    labels = c("Multi-exon sensitivity", "All conserved domains")
  ) +
  labs(
    x = "Domains terminal in at least one transcript (%)",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_figure() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 12.2),
    legend.position = "top",
    legend.justification = "left",
    legend.box.spacing = grid::unit(0, "pt"),
    legend.margin = margin(t = 0, r = 0, b = -12, l = 0),
    legend.box.margin = margin(t = 0, r = 0, b = -8, l = 0)
  )

# ==============================================================================
# PANEL C
# ==============================================================================
c_rows <- function(cont_key, fisher_key, analysis) {
  ct <- txt(cont_key) %>% filter(analysis == !!analysis)
  fi <- txt(fisher_key) %>% filter(analysis == !!analysis)
  if (nrow(fi) != 1L) stop("expected exactly one Fisher row for ", analysis, ", found ", nrow(fi))
  if (nrow(ct) != 2L || !setequal(ct$conservation_group, c("conserved", "lower_coverage"))) stop("expected one conserved and one lower_coverage row for ", analysis)
  g <- function(grp, col) ct[[col]][ct$conservation_group == grp]
  if (fi$conserved_terminal != g("conserved", "terminal") || fi$conserved_nonterminal != g("conserved", "nonterminal") ||
      fi$lower_coverage_terminal != g("lower_coverage", "terminal") || fi$lower_coverage_nonterminal != g("lower_coverage", "nonterminal"))
    stop("the Fisher table and the contingency table disagree for ", analysis)
  list(cont = ct, fisher = fi)
}
cd <- c_rows("C_cont", "C_fisher", SELECT$C_domain); cl <- c_rows("C_lcont", "C_lfisher", SELECT$C_locus)
num_cont <- function(d) d %>% mutate(across(c(terminal, nonterminal, total), as.integer), terminal_percent = as.numeric(terminal_percent))
num_fish <- function(d) d %>% mutate(across(c(odds_ratio, fisher_p_two_sided), as.numeric))
domain_contingency <- num_cont(cd$cont); locus_contingency <- num_cont(cl$cont)
domain_fisher <- num_fish(cd$fisher); locus_fisher <- num_fish(cl$fisher)

thr <- unique(c(cd$fisher$coverage_threshold, cl$fisher$coverage_threshold))
if (length(thr) != 1L || is.na(suppressWarnings(as.numeric(thr)))) stop("panel C: one coverage threshold expected across the Fisher tables")
LAB_LOW <- sprintf("Coverage <%.2f", as.numeric(thr)); LAB_CONS <- sprintf("Coverage \u2265%.2f", as.numeric(thr))   # historical "Coverage <0.40" / "Coverage >=0.40"

domain_primary_analysis <- SELECT$C_domain
locus_primary_analysis  <- SELECT$C_locus

domain_points <- domain_contingency %>%
  filter(analysis == domain_primary_analysis) %>%
  mutate(analysis_level = "All domain predictions", y_base = 1.85)

locus_points <- locus_contingency %>%
  filter(analysis == locus_primary_analysis) %>%
  mutate(analysis_level = "One representative per HERV locus", y_base = 1.15)

panelC_points <- bind_rows(domain_points, locus_points) %>%
  mutate(
    conservation_group = factor(conservation_group, levels = c("lower_coverage", "conserved")),
    point_label = sprintf("%.1f%% (%d/%d)", terminal_percent, terminal, total),
    label_colour = ifelse(conservation_group == "conserved", COL$conserved, COL$lower_coverage)
  )
in_limits(panelC_points$terminal_percent, c(55, 108), "a panel-C value")

domain_stats <- domain_fisher %>%
  transmute(
    y_base = 1.85,
    odds_ratio,
    fisher_p_two_sided
  )

locus_stats <- locus_fisher %>%
  transmute(
    y_base = 1.15,
    odds_ratio,
    fisher_p_two_sided
  )

panelC_stats <- bind_rows(domain_stats, locus_stats) %>%
  mutate(
    stats_label = sprintf(
      "OR %.2f\nP = %s",
      odds_ratio,
      vapply(fisher_p_two_sided, format_p, character(1), digits = 2)
    )
  )

panelC_segments <- panelC_points %>%
  select(y_base, conservation_group, terminal_percent) %>%
  tidyr::pivot_wider(names_from = conservation_group, values_from = terminal_percent)

pC <- ggplot() +
  geom_segment(
    data = panelC_segments,
    aes(x = lower_coverage, xend = conserved, y = y_base, yend = y_base),
    linewidth = 1.4,
    colour = "#B9BEC5"
  ) +
  geom_point(
    data = panelC_points %>% filter(conservation_group == "lower_coverage"),
    aes(x = terminal_percent, y = y_base, shape = LAB_LOW),
    size = 3.7,
    stroke = 1.0,
    fill = "white",
    colour = COL$lower_coverage,
    show.legend = TRUE
  ) +
  geom_point(
    data = panelC_points %>% filter(conservation_group == "conserved"),
    aes(x = terminal_percent, y = y_base, shape = LAB_CONS),
    size = 3.7,
    colour = COL$conserved,
    show.legend = TRUE
  ) +
  geom_text(
    data = panelC_points,
    aes(x = terminal_percent, y = y_base, label = point_label),
    nudge_y = 0.18,
    vjust = 0,
    size = 3.5,
    colour = panelC_points$label_colour,
    show.legend = FALSE
  ) +
  geom_text(
    data = panelC_stats,
    aes(x = 106.5, y = y_base, label = stats_label),
    hjust = 1,
    size = 3.6
  ) +
  scale_shape_manual(
    values = setNames(c(21, 16), c(LAB_LOW, LAB_CONS))
  ) +
  guides(
    shape = guide_legend(
      order = 1,
      override.aes = list(
        fill = c("white", COL$conserved),
        colour = c(COL$lower_coverage, COL$conserved)
      )
    )
  ) +
  scale_x_continuous(
    limits = c(55, 108),
    breaks = seq(60, 100, by = 10),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0.82, 2.18),
    breaks = c(1.15, 1.85),
    labels = c("One representative\nper HERV locus", "All domain\npredictions")
  ) +
  labs(
    x = "Terminal domains or representative loci (%)",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_figure() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 12.2),
    legend.position = "top",
    legend.justification = "left",
    legend.box.spacing = grid::unit(0, "pt"),
    legend.margin = margin(t = 0, r = 0, b = -12, l = 0),
    legend.box.margin = margin(t = 0, r = 0, b = -8, l = 0)
  )

# ==============================================================================
# PANEL D
# ==============================================================================
transcript_groups <- read_tsv(
  IN$D_groups$path,
  col_select = c(tx_id, is_single_exon_transcript, terminal_exon_length_bp, transcript_domain_group),
  col_types = cols(.default = col_character())
) %>% mutate(is_single_exon_transcript = as.logical(is_single_exon_transcript), terminal_exon_length_bp = as.numeric(terminal_exon_length_bp))
if (anyDuplicated(transcript_groups$tx_id)) stop("transcript_architecture_groups.tsv must hold one row per transcript")
if (anyNA(transcript_groups$is_single_exon_transcript)) stop("is_single_exon_transcript must be TRUE / FALSE")
transcript_summary <- txt("D_summary")

panelD_group_map <- c(
  "background_gene_without_conserved_domain" = "Background\ntranscripts",
  "nonterminal_conserved_domain" = "Non-terminal\ndomain transcripts",
  "terminal_conserved_domain" = "Terminal-domain\ntranscripts"
)
if (!identical(names(panelD_group_map), SELECT$D_groups)) stop("panel D group map out of sync")

panelD_group_levels <- c(
  "Background\ntranscripts",
  "Non-terminal\ndomain transcripts",
  "Terminal-domain\ntranscripts"
)

panelD_plot_data <- transcript_groups %>%
  filter(
    !is_single_exon_transcript,
    transcript_domain_group %in% names(panelD_group_map),
    !is.na(terminal_exon_length_bp),
    terminal_exon_length_bp > 0
  ) %>%
  mutate(
    group_label = unname(panelD_group_map[transcript_domain_group]),
    group_label = factor(group_label, levels = panelD_group_levels)
  )

panelD_summary_rows <- transcript_summary %>%
  filter(
    level == SELECT$D_summary[["level"]],
    scope == SELECT$D_summary[["scope"]],
    metric == SELECT$D_summary[["metric"]],
    group %in% names(panelD_group_map)
  )
if (nrow(panelD_summary_rows) != 3 || anyDuplicated(panelD_summary_rows$group)) stop("panel D needs one transcript-level summary row per group")
panelD_medians <- panelD_summary_rows %>%
  mutate(
    median = as.numeric(median),
    group_label = unname(panelD_group_map[group]),
    group_label = factor(group_label, levels = panelD_group_levels),
    display_label = sprintf("median %s bp", comma(median)),
    label_x = 15000
  )

panelD_counts <- panelD_plot_data %>%
  count(group_label, name = "n") %>%
  mutate(
    axis_label = paste0(as.character(group_label), "\n(n = ", comma(n), ")")
  )
chk <- panelD_counts %>% mutate(group = names(panelD_group_map)[match(as.character(group_label), panelD_group_map)]) %>%
  left_join(select(panelD_summary_rows, group, summary_n = n), by = "group")
if (nrow(chk) != 3 || any(as.character(chk$n) != chk$summary_n)) stop("panel D: plotted transcripts per group disagree with the order-34 summary n")

panelD_axis_labels <- setNames(
  panelD_counts$axis_label,
  as.character(panelD_counts$group_label)
)

panelD_colours <- c(
  "Background\ntranscripts" = COL$background,
  "Non-terminal\ndomain transcripts" = COL$nonterminal,
  "Terminal-domain\ntranscripts" = COL$terminal
)

pD <- ggplot(
  panelD_plot_data,
  aes(x = terminal_exon_length_bp, y = group_label, fill = group_label)
) +
  geom_boxplot(
    width = 0.62,
    outlier.shape = NA,
    linewidth = 0.6
  ) +
  geom_text(
    data = panelD_medians,
    aes(x = label_x, y = group_label, label = display_label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3.7,
    colour = "#222222"
  ) +
  scale_x_log10(
    breaks = c(100, 1000, 10000),
    labels = c("100", "1,000", "10,000"),
    expand = expansion(mult = c(0.06, 0.22))
  ) +
  scale_y_discrete(labels = panelD_axis_labels) +
  scale_fill_manual(values = panelD_colours, drop = FALSE) +
  labs(
    x = "Terminal-exon length (bp; log scale)",
    y = NULL
  ) +
  coord_cartesian(
    xlim = c(40, 20000),
    clip = "off"
  ) +
  theme_figure() +
  theme(
    legend.position = "none",
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text = element_text(size = 12.2),
    axis.text.y = element_text(margin = margin(r = 8)),
    axis.title.x = element_text(size = 11.2),
    plot.margin = margin(10, 60, 10, 28)
  )

# ==============================================================================
# Assemble (historical patchwork layout); panel E is added manually (MC04)
# ==============================================================================
figure_design <- "
AB
##
CD
"

figure2 <- pA + pB + pC + pD +
  plot_layout(
    design = figure_design,
    widths = c(1.00, 1.12),
    heights = c(1.00, CFG$row_gap, 1.00)
  ) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 15, colour = "#111111"),
    plot.tag.position = c(0, 1)
  )

# ------------------------------------------------------------------------------
# Output (both directories staged, then published)
# ------------------------------------------------------------------------------
stage_dir <- function(final) {
  parent <- dirname(normalizePath(final, mustWork = FALSE)); dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  s <- tempfile(pattern = ".stage_figure2.", tmpdir = parent); dir.create(s); s
}
ST_SRC <- stage_dir(SRC_OUT); ST_FIG <- stage_dir(FIG_OUT)
options(error = function() { unlink(c(ST_SRC, ST_FIG), recursive = TRUE); quit(save = "no", status = 1) })
ggsave(filename = file.path(ST_FIG, paste0(CFG$output_stem, ".svg")), plot = figure2, width = CFG$width_in, height = CFG$height_in,
       units = "in", device = svglite::svglite, bg = "white", limitsize = FALSE)
ggsave(filename = file.path(ST_FIG, paste0(CFG$output_stem, ".png")), plot = figure2, width = CFG$width_in, height = CFG$height_in,
       units = "in", dpi = CFG$dpi, device = ragg::agg_png, bg = "white", limitsize = FALSE)

# display statistics that ggplot2 itself draws in panel D (stat_boxplot on log10 bp, back-transformed); never used for text
box <- layer_data(pD, 1L)                               # horizontal boxplot: ggplot2 names the statistics x*
if (!all(c("xlower", "xmiddle", "xupper", "xmin", "xmax") %in% names(box))) stop("panel D: boxplot statistics not found in the layer data")
bt <- function(v) round(10^v, 3)
box_stats <- tibble(group_label = levels(panelD_plot_data$group_label)[box$group],
                    box_lower_bp = bt(box$xlower), box_median_bp = bt(box$xmiddle), box_upper_bp = bt(box$xupper),
                    whisker_low_bp = bt(box$xmin), whisker_high_bp = bt(box$xmax))
if (nrow(box_stats) != 3) stop("panel D: expected three boxplots")

src_a <- a_txt %>%
  mutate(panel_order = match(category, panelA_order), display_label = unname(panelA_labels[category]),
         denominator_herv_internal_regions = N_LOCI, flag_column = paste0(unname(SELECT$A_flag[category]), " & ", unname(SELECT$A_class[category])),
         value_label = panelA_plot_data$value_label[match(category, as.character(panelA_plot_data$category))],
         fill_colour = unname(panelA_colours[category]), mutually_exclusive = FALSE) %>%
  arrange(panel_order) %>%
  select(panel_order, category, display_label, n, percent, denominator_herv_internal_regions, flag_column, mutually_exclusive, value_label, fill_colour)
src_b <- b_txt %>%
  mutate(panel_row = match(analysis, panelB_analyses),
         display_label = c("All conserved domains", "Multi-exon sensitivity")[panel_row]) %>%
  left_join(select(panelB_data, analysis, observed_percent, expected_percent, lower_percent, upper_percent, observed_label, expected_label, p_label), by = "analysis") %>%
  arrange(panel_row) %>%
  select(panel_row, analysis, scope, scope_role, display_label, n_domains, observed_terminal_domains, observed_fraction, expected_terminal_domains,
         expected_fraction, exact_95pct_lower, exact_95pct_upper, exact_p_enrichment, observed_percent, expected_percent, lower_percent, upper_percent,
         observed_label, expected_label, p_label)
src_c <- bind_rows(mutate(cd$cont, level = "domain", display_label = "All domain predictions", source = "conservation_comparison"),
                   mutate(cl$cont, level = "locus", display_label = "One representative per HERV locus", source = "locus_representative")) %>%
  left_join(bind_rows(cd$fisher, cl$fisher) %>% select(analysis, coverage_threshold, odds_ratio, odds_ratio_95ci_lower, odds_ratio_95ci_upper, fisher_p_two_sided), by = "analysis") %>%
  mutate(panel_row = match(level, c("domain", "locus")),
         point_label = panelC_points$point_label[match(paste(analysis, conservation_group), paste(panelC_points$analysis, as.character(panelC_points$conservation_group)))],
         stats_label = gsub("\n", " | ", panelC_stats$stats_label[panel_row])) %>%
  arrange(panel_row, match(conservation_group, c("lower_coverage", "conserved"))) %>%
  select(panel_row, level, source, analysis, display_label, conservation_group, coverage_threshold, terminal, nonterminal, total, terminal_percent,
         odds_ratio, odds_ratio_95ci_lower, odds_ratio_95ci_upper, fisher_p_two_sided, point_label, stats_label)
src_d <- panelD_summary_rows %>%
  mutate(panel_row = match(group, names(panelD_group_map)), group_label = unname(panelD_group_map[group])) %>%
  left_join(box_stats, by = "group_label") %>%
  mutate(plotted_transcripts = as.integer(n), display_label = panelD_medians$display_label[match(group, panelD_medians$group)],
         axis_label = unname(panelD_axis_labels[group_label]),
         row_filter = "multi-exon (is_single_exon_transcript FALSE), terminal_exon_length_bp > 0, one of the three groups",
         group_label = gsub("\n", " ", group_label), axis_label = gsub("\n", " ", axis_label)) %>%
  arrange(panel_row) %>%
  select(panel_row, group, group_label, level, scope, metric, n, plotted_transcripts, median, q1, q3, display_label, axis_label,
         box_lower_bp, box_median_bp, box_upper_bp, whisker_low_bp, whisker_high_bp, row_filter)
contract <- tribble(
  ~panel, ~key, ~upstream_module, ~effective_order, ~row_unit, ~rows_used, ~display_transformation,
  "A", "A_counts", "03_transcript_context", "19 (D022; D024)", "gene-overlap category (HERV internal regions, same strand)", "4 of 4",
       paste0("bar length = percent (n / ", comma(N_LOCI), " loci); label sprintf('%.1f%%  (n = %s)'); non-exclusive categories"),
  "A", "A_loci", "03_transcript_context", "19 (D022; D024)", "HERV internal region (locus)", "consistency only",
       "row count = denominator; sum(is_X & overlap_Y) must equal n; not drawn",
  "B", "B_null", "04_domains_terminal_exons", "31", "analysis x scope of the transcript-aware positional null (unit: unique conserved genomic domain)",
       "2 of 12 (primary scope; any-transcript, all and multi-exon)", "percent = 100 x fraction; grey bar = exact 95% null range / n_domains; P = exact enrichment P",
  "B", "B_counts", "04_domains_terminal_exons", "30", "domain-centred count metric (conserved profile)", "consistency only",
       paste0(paste(cnt_n[SELECT$B_n], cnt_n[SELECT$B_analyses], sep = " / ", collapse = " and "), " must equal the null summary; not drawn"),
  "C", "C_cont", "04_domains_terminal_exons", "32", "coverage group x analysis (unit: exon-overlapping domain prediction)", "2 of 4 (all qualifying domains)",
       "point = terminal_percent; label sprintf('%.1f%% (%d/%d)')",
  "C", "C_fisher", "04_domains_terminal_exons", "32", "analysis (Fisher 2x2)", "1 of 2", "label OR %.2f, P format_p(digits = 2)",
  "C", "C_lcont", "04_domains_terminal_exons", "33", "coverage group x analysis (unit: HERV locus, one representative domain)", "2 of 4 (all transcripts)",
       "point = terminal_percent; label sprintf('%.1f%% (%d/%d)')",
  "C", "C_lfisher", "04_domains_terminal_exons", "33", "analysis (Fisher 2x2)", "1 of 2", "label OR %.2f, P format_p(digits = 2)",
  "D", "D_groups", "04_domains_terminal_exons", "34", "GENCODE transcript", paste0(nrow(panelD_plot_data), " of ", nrow(transcript_groups)),
       "boxplot of log10 terminal-exon length (ggplot2 stat_boxplot, outliers not drawn); x axis log10, view 40-20,000 bp",
  "D", "D_summary", "04_domains_terminal_exons", "34", "level x scope x group x metric", "3 of 48 (transcript, all_gene_types, terminal_exon_length_bp)",
       "median label sprintf('median %s bp'); n must equal the plotted transcripts"
) %>% mutate(upstream_file = vapply(key, function(k) IN[[k]]$rel, ""), sha256 = vapply(key, function(k) IN[[k]]$sha, "")) %>%
  select(panel, upstream_module, effective_order, upstream_file, sha256, row_unit, rows_used, display_transformation)
write_tsv(src_a, file.path(ST_SRC, "figure2A_source_data.tsv"))
write_tsv(src_b, file.path(ST_SRC, "figure2B_source_data.tsv"))
write_tsv(src_c, file.path(ST_SRC, "figure2C_source_data.tsv"))
write_tsv(src_d, file.path(ST_SRC, "figure2D_source_data.tsv"))
write_tsv(contract, file.path(ST_SRC, "figure2_source_contracts.tsv"))
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
message("[OK] Figure 2A: ", paste(src_a$category, src_a$n, collapse = "; "), " of ", N_LOCI, " loci; 2B: ",
        paste(src_b$observed_terminal_domains, src_b$n_domains, sep = "/", collapse = ", "), "; 2C: OR ",
        paste(unique(src_c$odds_ratio), collapse = ", "), "; 2D: n ", paste(src_d$plotted_transcripts, collapse = " / "))
message("[OK] wrote ", SRC_OUT, " and ", FIG_OUT)
