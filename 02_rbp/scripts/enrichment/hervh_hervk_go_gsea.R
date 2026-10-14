#!/usr/bin/env Rscript
# =============================================================================
# 02_rbp / enrichment : HERVH/HERVK profiles, GO ORA and pre-ranked GSEA (Stage 10)
# =============================================================================
# Refactor of scripts/rbp/10_hervh_hervk_composition_profiles_go_gsea.R.
# All scientific definitions are reproduced exactly; only paths, input
# resolution and validation changed.
#
#   Stage 9 subfamily x RBP enrichment tables
#     -> HERVH/HERVK profile classes
#     -> Figure 1A source data
#     -> GO BP over-representation (HERVH-enriched, HERVK-enriched)
#     -> pre-ranked GSEA, HERVH vs HERVK
#
# THIS STAGE RECOMPUTES NO STATISTICS FROM RAW FIMO. Stage 9 is the
# authoritative statistical input; this stage classifies, contrasts and
# functionally annotates what Stage 9 produced.
#
# SIGNIFICANCE FIELD: p_adj_global -- BH across all 9,672 Stage 9 subfamily x RBP
# tests. This is the same criterion Stage 9 used for its own `status` column and
# the one behind the manuscript's 20/40 and 15/37 counts. The
# p_adj_within_subfamily column is carried through but never decides anything.
#
# EFFECT SIZE: log2_observed_expected, from Stage 9.
#
# FIGURES ARE NOT RENDERED HERE. The release architecture keeps publication
# rendering in 08_figures_supplement and science modules emit source data only.
# Figure1A_heatmap_values.tsv is the canonical Figure 1A source; the historical
# script also wrote PDF/PNG/SVG, which module 08 will regenerate from this table.
#
# Usage:
#   hervh_hervk_go_gsea.R --stage9-dir D --out-dir D
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(clusterProfiler)
  library(fgsea)
  library(org.Hs.eg.db)
  library(GO.db)
  library(AnnotationDbi)
})

# ------------------------------------------------------------------ arguments
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  hit <- which(args == flag)
  if (length(hit) == 0) stop("missing required argument: ", flag)
  args[[hit + 1L]]
}

STAGE9_DIR <- get_arg("--stage9-dir")
OUT_DIR    <- get_arg("--out-dir")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Historical CFG, unchanged.
CFG <- list(
  fdr_cutoff = 0.05,
  heatmap_cap = 2,
  go_ontology = "BP",
  go_p_adjust_method = "BH",
  go_pvalue_cutoff = 1,        # keep every tested term; filter afterwards
  go_qvalue_cutoff = 1,
  go_min_gs_size = 10,
  go_max_gs_size = 500,
  go_significance_cutoff = 0.05,
  gsea_min_size = 5,
  gsea_max_size = 50,
  gsea_significance_cutoff = 0.05,
  gsea_seed = 260730           # fgseaMultilevel is stochastic; this pins it
)

# ---------------------------------------------------- package version guard
# org.Hs.eg.db and GO.db pin biological content, not just code: a different
# version changes GO term membership with byte-identical code. Refuse to run on
# anything but the historical versions.
required_versions <- c(
  clusterProfiler = "4.10.1", fgsea = "1.28.0",
  org.Hs.eg.db = "3.18.0", GO.db = "3.18.0", AnnotationDbi = "1.64.1"
)
for (pkg in names(required_versions)) {
  got <- as.character(packageVersion(pkg))
  if (!identical(got, required_versions[[pkg]])) {
    stop("package version mismatch: ", pkg, " is ", got,
         " but the release pins ", required_versions[[pkg]],
         ". GO/GSEA results are annotation-version dependent; refusing to run.")
  }
}

# ------------------------------------------------------------------ helpers
classify_status <- function(fdr, effect, cutoff = 0.05) {
  case_when(
    !is.na(fdr) & !is.na(effect) & fdr < cutoff & effect > 0 ~ "enriched",
    !is.na(fdr) & !is.na(effect) & fdr < cutoff & effect < 0 ~ "depleted",
    TRUE ~ "not_significant"
  )
}
make_safe_filename <- function(x) str_replace_all(x, "[^A-Za-z0-9._-]", "_")
safe_ratio <- function(numerator, denominator) {
  ifelse(is.na(denominator) | denominator == 0, NA_real_, numerator / denominator)
}
write_lines <- function(x, path) writeLines(enc2utf8(x), con = path, useBytes = TRUE)
fmt_int   <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_num   <- function(x, digits = 2) formatC(x, format = "f", digits = digits, big.mark = ",")
fmt_ratio <- function(x, digits = 3) formatC(x, format = "f", digits = digits)

# --------------------------------------------------------------------- input
p_hervh_hervk <- file.path(STAGE9_DIR,
  "rbp_HERVH_HERVK_composition_controlled_enrichment.tsv")
p_motif <- file.path(STAGE9_DIR, "rbp_motif_library_summary.tsv")
p_subfam <- file.path(STAGE9_DIR, "rbp_subfamily_length_gc_density_summary.tsv")

for (p in c(p_hervh_hervk, p_motif, p_subfam)) {
  if (!file.exists(p)) stop("required Stage 9 input missing: ", p)
}

hervh_hervk <- read_tsv(p_hervh_hervk, show_col_types = FALSE, progress = FALSE)
motif_library <- read_tsv(p_motif, show_col_types = FALSE, progress = FALSE)
subfamily_summary <- read_tsv(p_subfam, show_col_types = FALSE, progress = FALSE)

all_tested_rbps <- sort(unique(motif_library$rbp))
if (length(all_tested_rbps) != 93) {
  stop("expected 93 tested RBPs, found ", length(all_tested_rbps))
}
if (!setequal(unique(hervh_hervk$subfamily), c("HERVH", "HERVK"))) {
  stop("the HERVH/HERVK table does not contain exactly those two subfamilies")
}

# Observed/expected ratio for the Figure 1A export, back-transformed from the
# Stage 9 effect size: 2^log2_observed_expected. It therefore CARRIES THE 0.5
# PSEUDOCOUNT that Stage 9 applied inside the log2. It is deliberately NOT the
# raw density ratio observed_hits_per_kb / mean_shuffled_hits_per_kb -- the two
# differ (e.g. ANKHD1/HERVH: 0.8167 vs 0.8130), and the historical Figure 1A
# source data uses the back-transformed form.
hervh_hervk <- hervh_hervk %>%
  mutate(observed_expected_ratio = 2^log2_observed_expected)

# Stage 9 already carries `status`; recompute it here and require agreement, so
# a change in either stage's rule is caught rather than silently absorbed.
recomputed <- classify_status(hervh_hervk$p_adj_global,
                              hervh_hervk$log2_observed_expected,
                              CFG$fdr_cutoff)
if (!identical(as.character(recomputed), as.character(hervh_hervk$status))) {
  stop("recomputed status disagrees with the Stage 9 status column")
}

status_counts <- hervh_hervk %>%
  count(subfamily, status, name = "n_rbps") %>%
  complete(subfamily = c("HERVH", "HERVK"),
           status = c("enriched", "depleted", "not_significant"),
           fill = list(n_rbps = 0L)) %>%
  arrange(match(subfamily, c("HERVH", "HERVK")),
          match(status, c("enriched", "depleted", "not_significant")))

write_tsv(status_counts, file.path(OUT_DIR, "manuscript_RBP_status_counts.tsv"))
message("status counts:")
print(status_counts)

# Full Stage 9 schema is carried through, not a reduced projection.
rbp_lists <- hervh_hervk %>%
  arrange(subfamily,
          factor(status, levels = c("enriched", "depleted", "not_significant")),
          p_adj_global, p_adj_within_subfamily,
          desc(abs(log2_observed_expected)), rbp)

write_tsv(rbp_lists, file.path(OUT_DIR, "HERVH_HERVK_all_tested_RBPs.tsv"))
write_tsv(rbp_lists %>% filter(status == "enriched"),
          file.path(OUT_DIR, "HERVH_HERVK_enriched_RBPs.tsv"))
write_tsv(rbp_lists %>% filter(status == "depleted"),
          file.path(OUT_DIR, "HERVH_HERVK_depleted_RBPs.tsv"))

# ----------------------------------------- reciprocal / shared profile classes
profile_matrix <- hervh_hervk %>%
  dplyr::select(rbp, subfamily, status, log2_observed_expected,
                p_adj_global, p_adj_within_subfamily) %>%
  pivot_wider(
    names_from = subfamily,
    values_from = c(status, log2_observed_expected,
                    p_adj_global, p_adj_within_subfamily),
    names_sep = "__"
  ) %>%
  mutate(
    profile_class = case_when(
      status__HERVH == "enriched"        & status__HERVK == "depleted"        ~ "HERVH_enriched_HERVK_depleted",
      status__HERVH == "depleted"        & status__HERVK == "enriched"        ~ "HERVH_depleted_HERVK_enriched",
      status__HERVH == "enriched"        & status__HERVK == "enriched"        ~ "enriched_in_both",
      status__HERVH == "depleted"        & status__HERVK == "depleted"        ~ "depleted_in_both",
      status__HERVH == "enriched"        & status__HERVK == "not_significant" ~ "HERVH_only_enriched",
      status__HERVH == "not_significant" & status__HERVK == "enriched"        ~ "HERVK_only_enriched",
      status__HERVH == "depleted"        & status__HERVK == "not_significant" ~ "HERVH_only_depleted",
      status__HERVH == "not_significant" & status__HERVK == "depleted"        ~ "HERVK_only_depleted",
      status__HERVH == "not_significant" & status__HERVK == "not_significant" ~ "not_significant_in_either",
      TRUE ~ "other"
    )
  ) %>%
  arrange(profile_class, rbp)

write_tsv(profile_matrix, file.path(OUT_DIR, "HERVH_HERVK_profile_classes.tsv"))

profile_class_counts <- profile_matrix %>%
  count(profile_class, name = "n_rbps") %>%
  arrange(desc(n_rbps), profile_class)
write_tsv(profile_class_counts,
          file.path(OUT_DIR, "HERVH_HERVK_profile_class_counts.tsv"))

# ------------------------------------------------------ Figure 1A source data
# Selection: every RBP passing global BH FDR < 0.05 in at least one family.
# Display: the actual log2 O/E in BOTH families, whether or not that family
# itself passes FDR -- so the quantitative comparison is preserved rather than
# selected from the narrative.
heatmap_long <- hervh_hervk %>%
  group_by(rbp) %>%
  filter(any(p_adj_global < CFG$fdr_cutoff, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    heatmap_value_uncapped = log2_observed_expected,
    heatmap_value_capped = pmax(pmin(log2_observed_expected, CFG$heatmap_cap),
                                -CFG$heatmap_cap),
    significant = !is.na(p_adj_global) & p_adj_global < CFG$fdr_cutoff
  )

heatmap_export <- heatmap_long %>%
  dplyr::select(rbp, subfamily, heatmap_value_uncapped, heatmap_value_capped,
                observed_expected_ratio, p_value, p_adj_global,
                p_adj_within_subfamily, significant, status) %>%
  arrange(rbp, match(subfamily, c("HERVH", "HERVK")))

write_tsv(heatmap_export, file.path(OUT_DIR, "Figure1A_heatmap_values.tsv"))

# Rendering parameters recorded so module 08 can reproduce the historical panel
# exactly without re-deriving them from this script.
write_tsv(
  tibble(
    parameter = c("colour_low", "colour_mid", "colour_high", "cap_abs",
                  "breaks_n", "cluster_rows", "cluster_cols", "na_colour",
                  "fontsize", "fontsize_row", "fontsize_col",
                  "width_in", "height_in", "png_dpi"),
    value = c("#2166ac", "white", "#b2182b", as.character(CFG$heatmap_cap),
              "101", "TRUE", "FALSE", "grey90", "11", "7", "15", "6", "8", "300")
  ),
  file.path(OUT_DIR, "Figure1A_render_parameters.tsv")
)

# =============================================================================
# GO BP over-representation
# =============================================================================
# Unit of analysis is the RBP, not the HERV locus.
# Selected set : composition-controlled ENRICHED RBPs of that family.
# Universe     : the tested motif-library RBPs mappable in org.Hs.eg.db.
# Depleted RBPs are NOT submitted to ORA -- historical behaviour.
hervh_enriched_rbps <- hervh_hervk %>%
  filter(subfamily == "HERVH", p_adj_global < CFG$fdr_cutoff,
         log2_observed_expected > 0) %>%
  arrange(p_adj_global, desc(log2_observed_expected)) %>%
  pull(rbp) %>% unique()

hervk_enriched_rbps <- hervh_hervk %>%
  filter(subfamily == "HERVK", p_adj_global < CFG$fdr_cutoff,
         log2_observed_expected > 0) %>%
  arrange(p_adj_global, desc(log2_observed_expected)) %>%
  pull(rbp) %>% unique()

valid_orgdb_symbols <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL")

mapping_table <- tibble(
  rbp = all_tested_rbps,
  represented_in_motif_library = TRUE,
  mappable_in_org_Hs_eg_db = rbp %in% valid_orgdb_symbols,
  selected_for_HERVH_GO = rbp %in% hervh_enriched_rbps,
  selected_for_HERVK_GO = rbp %in% hervk_enriched_rbps
)
write_tsv(mapping_table, file.path(OUT_DIR, "GO_RBP_universe_and_mapping.tsv"))

effective_universe <- mapping_table %>%
  filter(mappable_in_org_Hs_eg_db) %>% pull(rbp)

hervh_go_input <- intersect(hervh_enriched_rbps, effective_universe)
hervk_go_input <- intersect(hervk_enriched_rbps, effective_universe)

message("tested RBPs: ", length(all_tested_rbps),
        "; mappable universe: ", length(effective_universe),
        "; unmapped: ", length(setdiff(all_tested_rbps, effective_universe)))

write_tsv(tibble(rbp = all_tested_rbps),
          file.path(OUT_DIR, "GO_background_all_tested_RBPs.tsv"))
write_tsv(tibble(rbp = effective_universe),
          file.path(OUT_DIR, "GO_background_mappable_RBPs.tsv"))
write_tsv(tibble(rbp = hervh_enriched_rbps),
          file.path(OUT_DIR, "GO_input_HERVH_enriched_RBPs.tsv"))
write_tsv(tibble(rbp = hervk_enriched_rbps),
          file.path(OUT_DIR, "GO_input_HERVK_enriched_RBPs.tsv"))

run_go_bp <- function(input_genes, label) {
  label_safe <- make_safe_filename(label)
  if (length(input_genes) == 0) {
    write_lines(paste0("No mappable enriched RBPs were available for ", label, "."),
                file.path(OUT_DIR, paste0("GO_", label_safe, "_NO_INPUT.txt")))
    return(NULL)
  }

  ego <- enrichGO(
    gene = input_genes,
    universe = effective_universe,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",            # symbols directly; no Entrez conversion
    ont = CFG$go_ontology,
    pAdjustMethod = CFG$go_p_adjust_method,
    pvalueCutoff = CFG$go_pvalue_cutoff,
    qvalueCutoff = CFG$go_qvalue_cutoff,
    minGSSize = CFG$go_min_gs_size,
    maxGSSize = CFG$go_max_gs_size,
    readable = FALSE
  )

  full_result <- as.data.frame(ego) %>% as_tibble() %>%
    arrange(p.adjust, pvalue, desc(Count), ID)

  significant_result <- full_result %>%
    filter(!is.na(p.adjust), p.adjust < CFG$go_significance_cutoff)

  write_tsv(full_result,
            file.path(OUT_DIR, paste0("GO_", label_safe, "_BP_complete.tsv")))
  write_tsv(significant_result,
            file.path(OUT_DIR, paste0("GO_", label_safe, "_BP_significant_FDR05.tsv")))

  if (nrow(significant_result) == 0) {
    write_lines(
      c(paste0("No GO Biological Process terms passed BH FDR < ",
               CFG$go_significance_cutoff, " for ", label, "."),
        paste0("Input RBPs: ", length(input_genes),
               "; effective mapped universe: ", length(effective_universe), "."),
        paste0("Parameters: minGSSize=", CFG$go_min_gs_size,
               ", maxGSSize=", CFG$go_max_gs_size, ".")),
      file.path(OUT_DIR, paste0("GO_", label_safe, "_NO_SIGNIFICANT_TERMS.txt"))
    )
  }
  list(full = full_result, significant = significant_result)
}

go_hervh <- run_go_bp(hervh_go_input, "HERVH_enriched_RBPs")
go_hervk <- run_go_bp(hervk_go_input, "HERVK_enriched_RBPs")

# =============================================================================
# Pre-ranked GSEA: HERVH versus HERVK
# =============================================================================
# Ranking statistic: log2(O/E)_HERVH - log2(O/E)_HERVK, over EVERY tested RBP
# (not just the significant ones). Effect-size based rather than p-value based,
# because HERVH and HERVK contain very different numbers of loci.
gsea_rank_table <- hervh_hervk %>%
  dplyr::select(rbp, subfamily, log2_observed_expected) %>%
  tidyr::pivot_wider(names_from = subfamily, values_from = log2_observed_expected) %>%
  dplyr::filter(!is.na(HERVH), !is.na(HERVK)) %>%
  dplyr::mutate(contrast_HERVH_minus_HERVK = HERVH - HERVK) %>%
  dplyr::arrange(dplyr::desc(contrast_HERVH_minus_HERVK), rbp)

if (nrow(gsea_rank_table) != length(all_tested_rbps)) {
  stop("the pre-ranked analysis expected ", length(all_tested_rbps),
       " tested RBPs but obtained ", nrow(gsea_rank_table), " contrasts.")
}

gsea_ranks_raw <- stats::setNames(gsea_rank_table$contrast_HERVH_minus_HERVK,
                                  gsea_rank_table$rbp)

# Deterministic <=1e-10 offset, purely to resolve exact ties reproducibly.
tie_offset <- (rank(names(gsea_ranks_raw), ties.method = "first") -
               mean(rank(names(gsea_ranks_raw), ties.method = "first"))) * 1e-12
gsea_ranks <- sort(gsea_ranks_raw + tie_offset, decreasing = TRUE)

gsea_rank_table <- gsea_rank_table %>%
  dplyr::mutate(rank_statistic_used = unname(gsea_ranks[rbp]),
                rank_position = rank(-rank_statistic_used, ties.method = "first")) %>%
  dplyr::arrange(rank_position)

write_tsv(gsea_rank_table, file.path(OUT_DIR, "GSEA_HERVH_vs_HERVK_RBP_ranking.tsv"))

# Gene sets are GO BP memberships from org.Hs.eg.db, restricted to the tested
# RBP universe. No external MSigDB collection is involved.
go_bp_mapping <- AnnotationDbi::select(
    org.Hs.eg.db, keys = names(gsea_ranks),
    keytype = "SYMBOL", columns = c("GO", "ONTOLOGY")) %>%
  tibble::as_tibble() %>%
  dplyr::filter(ONTOLOGY == "BP", !is.na(GO), SYMBOL %in% names(gsea_ranks)) %>%
  dplyr::distinct(GO, SYMBOL)

gsea_pathways <- split(go_bp_mapping$SYMBOL, go_bp_mapping$GO)
gsea_pathways <- lapply(gsea_pathways,
                        function(x) sort(unique(intersect(x, names(gsea_ranks)))))
gsea_pathway_sizes <- lengths(gsea_pathways)
gsea_pathways <- gsea_pathways[gsea_pathway_sizes >= CFG$gsea_min_size &
                               gsea_pathway_sizes <= CFG$gsea_max_size]

if (length(gsea_pathways) == 0) {
  stop("No GO BP gene sets remained after restricting to the tested RBP universe ",
       "and applying size limits ", CFG$gsea_min_size, "-", CFG$gsea_max_size, ".")
}

gsea_go_terms <- AnnotationDbi::select(GO.db, keys = names(gsea_pathways),
                                       keytype = "GOID", columns = "TERM") %>%
  tibble::as_tibble() %>% dplyr::filter(!is.na(TERM)) %>%
  dplyr::distinct(GOID, .keep_all = TRUE)

gsea_pathway_membership <- tibble::tibble(
    pathway = names(gsea_pathways),
    pathway_size_in_tested_universe = lengths(gsea_pathways),
    members = vapply(gsea_pathways, paste, collapse = ";", FUN.VALUE = character(1))
  ) %>%
  dplyr::left_join(gsea_go_terms, by = c("pathway" = "GOID")) %>%
  dplyr::relocate(pathway, TERM, pathway_size_in_tested_universe, members) %>%
  dplyr::arrange(pathway)

write_tsv(gsea_pathway_membership,
          file.path(OUT_DIR, "GSEA_HERVH_vs_HERVK_GO_BP_pathway_membership.tsv"))

set.seed(CFG$gsea_seed)
gsea_result_raw <- fgsea::fgseaMultilevel(
  pathways = gsea_pathways,
  stats = gsea_ranks,
  minSize = CFG$gsea_min_size,
  maxSize = CFG$gsea_max_size,
  scoreType = "std",
  eps = 0
)

gsea_result <- gsea_result_raw %>%
  tibble::as_tibble() %>%
  dplyr::left_join(gsea_go_terms, by = c("pathway" = "GOID")) %>%
  dplyr::mutate(
    direction = dplyr::case_when(NES > 0 ~ "HERVH", NES < 0 ~ "HERVK",
                                 TRUE ~ "neutral"),
    leading_edge_size = lengths(leadingEdge),
    leading_edge_RBPs = vapply(leadingEdge, paste, collapse = ";",
                               FUN.VALUE = character(1))
  ) %>%
  dplyr::select(pathway, TERM, direction, pval, padj, log2err, ES, NES, size,
                leading_edge_size, leading_edge_RBPs) %>%
  dplyr::arrange(padj, pval, dplyr::desc(abs(NES)), pathway)

write_tsv(gsea_result, file.path(OUT_DIR, "GSEA_HERVH_vs_HERVK_GO_BP_complete.tsv"))

gsea_significant <- gsea_result %>%
  dplyr::filter(!is.na(padj), padj < CFG$gsea_significance_cutoff)
write_tsv(gsea_significant,
          file.path(OUT_DIR, "GSEA_HERVH_vs_HERVK_GO_BP_significant_FDR05.tsv"))

# When nothing survives adjustment, the historical analysis recorded a note
# rather than leaving a silently empty file. In the primary profile one term is
# significant so this branch does not fire; in the strict profile none is.
if (nrow(gsea_significant) == 0) {
  write_lines(
    c(paste0("No pre-ranked GO Biological Process terms passed BH FDR < ",
             CFG$gsea_significance_cutoff, "."),
      paste0("The complete result and the top nominal terms are provided for ",
             "transparent inspection, but nominal terms should not be presented as ",
             "formal category-level enrichment."),
      paste0("Ranking statistic: log2(observed/expected)_HERVH - ",
             "log2(observed/expected)_HERVK."),
      paste0("Pathway-size limits within the tested RBP universe: ",
             CFG$gsea_min_size, "-", CFG$gsea_max_size, ".")),
    file.path(OUT_DIR, "GSEA_HERVH_vs_HERVK_NO_SIGNIFICANT_TERMS.txt")
  )
}

gsea_top_nominal <- gsea_result %>%
  dplyr::filter(!is.na(pval)) %>%
  dplyr::arrange(pval, dplyr::desc(abs(NES)), pathway) %>%
  head(25)
write_tsv(gsea_top_nominal,
          file.path(OUT_DIR, "GSEA_HERVH_vs_HERVK_GO_BP_top25_nominal.tsv"))

# ---------------------------------------------------------- manuscript support
write_tsv(motif_library %>% arrange(desc(n_motifs), rbp),
          file.path(OUT_DIR, "motif_library_complete.tsv"))

write_tsv(
  tibble(n_unique_rbps = n_distinct(motif_library$rbp),
         n_total_motifs = sum(motif_library$n_motifs, na.rm = TRUE),
         minimum_motifs_per_rbp = min(motif_library$n_motifs, na.rm = TRUE),
         maximum_motifs_per_rbp = max(motif_library$n_motifs, na.rm = TRUE)),
  file.path(OUT_DIR, "manuscript_motif_library_scope.tsv")
)

hervh_hervk_global <- subfamily_summary %>%
  filter(subfamily %in% c("HERVH", "HERVK")) %>%
  mutate(
    calculated_density_ratio = safe_ratio(observed_hits_per_kb,
                                          mean_shuffled_hits_per_kb),
    calculated_total_ratio = safe_ratio(observed_total_hits,
                                        mean_shuffled_total_hits)
  ) %>%
  arrange(match(subfamily, c("HERVH", "HERVK")))

if (nrow(hervh_hervk_global) != 2) {
  stop("Could not retrieve exactly one HERVH row and one HERVK row.")
}
# Independent cross-check of Stage 9's ratio against a recomputation here.
if (nrow(hervh_hervk_global %>%
         mutate(d = abs(calculated_density_ratio - global_density_ratio)) %>%
         filter(d > 1e-6)) > 0) {
  stop("global_density_ratio does not match the observed/shuffled density ratio.")
}
write_tsv(hervh_hervk_global,
          file.path(OUT_DIR, "manuscript_HERVH_HERVK_global_summary.tsv"))

# Descriptive profile contrast. NOT a formal interaction test.
effect_wide <- hervh_hervk %>%
  dplyr::select(rbp, subfamily, log2_observed_expected) %>%
  pivot_wider(names_from = subfamily, values_from = log2_observed_expected)

cap <- function(x) pmax(pmin(x, CFG$heatmap_cap), -CFG$heatmap_cap)
profile_contrast <- tibble(
  comparison = "HERVH_vs_HERVK",
  n_tested_rbps = nrow(effect_wide),
  pearson_correlation = cor(effect_wide$HERVH, effect_wide$HERVK,
                            method = "pearson", use = "complete.obs"),
  spearman_correlation = cor(effect_wide$HERVH, effect_wide$HERVK,
                             method = "spearman", use = "complete.obs"),
  euclidean_distance_uncapped =
    sqrt(sum((effect_wide$HERVH - effect_wide$HERVK)^2)),
  euclidean_distance_capped_at_heatmap_limit =
    sqrt(sum((cap(effect_wide$HERVH) - cap(effect_wide$HERVK))^2))
)
write_tsv(profile_contrast,
          file.path(OUT_DIR, "HERVH_HERVK_descriptive_profile_contrast.tsv"))

# ------------------------------------------------- manuscript support report
# A narrative restatement of every manuscript-facing number produced above,
# plus the interpretive guardrails the historical analysis recorded alongside
# them. Reproduced verbatim because those guardrails constrain how these
# results may be worded.
get_global_row <- function(subfamily_name) {
  out <- hervh_hervk_global[hervh_hervk_global$subfamily == subfamily_name, , drop = FALSE]
  if (nrow(out) != 1L) {
    stop("Expected exactly one global-summary row for ", subfamily_name,
         ", but found ", nrow(out), ".")
  }
  out
}
get_status_n <- function(sf, st) {
  status_counts %>% filter(subfamily == sf, status == st) %>% pull(n_rbps)
}

hervh_global <- get_global_row("HERVH")
hervk_global <- get_global_row("HERVK")
motif_scope <- tibble(
  n_unique_rbps = n_distinct(motif_library$rbp),
  n_total_motifs = sum(motif_library$n_motifs, na.rm = TRUE),
  minimum_motifs_per_rbp = min(motif_library$n_motifs, na.rm = TRUE),
  maximum_motifs_per_rbp = max(motif_library$n_motifs, na.rm = TRUE)
)

n_reciprocal_hervh <- profile_matrix %>% filter(profile_class == "HERVH_enriched_HERVK_depleted") %>% nrow()
n_reciprocal_hervk <- profile_matrix %>% filter(profile_class == "HERVH_depleted_HERVK_enriched") %>% nrow()
n_shared_enriched  <- profile_matrix %>% filter(profile_class == "enriched_in_both") %>% nrow()

hervh_go_n_sig <- if (is.null(go_hervh)) 0L else nrow(go_hervh$significant)
hervk_go_n_sig <- if (is.null(go_hervk)) 0L else nrow(go_hervk$significant)

gsea_n_sig   <- nrow(gsea_significant)
gsea_n_hervh <- sum(gsea_significant$NES > 0, na.rm = TRUE)
gsea_n_hervk <- sum(gsea_significant$NES < 0, na.rm = TRUE)
gsea_hervh_indices <- which(gsea_significant$NES > 0)
gsea_hervk_indices <- which(gsea_significant$NES < 0)
gsea_top_hervh <- if (length(gsea_hervh_indices) > 0) gsea_significant[gsea_hervh_indices[[1]], , drop = FALSE] else gsea_significant[0, , drop = FALSE]
gsea_top_hervk <- if (length(gsea_hervk_indices) > 0) gsea_significant[gsea_hervk_indices[[1]], , drop = FALSE] else gsea_significant[0, , drop = FALSE]

top_label <- function(tb) {
  if (nrow(tb) == 0) return("none")
  nm <- if (!is.na(tb$TERM[[1]])) tb$TERM[[1]] else tb$pathway[[1]]
  paste0(nm, " (NES=", fmt_num(tb$NES[[1]], 2), ", FDR=",
         formatC(tb$padj[[1]], format = "g", digits = 3), ")")
}
gsea_top_hervh_label <- top_label(gsea_top_hervh)
gsea_top_hervk_label <- top_label(gsea_top_hervk)

support_report <- c(
  "COMPOSITION-CONTROLLED HERVH/HERVK RBP ANALYSIS",
  "================================================",
  "",
  paste0("Motif library: ", motif_scope$n_total_motifs, " motifs representing ",
         motif_scope$n_unique_rbps, " RBPs (", motif_scope$minimum_motifs_per_rbp,
         "-", motif_scope$maximum_motifs_per_rbp, " motifs per RBP)."),
  "",
  "GLOBAL SUBFAMILY BURDEN",
  paste0("HERVH: ", fmt_int(hervh_global$observed_total_hits), " observed hits versus ",
         fmt_int(hervh_global$mean_shuffled_total_hits), " shuffled expectation; ",
         fmt_num(hervh_global$observed_hits_per_kb, 2), " versus ",
         fmt_num(hervh_global$mean_shuffled_hits_per_kb, 2),
         " hits/kb; observed/expected ratio = ",
         fmt_ratio(hervh_global$global_density_ratio, 3), "."),
  paste0("HERVK: ", fmt_int(hervk_global$observed_total_hits), " observed hits versus ",
         fmt_int(hervk_global$mean_shuffled_total_hits), " shuffled expectation; ",
         fmt_num(hervk_global$observed_hits_per_kb, 2), " versus ",
         fmt_num(hervk_global$mean_shuffled_hits_per_kb, 2),
         " hits/kb; observed/expected ratio = ",
         fmt_ratio(hervk_global$global_density_ratio, 3), "."),
  "",
  "RBP-SPECIFIC RESULTS",
  paste0("HERVH: ", get_status_n("HERVH", "enriched"), " enriched and ",
         get_status_n("HERVH", "depleted"), " depleted RBPs at global BH FDR < ",
         CFG$fdr_cutoff, "."),
  paste0("HERVK: ", get_status_n("HERVK", "enriched"), " enriched and ",
         get_status_n("HERVK", "depleted"), " depleted RBPs at global BH FDR < ",
         CFG$fdr_cutoff, "."),
  paste0("Reciprocal classes: ", n_reciprocal_hervh,
         " HERVH-enriched/HERVK-depleted; ", n_reciprocal_hervk,
         " HERVH-depleted/HERVK-enriched."),
  paste0("Shared enrichment: ", n_shared_enriched, " RBPs enriched in both families."),
  "",
  "GO BIOLOGICAL PROCESS ENRICHMENT",
  paste0("HERVH input: ", length(hervh_enriched_rbps), " enriched RBPs; ",
         length(hervh_go_input), " mapped to org.Hs.eg.db; ", hervh_go_n_sig,
         " significant BP terms at BH FDR < ", CFG$go_significance_cutoff, "."),
  paste0("HERVK input: ", length(hervk_enriched_rbps), " enriched RBPs; ",
         length(hervk_go_input), " mapped to org.Hs.eg.db; ", hervk_go_n_sig,
         " significant BP terms at BH FDR < ", CFG$go_significance_cutoff, "."),
  "",
  "PRE-RANKED HERVH-VERSUS-HERVK GO BP ANALYSIS",
  paste0("Ranking: log2(observed/expected)_HERVH - ",
         "log2(observed/expected)_HERVK across ", nrow(gsea_rank_table),
         " tested RBPs."),
  paste0("Tested GO BP sets after restricting membership to the motif-library ",
         "universe and applying size limits ", CFG$gsea_min_size, "-",
         CFG$gsea_max_size, ": ", nrow(gsea_result), "."),
  paste0("Significant terms at BH FDR < ", CFG$gsea_significance_cutoff, ": ",
         gsea_n_sig, " total; ", gsea_n_hervh, " HERVH-biased and ",
         gsea_n_hervk, " HERVK-biased."),
  paste0("Top significant HERVH-biased term: ", gsea_top_hervh_label, "."),
  paste0("Top significant HERVK-biased term: ", gsea_top_hervk_label, "."),
  "",
  "INTERPRETIVE GUARDRAILS",
  paste0("- The phrase 'strongest contrast across subfamilies' is not established ",
         "by this analysis and should not be used without a predefined, formal ",
         "cross-subfamily contrast statistic."),
  paste0("- The heatmap supports 'contrasting' or 'partially reciprocal' ",
         "quantitative profiles."),
  paste0("- Functional category claims should be based on the complete GO result ",
         "tables generated here, not on a manually selected list of RBPs."),
  paste0("- Nominal pre-ranked GO terms that do not pass adjusted significance ",
         "may be described only as exploratory trends, not as formal support for ",
         "a developmental enrichment claim."),
  paste0("- These are predicted motif enrichments, not evidence of physical RBP ",
         "binding.")
)

write_lines(support_report, file.path(OUT_DIR, "manuscript_support_report.txt"))

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

cat("\n=== PROFILE CLASS COUNTS ===\n"); print(profile_class_counts, n = Inf)
cat("\n=== GO ===\n")
cat("HERVH input:", length(hervh_go_input), " terms tested:",
    if (is.null(go_hervh)) 0 else nrow(go_hervh$full), "\n")
cat("HERVK input:", length(hervk_go_input), " terms tested:",
    if (is.null(go_hervk)) 0 else nrow(go_hervk$full), "\n")
cat("\n=== GSEA ===\n")
cat("ranked RBPs:", length(gsea_ranks), " pathways:", length(gsea_pathways),
    " significant:", nrow(gsea_significant), "\n")
cat("\nResults written to:\n", OUT_DIR, "\n")
