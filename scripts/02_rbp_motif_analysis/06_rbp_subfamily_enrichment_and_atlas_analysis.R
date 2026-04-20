# ==============================================================================
# 02_RBP_regulatory_landscape_and_enrichment.R
# ==============================================================================

# Purpose
# Characterize the RNA regulatory landscape of HERV internal regions by integrating:
#   (1) motif-level RBP binding predictions (FIMO)
#   (2) subfamily-level enrichment of RBP motifs
#   (3) locus-level structural and regulatory features from the HERV atlas
#
# Biological rationale
# RNA-binding proteins (RBPs) play central roles in post-transcriptional regulation,
# including splicing, mRNA stability, localization, and translation.
# HERV sequences embedded in transcripts may contribute regulatory modules
# that recruit RBPs in a subfamily-specific and locus-specific manner.
#
# Key idea
# HERV loci are treated as RNA regulatory elements whose functional potential is shaped by:
#   - motif abundance (RBP burden)
#   - motif diversity (unique RBPs per locus)
#   - subfamily-specific enrichment of RBP motifs
#   - integration with transcriptomic and structural context
#
# Pipeline overview
#   1. Load FIMO motif predictions and parse RBP annotations
#   2. Summarize RBP motif burden per subfamily:
#        - total motif counts
#        - top RBPs by abundance and statistical strength
#   3. Perform enrichment analysis (Fisher’s exact test):
#        - identify RBPs preferentially associated with each subfamily
#   4. Visualize enrichment patterns:
#        - heatmaps (e.g., HERVH vs HERVK comparison)
#   5. Integrate with HERV atlas:
#        - compute locus-level features (RBP burden, length, domains)
#        - summarize subfamily-level regulatory profiles
#   6. Explore locus-level relationships:
#        - RBP burden vs locus length
#        - identification of highly RBP-loaded loci
#
# Outputs
#   - RBP counts and rankings per subfamily
#   - RBP enrichment statistics (odds ratios, p-values, FDR)
#   - Heatmaps of subfamily-specific RBP enrichment
#   - Atlas summaries of regulatory features
#   - Scatter plots of RBP burden vs locus length
#
# Notes
#   - Motif hits are not collapsed here (burden reflects total motif signal)
#   - q-value filtering is optional (controlled by CFG$q_cut)
#   - Enrichment is computed using Fisher’s exact test on motif counts
#   - Interpretation:
#        odds_ratio > 1 → enrichment
#        odds_ratio < 1 → depletion
#
# Position in pipeline
#   This script corresponds to the "RBP / RNA regulation" module.
#   It complements:
#     - transcript-context analysis (structural integration)
#     - LTR regulatory architecture (IFN motif clustering)
#
#   Together, these layers define the RNA-level functional potential of HERVs.
# ==============================================================================

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(openxlsx)
})

# ============================================================
# Config
# ============================================================
CFG <- list(
  fimo_in    = "../../../data/fimo_sense_parsed.tsv",
  atlas_in   = "../../../data/internal_integrated.with_hervs_id.v2.tsv",
  out_dir    = "../../../results/RBP_subfamilies/",
  top_n      = 10,
  q_cut      = 1,     # set 0.05 to restrict to significant
  make_plots = TRUE
)

dir.create(CFG$out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Helpers
# ============================================================
as_num <- function(x) suppressWarnings(as.numeric(x))
as_int <- function(x) suppressWarnings(as.integer(x))
prev  <- function(x) mean(x, na.rm = TRUE)

safe_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  suppressWarnings(min(x, na.rm = TRUE))
}

safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  suppressWarnings(max(x, na.rm = TRUE))
}

infer_subfamily <- function(x) {
  # "HERVK14-int_00001" -> "HERVK14"
  # "ERVL-E-int_00015"  -> "ERVL-E"
  str_remove(x, "-int_.*$")
}

infer_class <- function(subfamily) {
  case_when(
    str_detect(subfamily, "^HERV") ~ "ERV",
    TRUE ~ "Other"
  )
}

# ============================================================
# 1) FIMO: load + parse
# ============================================================
load_fimo <- function(path) {
  fimo_raw <- read_tsv(path, show_col_types = FALSE)
  
  fimo_parsed <- fimo_raw %>%
    mutate(
      subfamily = infer_subfamily(sequence_name),
      rbp       = str_extract(motif_alt_id, "^[^|]+"),
      `p-value` = as_num(`p-value`),
      `q-value` = as_num(`q-value`),
      score     = as_num(score)
    )
  
  if (any(is.na(fimo_parsed$subfamily))) {
    warning("Some rows have NA subfamily. Check sequence_name format for '-int'.")
  }
  if (any(is.na(fimo_parsed$rbp))) {
    warning("Some rows have NA rbp. Check motif_alt_id format for '|'.")
  }
  
  list(raw = fimo_raw, parsed = fimo_parsed)
}

# ============================================================
# 2) FIMO summary tables
# ============================================================
summarise_fimo_by_subfamily <- function(fimo_parsed, q_cut, top_n) {
  fimo_sig <- fimo_parsed %>%
    filter(!is.na(subfamily), !is.na(rbp)) %>%
    filter(is.na(`q-value`) | `q-value` <= q_cut)
  
  counts <- fimo_sig %>%
    count(subfamily, rbp, name = "n_hits") %>%
    arrange(subfamily, desc(n_hits), rbp)
  
  top_by_count <- counts %>%
    group_by(subfamily) %>%
    slice_max(order_by = n_hits, n = top_n, with_ties = FALSE) %>%
    ungroup()
  
  best_stats <- fimo_sig %>%
    group_by(subfamily, rbp) %>%
    summarise(
      n_hits = n(),
      best_q = safe_min(`q-value`),
      best_p = safe_min(`p-value`),
      best_score = safe_max(score),
      .groups = "drop"
    )
  
  top_combined <- best_stats %>%
    arrange(subfamily, desc(n_hits), best_q, desc(best_score)) %>%
    group_by(subfamily) %>%
    slice_head(n = top_n) %>%
    ungroup()
  
  list(
    fimo_sig     = fimo_sig,
    counts       = counts,
    top_by_count = top_by_count,
    best_stats   = best_stats,
    top_combined = top_combined
  )
}

# ============================================================
# 3) Enrichment analysis (ALL RBPs by default)
# ============================================================
compute_enrichment_safe <- function(rbp_i, subfam_i, data) {
  a <- sum(data$rbp == rbp_i & data$subfamily == subfam_i, na.rm = TRUE)
  b <- sum(data$rbp == rbp_i & data$subfamily != subfam_i, na.rm = TRUE)
  c <- sum(data$rbp != rbp_i & data$subfamily == subfam_i, na.rm = TRUE)
  d <- sum(data$rbp != rbp_i & data$subfamily != subfam_i, na.rm = TRUE)
  
  if (any(!is.finite(c(a, b, c, d))) || (a + b + c + d == 0)) {
    return(tibble(
      rbp = rbp_i, subfamily = subfam_i, hits = a,
      odds_ratio = NA_real_, p_value = NA_real_
    ))
  }
  
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
  
  tibble(
    rbp = rbp_i,
    subfamily = subfam_i,
    hits = a,
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value
  )
}

run_enrichment_fast <- function(fimo_sig, rbps = NULL, subfamilies = NULL, full_grid = TRUE) {
  dat <- fimo_sig %>%
    dplyr::filter(!is.na(rbp), !is.na(subfamily))
  
  if (!is.null(rbps))       dat <- dat %>% dplyr::filter(rbp %in% rbps)
  if (!is.null(subfamilies)) dat <- dat %>% dplyr::filter(subfamily %in% subfamilies)
  
  # a = observed hits per (rbp, subfamily)
  a_tbl <- dat %>%
    dplyr::count(rbp, subfamily, name = "a")
  
  # margins
  rbp_tot <- dat %>% dplyr::count(rbp, name = "rbp_total")
  sub_tot <- dat %>% dplyr::count(subfamily, name = "subfam_total")
  N <- nrow(dat)
  
  # grid of pairs
  if (full_grid) {
    grid <- tidyr::expand_grid(rbp = rbp_tot$rbp, subfamily = sub_tot$subfamily)
  } else {
    # much faster + smaller result, but only pairs with at least one hit (a>0)
    grid <- a_tbl %>% dplyr::select(rbp, subfamily)
  }
  
  # join and compute b,c,d without rescanning dat
  tab <- grid %>%
    dplyr::left_join(a_tbl, by = c("rbp", "subfamily")) %>%
    dplyr::mutate(a = dplyr::coalesce(a, 0L)) %>%
    dplyr::left_join(rbp_tot, by = "rbp") %>%
    dplyr::left_join(sub_tot, by = "subfamily") %>%
    dplyr::mutate(
      b = rbp_total    - a,
      c = subfam_total - a,
      d = N - a - b - c
    )
  
  # Fisher test per row (still the heavy part, but now much cheaper per row)
  # Use mapply to avoid rowwise() overhead
  ft <- mapply(
    FUN = function(a, b, c, d) {
      mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
      out <- stats::fisher.test(mat)
      c(odds_ratio = unname(out$estimate), p_value = out$p.value)
    },
    a = tab$a, b = tab$b, c = tab$c, d = tab$d,
    SIMPLIFY = TRUE
  )
  
  tab$odds_ratio <- as.numeric(ft["odds_ratio", ])
  tab$p_value    <- as.numeric(ft["p_value", ])
  tab$p_adj      <- p.adjust(tab$p_value, method = "BH")
  
  tab %>%
    dplyr::select(rbp, subfamily, hits = a, odds_ratio, p_value, p_adj) %>%
    dplyr::arrange(subfamily, p_adj)
}

# ============================================================
# 4) Heatmap for HERVH vs HERVK
# ============================================================
plot_heatmap_hervh_hervk <- function(enrichment_results, subfamilies, out_file = NULL) {
  
  heatmap_df <- enrichment_results %>%
    filter(subfamily %in% subfamilies) %>%
    group_by(rbp) %>%
    filter(any(p_adj < 0.05, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(log2_or = log2(odds_ratio)) %>%
    select(rbp, subfamily, log2_or) %>%
    pivot_wider(names_from = subfamily, values_from = log2_or) %>%
    tibble::column_to_rownames("rbp") %>%
    as.data.frame()
  
  # Ensure column order
  if (all(subfamilies %in% colnames(heatmap_df))) {
    heatmap_df <- heatmap_df[, subfamilies, drop = FALSE]
  }
  
  # Cap extreme values
  cap <- 2
  heatmap_df_cap <- heatmap_df %>%
    mutate(across(everything(), ~ pmax(pmin(.x, cap), -cap))) %>%
    as.data.frame()
  
  # --- Open device if requested
  if (!is.null(out_file)) {
    ext <- tolower(tools::file_ext(out_file))
    
    if (ext == "pdf") {
      pdf(out_file, width = 6, height = 8)
      
    } else if (ext == "png") {
      png(out_file, width = 6, height = 8, units = "in", res = 300)
      
    } else if (ext == "svg") {
      svg(out_file, width = 6, height = 8)
      
    } else {
      stop("Unsupported format: use .pdf, .png or .svg")
    }
  }
  
  # --- Plot
  pheatmap::pheatmap(
    as.matrix(heatmap_df_cap),
    scale = "none",
    color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
    na_col = "grey90",
    cluster_rows = TRUE,   # ✅ removed row clustering
    cluster_cols = FALSE,
    fontsize = 11,           # ✅ smaller global font
    fontsize_row = 7,       # optional: even smaller row labels
    fontsize_col = 15
  )
  
  # --- Close device
  if (!is.null(out_file)) {
    dev.off()
  }
  
  invisible(heatmap_df_cap)
}

# ============================================================
# 5) Atlas: load + derive
# ============================================================
load_atlas <- function(path) {
  atlas_raw <- read_tsv(path, col_types = cols(.default = col_character()), progress = FALSE)
  stopifnot("HERV_id" %in% names(atlas_raw))
  
  atlas <- atlas_raw %>%
    mutate(
      subfamily = infer_subfamily(HERV_id),
      class = infer_class(subfamily),
      
      length_num       = as_num(length),
      domain_count_num = as_num(domain_count),
      rbp_burden_num   = as_num(rbp_burden),
      rbp_unique_num   = as_num(rbp_unique),
      ltr5_len_num     = as_num(ltr5_length),
      ltr3_len_num     = as_num(ltr3_length),
      
      mirna_pot_num    = as_int(mirna_potential),
      dsrna_pot_num    = as_int(dsrna_potential),
      dsrna_struct_num = as_int(dsrna_structural),
      dsrna_bidir_num  = as_int(dsrna_bidirectional),
      
      has_domain  = !is.na(domain_count_num) & domain_count_num > 0,
      has_ltr5    = !is.na(ltr5_len_num) & ltr5_len_num > 0,
      has_ltr3    = !is.na(ltr3_len_num) & ltr3_len_num > 0,
      has_ltr_any = has_ltr5 | has_ltr3,
      
      has_rbp   = !is.na(rbp_burden_num) & rbp_burden_num > 0,
      has_mirna = (mirna_pot_num == 1),
      
      has_dsrna_any    = (dsrna_pot_num == 1),
      has_dsrna_struct = (dsrna_struct_num == 1),
      has_dsrna_bidir  = (dsrna_bidir_num == 1)
    )
  
  list(raw = atlas_raw, derived = atlas)
}

atlas_subfamily_summary <- function(atlas) {
  atlas %>%
    group_by(subfamily) %>%
    summarise(
      n = n(),
      p_rbp   = prev(has_rbp),
      p_mirna = prev(has_mirna),
      p_dsrna_any    = prev(has_dsrna_any),
      p_dsrna_struct = prev(has_dsrna_struct),
      p_dsrna_bidir  = prev(has_dsrna_bidir),
      p_domain  = prev(has_domain),
      p_ltr_any = prev(has_ltr_any),
      p_intron     = prev(feature_overlap == "intron"),
      p_intergenic = prev(feature_overlap == "intergenic"),
      med_len        = median(length_num, na.rm = TRUE),
      med_rbp_burden = median(rbp_burden_num, na.rm = TRUE),
      med_rbp_unique = median(rbp_unique_num, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n))
}

# ============================================================
# 6) Plot: RBP burden vs length (labels: custom IDs OR exact top_n)
# ============================================================
plot_rbp_vs_length <- function(atlas, out_png = NULL, out_svg = NULL, label_ids = NULL, top_n = 10, log_x_axis = FALSE) {
  plot_df <- atlas %>%
    transmute(
      HERV_id, locid, subfamily,
      length = length_num,
      rbp_burden = rbp_burden_num
    ) %>%
    filter(!is.na(length), !is.na(rbp_burden), length > 0)
  
  if (!is.null(label_ids) && length(label_ids) > 0) {
    plot_df <- plot_df %>%
      mutate(is_label = HERV_id %in% label_ids)
    subtitle_extra <- paste0("Labeling provided IDs (n=", sum(plot_df$is_label), ").")
  } else {
    top_ids <- plot_df %>%
      arrange(desc(rbp_burden), desc(length), HERV_id) %>% # deterministic tiebreak
      slice_head(n = top_n) %>%
      pull(HERV_id)
    
    plot_df <- plot_df %>%
      mutate(is_label = HERV_id %in% top_ids)
    subtitle_extra <- paste0("Labeling top ", top_n, " loci by RBP burden.")
  }
  
  q <- quantile(plot_df$rbp_burden, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- as.numeric(q[2] - q[1])
  thr_extreme <- as.numeric(q[2] + 10 * iqr)
  
  
  if (log_x_axis) {
    p <- ggplot(plot_df, aes(x = length, y = rbp_burden)) +
      geom_point(alpha = 0.35, size = 2) +
      geom_point(data = subset(plot_df, is_label), size = 3) +
      ggrepel::geom_text_repel(
        data = subset(plot_df, is_label),
        aes(label = HERV_id),
        size = 5,
        max.overlaps = Inf
      ) +
      scale_x_log10() +
      labs(
        x = "Internal region length",
        y = "RBP burden",
      ) +
      theme_classic() +
      theme(
        axis.title = element_text(size = 16),   # axis titles (x/y labels)
        axis.text  = element_text(size = 14)    # axis tick numbers
      )
  } 
  else {
    p <- ggplot(plot_df, aes(x = length, y = rbp_burden)) +
      geom_point(alpha = 0.35, size = 2) +
      geom_point(data = subset(plot_df, is_label), size = 3) +
      ggrepel::geom_text_repel(
        data = subset(plot_df, is_label),
        aes(label = HERV_id),
        size = 5,
        max.overlaps = Inf
      ) +
      labs(
        x = "Internal region length",
        y = "RBP burden",
        title = "RBP burden vs locus length",
        subtitle = paste0(
          subtitle_extra,
          " Extreme threshold reference: rbp_burden ≥ ", round(thr_extreme, 2), " (Q3 + 10×IQR)."
        )
      ) +
      theme_classic() +
      theme(
        axis.title = element_text(size = 16),   # axis titles (x/y labels)
        axis.text  = element_text(size = 14)    # axis tick numbers
      )
  }
  
  
  if (!is.null(out_png)) ggsave(out_png, p, width = 12, height = 9, dpi = 300)
  if (!is.null(out_svg)) ggsave(out_svg, p, width = 12, height = 9, dpi = 300)
  p
}

# ============================================================
# 7) Convenience: pull FIMO rows for a given HERV_id (via locid mapping)
# ============================================================
get_fimo_for_herv_id <- function(fimo_parsed, atlas_derived, herv_id) {
  loc <- atlas_derived$locid[atlas_derived$HERV_id == herv_id]
  if (length(loc) == 0 || is.na(loc)) return(tibble())
  fimo_parsed %>% filter(sequence_name == loc)
}

# ============================================================
# Main
# ============================================================

# --- FIMO
fimo <- load_fimo(CFG$fimo_in)
fimo_summ <- summarise_fimo_by_subfamily(fimo$parsed, q_cut = CFG$q_cut, top_n = CFG$top_n)

write_tsv(fimo_summ$counts,       file.path(CFG$out_dir, "rbp_counts_by_subfamily.tsv"))
write_tsv(fimo_summ$top_by_count, file.path(CFG$out_dir, "top_rbps_by_subfamily_count.tsv"))
write_tsv(fimo_summ$best_stats,   file.path(CFG$out_dir, "rbp_beststats_by_subfamily.tsv"))
write_tsv(fimo_summ$top_combined, file.path(CFG$out_dir, "top_rbps_by_subfamily_combined.tsv"))

cat("\n=== FIMO ===\n")
cat("Loaded rows:", nrow(fimo$parsed), " | After q-cut filter:", nrow(fimo_summ$fimo_sig), "\n")

# --- Enrichment (ALL RBPs, ALL subfamilies)
enrichment_results <- run_enrichment_fast(fimo_summ$fimo_sig)
write_tsv(enrichment_results, file.path(CFG$out_dir, "rbp_enrichment_all.tsv"))

# --- Export HERVH + HERVK to one Excel (2 tabs)
out_xlsx <- file.path(CFG$out_dir, "rbp_enrichment_HERVH_HERVK.xlsx")

hervhk <- enrichment_results %>%
  filter(subfamily %in% c("HERVH", "HERVK")) %>%
  arrange(subfamily, p_adj, desc(hits))

df_hervh <- hervhk %>% filter(subfamily == "HERVH" & p_adj < 0.05)
df_hervk <- hervhk %>% filter(subfamily == "HERVK" & p_adj < 0.05)

wb <- createWorkbook()

addWorksheet(wb, "HERVH")
writeDataTable(wb, "HERVH", df_hervh)

addWorksheet(wb, "HERVK")
writeDataTable(wb, "HERVK", df_hervk)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

if (CFG$make_plots) {
  # Save as PDF
  plot_heatmap_hervh_hervk(
    enrichment_results,
    subfamilies = c("HERVH", "HERVK"),
    out_file = file.path(CFG$out_dir, "heatmap_HERVH_HERVK.png")
  )
  
  # Save as PNG
  plot_heatmap_hervh_hervk(
    enrichment_results,
    subfamilies = c("HERVH", "HERVK"),
    out_file = file.path(CFG$out_dir, "heatmap_HERVH_HERVK.png")
  )
}

if (CFG$make_plots) {
  plot_heatmap_hervh_hervk(enrichment_results, subfamilies = c("HERVH", "HERVK", "HERV17", "HERVK9"))
  # Save as PDF
  plot_heatmap_hervh_hervk(
    enrichment_results,
    subfamilies = c("HERVH", "HERVK", "HERV17", "HERVK9"),
    out_file = file.path(CFG$out_dir, "heatmap_HERVH_HERVK_HERV17_HERVK9.png")
  )
  
  # Save as PNG
  plot_heatmap_hervh_hervk(
    enrichment_results,
    subfamilies = c("HERVH", "HERVK", "HERV17", "HERVK9"),
    out_file = file.path(CFG$out_dir, "heatmap_HERVH_HERVK_HERV17_HERVK9.png")
  )
}

# --- Atlas
atlas <- load_atlas(CFG$atlas_in)
subfam_summary <- atlas_subfamily_summary(atlas$derived)
write_tsv(subfam_summary, file.path(CFG$out_dir, "atlas_subfamily_summary.tsv"))

if (CFG$make_plots) {
  p <- plot_rbp_vs_length(
    atlas$derived,
    top_n = CFG$top_n,
    out_png = file.path(CFG$out_dir, "rbp_burden_vs_length.png"),
    out_svg = file.path(CFG$out_dir, "rbp_burden_vs_length.svg")
  )
  print(p)
}

# --- Example: specific subfamily plot
atlas_hervip10fh <- atlas$derived %>% filter(subfamily == "HERVIP10FH")
if (CFG$make_plots) {
  p <- plot_rbp_vs_length(
    atlas_hervip10fh,
    top_n = 4,
    out_png = file.path(CFG$out_dir, "rbp_burden_vs_length_HERVIP10FH.png"),
    out_svg = file.path(CFG$out_dir, "rbp_burden_vs_length_HERVIP10FH.svg"),
    log_x_axis = TRUE
  )
  print(p)
}

# --- Example: inspect one locus
herv_id_focus <- "HERVIP10FH-int_00257"
fimo_focus <- get_fimo_for_herv_id(fimo$parsed, atlas$derived, herv_id_focus)
write_tsv(fimo_focus, file.path(CFG$out_dir, paste0("fimo_", herv_id_focus, ".tsv")))

cat("\nWrote outputs to:", CFG$out_dir, "\n")




### GENERAL STATISTICS ###

# % with >= 1 RBP motif
sum(atlas$derived$rbp_burden > 0) # 47076
sum(atlas$derived$rbp_burden > 0)/ nrow(as.data.frame(atlas$derived)) # 0.8852866

mean(as.numeric(atlas$derived$rbp_burden))


