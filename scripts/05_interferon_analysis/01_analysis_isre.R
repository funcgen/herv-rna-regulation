# Interferon-associated LTR motif architecture from FIMO scans (STAT1::STAT2 + IRF)
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

# -----------------------
# Paths / configuration
# -----------------------
fimo_all_file  <- "../../../data/merged_fimo_parsed_v4.tsv"
fimo_stat_file <- "../../../data/fimo_parsed_v4_STAT1_STAT2.tsv"
out_dir <- "../../../results/interferon/"
out_file       <- "../../../results/interferon/fimo_parsed_v4_STAT1_STAT2_only.tsv"
stat_counts_out <- "../../../results/interferon/LTR_IFN_STAT1STAT2_IRF_summary.tsv"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Define the STAT-layer motif_alt_id to use everywhere
STAT_TF_ALTID <- "STAT1::STAT2"  # <- change if your motif_alt_id differs

# -----------------------
# Load tables
# -----------------------
fimo_all  <- fread(fimo_all_file)
fimo_stat <- fread(fimo_stat_file)

if (ncol(fimo_stat) == ncol(fimo_all)) {
  setnames(fimo_stat, names(fimo_all))
}

# -----------------------
# ==============================================================================
# 01_LTR_interferon_motif_architecture.R
# ==============================================================================

# Purpose
# Identify and classify interferon-responsive LTRs based on local motif architecture.
#
# This script integrates FIMO motif predictions to quantify the spatial
# organization of STAT (STAT1 or STAT1::STAT2) and IRF binding sites within LTRs,
# and derives an interpretable interferon stimulation potential score.
#
# Biological rationale
# Interferon-responsive enhancers are characterized not only by the presence of
# STAT and IRF motifs, but by their co-localization within short genomic windows.
# Clusters of these motifs enable cooperative binding and amplify transcriptional response.
#
# Key idea
# LTRs are treated as candidate regulatory modules whose IFN potential depends on:
#   - motif abundance (STAT / IRF burden)
#   - local density (clustering within 75 bp windows)
#   - combinatorial architecture (STAT + IRF co-occurrence)
#
# Pipeline overview
#   1. Load full FIMO motif scan (all motifs, all LTRs)
#   2. Restrict to LTRs containing ≥1 STAT motif
#   3. Collapse overlapping motif hits into unique binding sites
#   4. Compute per-LTR motif burden and density
#   5. Identify densest 75 bp windows for:
#        - STAT sites
#        - IRF sites
#        - combined STAT + IRF sites
#   6. Quantify cooperative STAT–IRF architecture
#   7. Assign categorical IFN stimulation potential per LTR
#
# IFN stimulation potential classes
#   1_Very High  : dense STAT + IRF cluster (≥3 sites in same window)
#   2_High       : STAT + IRF pair (co-localized but sparse)
#   3_High       : strong single-factor cluster (STAT or IRF ≥3)
#   4_Medium     : motif pair (same factor)
#   5_Low        : isolated site
#   6_Unknown    : fallback / edge cases
#
# Inputs
#   - merged_fimo_parsed_v4.tsv         (full motif scan)
#   - fimo_parsed_v4_STAT*.tsv          (STAT-containing LTRs)
#
# Outputs
#   - Filtered FIMO table (STAT-positive LTRs)
#   - Per-LTR summary table with motif architecture and IFN classification
#
# Notes
#   - Overlapping motif hits are collapsed (≥50% overlap) to define unique sites
#   - Motif strand is ignored (TF binding is strand-independent)
#   - Fixed p-value thresholds are used (FDR not comparable across subfamily scans)
#   - All coordinates are LTR-relative
#
# Position in pipeline
#   This is the first script of the "LTR regulatory architecture" module.
#   It produces the core IFN annotation layer used in downstream analyses.
# ==============================================================================


# Restrict to LTRs that have ≥1 STAT1::STAT2 site
# -----------------------
if (!"sequence_name" %in% names(fimo_all))  stop("fimo_all: missing sequence_name")
if (!"sequence_name" %in% names(fimo_stat)) stop("fimo_stat: missing sequence_name")

stat_ltr_set <- unique(fimo_stat$sequence_name)

filtered <- fimo_all[sequence_name %in% stat_ltr_set]
fwrite(filtered, out_file, sep = "\t")

cat("FIMO all rows:              ", nrow(fimo_all), "\n")
cat("LTRs with STAT1::STAT2 (≥1):", length(stat_ltr_set), "\n")
cat("Filtered rows (all motifs): ", nrow(filtered), "\n\n")

# -----------------------
# Helpers
# -----------------------
collapse_overlaps_minfrac_generic <- function(dt, min_frac = 0.5) {
  if (nrow(dt) == 0) return(dt)
  
  d <- copy(dt)
  d[, `:=`(start = as.integer(start), stop = as.integer(stop))]
  d[, width := stop - start + 1L]
  setorder(d, start, stop)
  
  groups <- integer(nrow(d))
  g <- 1L
  groups[1] <- g
  
  cur_start <- d$start[1]
  cur_stop  <- d$stop[1]
  cur_width <- d$width[1]
  
  if (nrow(d) >= 2) {
    for (i in 2:nrow(d)) {
      s <- d$start[i]; e <- d$stop[i]; w <- d$width[i]
      inter <- max(0L, min(cur_stop, e) - max(cur_start, s) + 1L)
      overlap_frac <- if (inter > 0L) inter / min(cur_width, w) else 0
      
      if (overlap_frac >= min_frac) {
        groups[i] <- g
        cur_start <- min(cur_start, s)
        cur_stop  <- max(cur_stop,  e)
        cur_width <- cur_stop - cur_start + 1L
      } else {
        g <- g + 1L
        groups[i] <- g
        cur_start <- s
        cur_stop  <- e
        cur_width <- w
      }
    }
  }
  
  d[, site_id := groups]
  
  # pick ONE representative hit per site (best p-value; break ties by score)
  d[order(`p-value`, -score), .SD[1], by = site_id][, site_id := NULL][]
}

window_bp <- 75L

# -----------------------
# Summarize STAT1::STAT2 burden per LTR (deduplicated)
# -----------------------
stat_hits <- filtered[motif_alt_id == STAT_TF_ALTID]
# If needed instead of exact match:
# stat_hits <- filtered[grepl("^STAT1::STAT2", motif_alt_id)]

stat_hits_u <- stat_hits[, collapse_overlaps_minfrac_generic(.SD, min_frac = 0.5), by = sequence_name]

stat_hits_u[, pos := as.integer(floor((start + stop) / 2))]
setorder(stat_hits_u, sequence_name, pos)

stat_counts <- stat_hits_u[, .(
  n_stat1stat2          = .N,
  best_p_stat1stat2     = min(`p-value`),
  best_score_stat1stat2 = max(score),
  chrom            = first(chrom),
  ltr_start        = first(ltr_start),
  ltr_end          = first(ltr_end),
  ltr_strand       = first(ltr_strand),
  subfamily        = first(subfamily)
), by = sequence_name]

stat_counts[, ltr_len := ltr_end - ltr_start + 1L]
stat_counts[, stat1stat2_sites_per_kb := n_stat1stat2 / (ltr_len / 1000)]
stat_counts[, log_ltr_len := log10(ltr_len)]

cluster_stats <- stat_hits_u[, {
  p <- pos
  n <- length(p)
  j <- 1L
  max_in_win <- 0L
  max_start <- NA_integer_
  max_end   <- NA_integer_
  
  for (i in seq_len(n)) {
    if (j < i) j <- i
    while (j <= n && p[j] <= p[i] + window_bp) j <- j + 1L
    k <- j - i
    if (k > max_in_win) {
      max_in_win <- k
      max_start <- p[i]
      max_end   <- p[i] + window_bp
    }
  }
  
  .(
    max_stat1stat2_in_75bp = max_in_win,
    densest_stat1stat2_win_start = max_start,
    densest_stat1stat2_win_end   = max_end,
    prop_stat1stat2_in_densest_75bp = max_in_win / n
  )
}, by = sequence_name]

stat_counts <- merge(stat_counts, cluster_stats, by = "sequence_name", all.x = TRUE)

# -----------------------
# IRF support
# -----------------------
irf_hits <- filtered[grepl("^IRF", motif_alt_id)]

irf_sites <- irf_hits[
  , collapse_overlaps_minfrac_generic(.SD, min_frac = 0.5),
  by = sequence_name
]

irf_flag <- irf_sites[, .(
  has_irf        = .N > 0,
  n_irf          = .N,
  best_p_irf     = min(`p-value`),
  best_score_irf = max(score)
), by = sequence_name]

stat_counts <- merge(stat_counts, irf_flag, by = "sequence_name", all.x = TRUE)
stat_counts[is.na(has_irf), has_irf := FALSE]
stat_counts[is.na(n_irf), n_irf := 0L]

irf_sites[, irf_pos := as.integer(floor((start + stop) / 2))]
setorder(irf_sites, sequence_name, irf_pos)

irf_cluster <- irf_sites[, {
  p <- irf_pos
  n <- length(p)
  j <- 1L
  
  max_in_win <- 0L
  max_start <- NA_integer_
  max_end   <- NA_integer_
  
  for (i in seq_len(n)) {
    if (j < i) j <- i
    while (j <= n && p[j] <= p[i] + window_bp) j <- j + 1L
    k <- j - i
    if (k > max_in_win) {
      max_in_win <- k
      max_start <- p[i]
      max_end   <- p[i] + window_bp
    }
  }
  
  .(
    max_irf_in_75bp = max_in_win,
    densest_irf_win_start = max_start,
    densest_irf_win_end   = max_end,
    prop_irf_in_densest_75bp = if (n > 0) max_in_win / n else 0
  )
}, by = sequence_name]

stat_counts <- merge(stat_counts, irf_cluster, by = "sequence_name", all.x = TRUE)
stat_counts[is.na(max_irf_in_75bp), `:=`(
  max_irf_in_75bp = 0L,
  prop_irf_in_densest_75bp = 0
)]

# -----------------------
# STAT1::STAT2 + IRF joint clusters (75 bp window)
# -----------------------
stat_sites_pos <- stat_hits_u[, .(
  sequence_name,
  pos = as.integer(floor((start + stop) / 2)),
  tf = "STAT1::STAT2"
)]

irf_sites_pos <- irf_sites[, .(
  sequence_name,
  pos = as.integer(floor((start + stop) / 2)),
  tf = "IRF"
)]

stat_irf_sites <- rbind(stat_sites_pos, irf_sites_pos)
setorder(stat_irf_sites, sequence_name, pos)

stat_irf_cluster <- stat_irf_sites[, {
  
  p  <- pos
  tt <- tf
  n <- length(p)
  j <- 1L
  
  max_total <- 0L
  best_stat <- 0L
  best_irf  <- 0L
  best_start <- NA_integer_
  best_end   <- NA_integer_
  
  for (i in seq_len(n)) {
    if (j < i) j <- i
    while (j <= n && p[j] <= p[i] + window_bp) j <- j + 1L
    
    idx <- i:(j - 1L)
    k <- length(idx)
    
    if (k > max_total) {
      max_total <- k
      best_stat <- sum(tt[idx] == "STAT1::STAT2")
      best_irf  <- sum(tt[idx] == "IRF")
      best_start <- p[i]
      best_end   <- p[i] + window_bp
    }
  }
  
  .(
    max_stat1stat2_irf_in_75bp    = max_total,
    stat1stat2_in_best_ifn_window = best_stat,
    irf_in_best_ifn_window        = best_irf,
    best_ifn_win_start            = best_start,
    best_ifn_win_end              = best_end
  )
  
}, by = sequence_name]

stat_counts <- merge(stat_counts, stat_irf_cluster, by = "sequence_name", all.x = TRUE)

# In principle these should exist for all LTRs (they all have STAT1::STAT2 by construction),
# but keep safe defaults anyway.
stat_counts[is.na(max_stat1stat2_irf_in_75bp), `:=`(
  max_stat1stat2_irf_in_75bp = 0L,
  stat1stat2_in_best_ifn_window = 0L,
  irf_in_best_ifn_window = 0L
)]

# -----------------------
# Categorical label (updated for STAT1::STAT2)
# -----------------------
stat_counts[, ifn_stimulation_potential := fcase(
  stat1stat2_in_best_ifn_window >= 1 & irf_in_best_ifn_window >= 1 & max_stat1stat2_irf_in_75bp >= 3,
  "1_Very High (Synergistic Cluster)",
  
  stat1stat2_in_best_ifn_window >= 1 & irf_in_best_ifn_window >= 1 & max_stat1stat2_irf_in_75bp < 3,
  "2_High (Synergistic Pair)",
  
  stat1stat2_in_best_ifn_window >= 3, "3_High (STAT1::STAT2 Cluster)",
  irf_in_best_ifn_window >= 3,       "3_High (IRF Cluster)",
  
  max_stat1stat2_irf_in_75bp == 2 & stat1stat2_in_best_ifn_window == 2, "4_Medium (STAT1::STAT2 Pair)",
  max_stat1stat2_irf_in_75bp == 2 & irf_in_best_ifn_window == 2,       "4_Medium (IRF Pair)",
  
  max_stat1stat2_irf_in_75bp == 1, "5_Low (Isolated Site)",
  
  default = "6_Unknown"
)]

stat_counts[, .N, by = ifn_stimulation_potential][order(ifn_stimulation_potential)]

# -----------------------
# Export per-LTR IFN summary table
# -----------------------
fwrite(stat_counts, stat_counts_out, sep = "\t")
cat("Exported stat_counts to:", stat_counts_out, "\n")
