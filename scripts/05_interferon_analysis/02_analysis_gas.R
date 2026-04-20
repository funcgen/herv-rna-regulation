# Interferon-associated LTR motif architecture from FIMO scans (STAT1 + IRF)
# ------------------------------------------------------------------------
# Goal:
#   1) Start from the full parsed FIMO table (all motifs, all LTR instances).
#   2) Restrict to LTR instances (sequence_name) that contain ≥1 STAT1 motif hit.
#      (STAT1 PWMs are used here as a proxy for GAS-like interferon-responsive sites.)
#   3) For those LTRs, quantify STAT1 motif burden per LTR using deduplicated sites:
#        - collapse palindromic/opposite-strand or highly overlapping hits (≥50% overlap)
#        - compute per-LTR counts, best p-value/score, motif density (sites/kb)
#   4) Compute local clustering of STAT1 sites within each LTR:
#        - max number of STAT1 sites in any 75 bp window
#        - coordinates of the densest 75 bp window and fraction of sites within it
#   5) Add IRF-family support (IRF1/2/3/7/9 etc.):
#        - collapse overlapping IRF hits (≥50% overlap) to define unique IRF sites
#        - compute per-LTR IRF counts and best p-value/score
#        - compute local IRF clustering within 75 bp windows
#   6) Compute combined STAT1+IRF “synergy” metrics:
#        - densest 75 bp window across both TF classes
#        - number of STAT1 and IRF sites in that best window
#   7) Compute a joint STAT1+IRF clustering table (stat_irf_cluster) and merge it back:
#        - Build a combined per-site table with STAT1 and IRF positions (LTR-relative midpoints).
#        - Scan 75 bp windows to find the densest "IFN window" across both TF classes.
#        - For each LTR, record:
#            * max_stat1_irf_in_75bp      : max total (STAT1+IRF) sites in any 75 bp window
#            * stat1_in_best_ifn_window   : number of STAT1 sites in that best window
#            * irf_in_best_ifn_window     : number of IRF sites in that best window
#            * best_ifn_win_start/end     : LTR-relative coordinates of the best 75 bp window
#        - Merge these metrics into stat_counts (by sequence_name) and set missing values to 0
#          for LTRs without IRF sites (or without any joint window signal).
#
#   8) Assign an interpretable categorical label (ifn_stimulation_potential) AFTER the merge:
#        - This label summarizes the local architecture of STAT1/IRF sites within 75 bp:
#            1_Very High : ≥1 STAT1 AND ≥1 IRF within the same window, and total ≥3 sites
#                          (dense synergistic module / potential enhancer-like cluster)
#            2_High      : ≥1 STAT1 AND ≥1 IRF within the same window, but total <3
#                          (synergistic pair; weaker but still cooperative architecture)
#            3_High      : strong single-factor cluster (≥3 STAT1 or ≥3 IRF in the best window)
#            4_Medium    : exactly a pair of the same factor in the best window
#            5_Low       : only one site in the best window
#            6_Unknown   : fallback label (should be rare; indicates unexpected edge case)
#        - IMPORTANT: ifn_stimulation_potential depends on the joint-window columns created
#          in stat_irf_cluster, so it must be computed only after stat_irf_cluster is merged.
#
# Notes:
#   - q-values are not used because scanning was parallelized by subfamily, so FDR is
#     not globally comparable. We rely on the fixed p-value cutoff used in the FIMO runs.
#   - Motif strand is not used for filtering/deduplication; TFs bind dsDNA. Strand is
#     only relevant later when mapping LTRs to loci and interpreting transcription direction.
#   - Window coordinates (e.g., densest_*_win_start/end) refer to the 75 bp window in
#     LTR-relative coordinates, not the motif span itself.
#
# Inputs:
#   - merged_fimo_parsed_v4.tsv: full set of motif hits (all motifs, all LTR instances)
#   - fimo_parsed_v4_STAT1.tsv: helper file used to define the set of LTRs with ≥1 STAT1 hit
#
# Outputs:
#   - fimo_parsed_v4_STAT1_only.tsv: subset of the full table restricted to LTRs with ≥1 STAT1 hit
#   - LTR_IFN_STAT1_IRF_summary.tsv (stat_counts_out): one row per LTR with STAT1/IRF burden,
#     local clustering metrics, joint STAT1+IRF window metrics, and ifn_stimulation_potential


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
fimo_stat_file <- "../../../data/fimo_parsed_v4_STAT1.tsv"
out_dir <- "../../../results/interferon/"
out_file       <- "../../../results/interferon/fimo_parsed_v4_STAT1_only.tsv"
stat_counts_out <- "../../../results/interferon/LTR_IFN_STAT1_summary.tsv"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------
# Load tables
# -----------------------
fimo_all  <- fread(fimo_all_file)
fimo_stat <- fread(fimo_stat_file)

# Ensure fimo_stat has the same columns as fimo_all (common when fimo_stat is a subset)
# If fimo_stat only has a subset of columns but contains sequence_name, that's OK too.
if (ncol(fimo_stat) == ncol(fimo_all)) {
  setnames(fimo_stat, names(fimo_all))
}

# -----------------------
# Restrict to LTRs that have ≥1 STAT1 site
# -----------------------
if (!"sequence_name" %in% names(fimo_all)) stop("fimo_all: missing sequence_name")
if (!"sequence_name" %in% names(fimo_stat)) stop("fimo_stat: missing sequence_name")

stat_ltr_set <- unique(fimo_stat$sequence_name)

# Keep *all motif hits* for those LTRs (useful later to compute IRF/NFKB/AP1 support etc.)
filtered <- fimo_all[sequence_name %in% stat_ltr_set]

# Export the subset (optional but convenient for reproducibility / HPC re-use)
fwrite(filtered, out_file, sep = "\t")

# Sanity checks
cat("FIMO all rows:              ", nrow(fimo_all), "\n")
cat("LTRs with STAT1 (≥1):", length(stat_ltr_set), "\n")
cat("Filtered rows (all motifs): ", nrow(filtered), "\n\n")

# -----------------------
# Summarize STAT1 burden per LTR (deduplicated)
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

# Raw STAT1 hits
stat_hits <- filtered[motif_alt_id == "STAT1"]

# Deduplicate palindromic / overlapping calls within each LTR
stat_hits_u <- stat_hits[, collapse_overlaps_minfrac_generic(.SD, min_frac = 0.5), by = sequence_name]

# Position (midpoint) for clustering
stat_hits_u[, pos := as.integer(floor((start + stop) / 2))]
setorder(stat_hits_u, sequence_name, pos)

# Per-LTR summary (IMPORTANT: use stat_hits_u, not stat_hits)
stat_counts <- stat_hits_u[, .(
  n_stat1          = .N,
  best_p_stat1     = min(`p-value`),
  best_score_stat1 = max(score),
  chrom            = first(chrom),
  ltr_start        = first(ltr_start),
  ltr_end          = first(ltr_end),
  ltr_strand       = first(ltr_strand),
  subfamily        = first(subfamily)
), by = sequence_name]

# LTR length + motif density
stat_counts[, ltr_len := ltr_end - ltr_start + 1L]
stat_counts[, stat_sites_per_kb := n_stat1 / (ltr_len / 1000)]
stat_counts[, log_ltr_len := log10(ltr_len)]

# Local clustering (on deduplicated hits)
window_bp <- 75L

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
    max_stat1_in_75bp = max_in_win,
    densest_win_start = max_start,
    densest_win_end   = max_end,
    prop_in_densest_75bp = max_in_win / n
  )
}, by = sequence_name]

# Merge clustering into stat_counts
stat_counts <- merge(stat_counts, cluster_stats, by = "sequence_name", all.x = TRUE)





# -----------------------
# IRF support (cooperative IFN program evidence)
# -----------------------
# -----------------------
# IRF support (cooperative IFN program evidence)
# -----------------------
irf_hits <- filtered[grepl("^IRF", motif_alt_id)]

# Collapse overlapping IRF hits per LTR (unique sites)
irf_sites <- irf_hits[
  , collapse_overlaps_minfrac_generic(.SD, min_frac = 0.5),
  by = sequence_name
]

# Per-LTR IRF summary (unique-site based)
irf_flag <- irf_sites[, .(
  has_irf        = .N > 0,
  n_irf          = .N,                 # unique IRF sites
  best_p_irf     = min(`p-value`),
  best_score_irf = max(score)
), by = sequence_name]

# Merge + fill
stat_counts <- merge(stat_counts, irf_flag, by = "sequence_name", all.x = TRUE)
stat_counts[is.na(has_irf), has_irf := FALSE]
stat_counts[is.na(n_irf), n_irf := 0L]
# (best_p_irf / best_score_irf can stay NA if no IRFs)

# -----------------------
# IRF clustering in 75 bp windows (using unique sites)
# -----------------------
window_bp <- 75L

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

# Merge + fill
stat_counts <- merge(stat_counts, irf_cluster, by = "sequence_name", all.x = TRUE)
stat_counts[is.na(max_irf_in_75bp), `:=`(
  max_irf_in_75bp = 0L,
  prop_irf_in_densest_75bp = 0
)]






## STAT1 + IRF clusters
# Build a joint table
# Use deduplicated STAT1 sites
stat_sites_pos <- stat_hits_u[, .(
  sequence_name,
  pos = as.integer(floor((start + stop) / 2)),
  tf = "STAT1"
)]

# Use deduplicated IRF sites (irf_sites) returned by the generic collapse function
irf_sites_pos <- irf_sites[, .(
  sequence_name,
  pos = as.integer(floor((start + stop) / 2)),
  tf = "IRF"
)]

# Combine
stat_irf_sites <- rbind(stat_sites_pos, irf_sites_pos)
setorder(stat_irf_sites, sequence_name, pos)

# Cluster in 75 bp window
window_bp <- 75L


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
      best_stat <- sum(tt[idx] == "STAT1")
      best_irf  <- sum(tt[idx] == "IRF")
      best_start <- p[i]
      best_end   <- p[i] + window_bp
    }
  }
  
  .(
    max_stat1_irf_in_75bp    = max_total,
    stat1_in_best_ifn_window = best_stat,
    irf_in_best_ifn_window   = best_irf,
    best_ifn_win_start       = best_start,
    best_ifn_win_end         = best_end
  )
  
}, by = sequence_name]


# Merge joint-window metrics into stat_counts
stat_counts <- merge(stat_counts, stat_irf_cluster, by = "sequence_name", all.x = TRUE)

# Fill missing (LTRs with no IRF sites will still have STAT1-only windows; but safe to fill)
stat_counts[is.na(max_stat1_irf_in_75bp), `:=`(
  max_stat1_irf_in_75bp = 0L,
  stat1_in_best_ifn_window = 0L,
  irf_in_best_ifn_window = 0L
)]


# Apply the classification logic
stat_counts[, ifn_stimulation_potential := fcase(
  # 1. Very High: >=1 STAT1, >=1 IRF, and total >= 3
  stat1_in_best_ifn_window >= 1 & irf_in_best_ifn_window >= 1 & max_stat1_irf_in_75bp >= 3, '1_Very High (Synergistic Cluster)',
  
  # 2. High: Exactly 1 STAT1 and 1 IRF (total will be 2)
  stat1_in_best_ifn_window >= 1 & irf_in_best_ifn_window >= 1 & max_stat1_irf_in_75bp < 3,  '2_High (Synergistic Pair)',
  
  # 3. High Clusters (Pure STAT1 or Pure IRF)
  stat1_in_best_ifn_window >= 3, '3_High (STAT1 Cluster)',
  irf_in_best_ifn_window >= 3,   '3_High (IRF Cluster)',    
  
  # 4. Medium: Exactly 2 of the same motif
  max_stat1_irf_in_75bp == 2 & stat1_in_best_ifn_window == 2, '4_Medium (STAT1 Pair)',
  max_stat1_irf_in_75bp == 2 & irf_in_best_ifn_window == 2,   '4_Medium (IRF Pair)',
  
  # 5. Low: Only 1 site in the 75bp window
  max_stat1_irf_in_75bp == 1, '5_Low (Isolated Site)',
  
  # 6. Fallback
  default = '6_Unknown'
)]

# Check the distribution again!
stat_counts[, .N, by = ifn_stimulation_potential][order(ifn_stimulation_potential)]

# Check the distribution of your new column
stat_counts[, .N, by = ifn_stimulation_potential][order(ifn_stimulation_potential)]



# -----------------------
# Export per-LTR IFN summary table
# -----------------------
# This table is the main reusable output of the IFN-LTR analysis:
# one row per LTR, with motif burden, density, IRF support, and IFN class.
fwrite(
  stat_counts,
  stat_counts_out,
  sep = "\t"
)

cat("Exported stat_counts to:", stat_counts_out, "\n")
