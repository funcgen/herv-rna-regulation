#!/usr/bin/env Rscript
# ==============================================================================
# 07_sparcs step 2 (frozen global order 43, effective order 44), LOCAL, offline:
# per-LTR interferon STAT / IRF motif architecture from the published LTR FIMO scan (SRC09).
#
# REFACTOR_SCRIPT of D_interferon_and_motif_architectur/01_analysis_isre.R (sha256 4a5e8311...,
# STAT1::STAT2 + IRF, "ISGF3" / alpha-beta proxy) and 02_analysis_gas.R (sha256 ca64bcce..., STAT1 + IRF,
# "GAF" / gamma proxy), unified into one parameterised step. The two historical scripts differ only in
# the STAT motif and in output column names / labels; the scientific code below is theirs, verbatim:
#
#   source       SRC09 merged_fimo_parsed_v4.tsv (JASPAR CORE matrices scanned by FIMO over every SRC07
#                LTR; motif start/stop 1-based inclusive, LTR-relative; ltr_start = BED start + 1). No
#                p-value, q-value or score filter here: the upstream scan's fixed threshold is used as is.
#   LTR set      LTRs (sequence_name) with >= 1 row whose motif_alt_id EQUALS the profile's STAT motif
#                (historically read from a helper file byte-identical to exactly that SRC09 subset)
#   STAT hits    motif_alt_id == "STAT1::STAT2" (ISGF3) or == "STAT1" (GAF): exact, case-sensitive
#   IRF hits     grepl("^IRF", motif_alt_id, ignore.case = TRUE): IRF1-IRF9. D040: the historical rule was
#                case-sensitive and dropped the JASPAR IRF1 matrix, named "Irf1" (see PROVENANCE.md)
#   collapse     per LTR and per class (the STAT motif; all IRF motifs pooled), strand ignored: hits sorted
#                by (start, stop); a hit joins the running group when
#                  inter / min(width(group span), width(hit)) >= sparcs.ifn_motif_overlap_collapse_fraction
#                with inclusive widths (stop - start + 1) and inter = min(stops) - max(starts) + 1; the
#                group span then grows to the union, so chaining is transitive through the span.
#                One representative per group: lowest p-value, then highest score, then first in
#                (start, stop) order.
#   window       site position = floor((start + stop) / 2) of the representative; for each site i the
#                window holds every site j with p_i <= p_j <= p_i + sparcs.ifn_cluster_window_bp
#                (inclusive); the densest window is the first one reaching the maximum count. Computed
#                for STAT sites, IRF sites and STAT + IRF sites together (the joint window gives the
#                STAT and IRF counts used by ifn_stimulation_potential). No cluster table exists.
#   output unit  one row per LTR of the profile's LTR set (merge by sequence_name, so C-locale order)
#
# Release changes (declared before production):
#   * the STAT LTR set is derived from SRC09 directly instead of the two unregistered helper files
#     (proved byte-identical to the exact-match subsets); the 1.0-1.4 GB filtered copies of SRC09
#     (fimo_parsed_v4_STAT1*_only.tsv) had no consumer and are not written
#   * 0.50 and 75 come from config/parameters.yaml (sparcs:), never literals
#   * SRC09 is content-verified (md5 of the registered member), schema, coordinate and identity guards
#   * a site table (every collapsed site with its group size and span) and a QC table are added;
#     the two per-LTR summaries keep their historical names, columns and labels
#   * D040 (2026-10-09): IRF motifs are matched case-insensitively, so "Irf1" is included; the per-LTR
#     summaries differ from the historical files in the IRF-dependent columns only
#   * staged, atomic output with checksums.sha256
# Usage: summarise_ifn_motifs.R --fimo <SRC09 member> --archive-members 01_sources/manifests/archive_members.tsv
#                               --params config/parameters.yaml --outdir <dir>
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("missing argument ", flag)
  args[[i + 1]]
}
FIMO <- get_arg("--fimo")
MEMBERS <- get_arg("--archive-members")
PARAMS <- get_arg("--params")
FINAL_OUT <- get_arg("--outdir")

HEADER <- c("subfamily", "motif_id", "motif_alt_id", "sequence_name", "start", "stop", "strand", "score",
            "p-value", "q-value", "matched_sequence", "chrom", "ltr_start", "ltr_end", "ltr_strand")
# the two historical profiles: STAT motif, output file and the GAF column / label names
PROFILES <- list(
  ISGF3 = list(stat = "STAT1::STAT2", file = "LTR_IFN_STAT1STAT2_IRF_summary.tsv", rename = NULL),
  GAF   = list(stat = "STAT1", file = "LTR_IFN_STAT1_summary.tsv", rename = c(
    n_stat1stat2 = "n_stat1", best_p_stat1stat2 = "best_p_stat1", best_score_stat1stat2 = "best_score_stat1",
    stat1stat2_sites_per_kb = "stat_sites_per_kb", max_stat1stat2_in_75bp = "max_stat1_in_75bp",
    densest_stat1stat2_win_start = "densest_win_start", densest_stat1stat2_win_end = "densest_win_end",
    prop_stat1stat2_in_densest_75bp = "prop_in_densest_75bp", max_stat1stat2_irf_in_75bp = "max_stat1_irf_in_75bp",
    stat1stat2_in_best_ifn_window = "stat1_in_best_ifn_window"))
)
IRF_PATTERN <- "^IRF"

# ------------------------------------------------------------------------------
# RELEASE: parameters and input guards
# ------------------------------------------------------------------------------
read_sparcs_param <- function(path, key) {
  txt <- readLines(path, warn = FALSE)
  s <- grep("^sparcs:", txt)
  if (length(s) != 1) stop("no sparcs: block in ", path)
  blk <- txt[(s + 1):length(txt)]
  end <- grep("^\\S", blk)
  if (length(end)) blk <- blk[seq_len(end[1] - 1)]
  hit <- grep(paste0("^\\s+", key, ":\\s*[0-9.eE+-]+\\s*(#.*)?$"), blk, value = TRUE)
  if (length(hit) != 1) stop("parameter sparcs.", key, " not found")
  as.numeric(sub(paste0("^\\s+", key, ":\\s*([0-9.eE+-]+).*$"), "\\1", hit))
}
MIN_FRAC <- read_sparcs_param(PARAMS, "ifn_motif_overlap_collapse_fraction")
WINDOW_BP <- read_sparcs_param(PARAMS, "ifn_cluster_window_bp")
if (!(is.finite(MIN_FRAC) && MIN_FRAC > 0 && MIN_FRAC <= 1)) stop("ifn_motif_overlap_collapse_fraction must be in (0, 1]")
if (!(is.finite(WINDOW_BP) && WINDOW_BP >= 1 && WINDOW_BP == round(WINDOW_BP))) stop("ifn_cluster_window_bp must be a positive integer")
WINDOW_BP <- as.integer(WINDOW_BP)
message("[ok] parameters: overlap collapse fraction ", MIN_FRAC, "; window ", WINDOW_BP, " bp")

m <- fread(MEMBERS, sep = "\t", colClasses = "character")
want_md5 <- m[get(names(m)[1]) == "SRC09"][[5]]
if (length(want_md5) != 1) stop("SRC09 member md5 not found in ", MEMBERS)
got_md5 <- unname(tools::md5sum(FIMO))
if (!identical(got_md5, want_md5)) stop(FIMO, ": md5 ", got_md5, " is not the registered SRC09 member md5 ", want_md5)
message("[ok] SRC09 member md5 ", got_md5, " (", MEMBERS, ")")
if (!identical(names(fread(FIMO, nrows = 0)), HEADER)) stop(FIMO, ": not the SRC09 15-column FIMO schema")

# matched_sequence and q-value are never used; everything else is read with fread's own type detection
fimo_all <- fread(FIMO, select = setdiff(HEADER, c("matched_sequence", "q-value")))
n_src_rows <- nrow(fimo_all)
n_src_ltrs <- uniqueN(fimo_all$sequence_name)
alt_counts <- fimo_all[, .N, by = motif_alt_id]
stat_ids <- vapply(PROFILES, `[[`, "", "stat")
ifn <- fimo_all[motif_alt_id %in% stat_ids | grepl(IRF_PATTERN, motif_alt_id, ignore.case = TRUE)]   # SRC09 row order kept
rm(fimo_all); invisible(gc())
message("[ok] SRC09: ", n_src_rows, " rows, ", n_src_ltrs, " LTRs; STAT/IRF rows ", nrow(ifn))

bad <- ifn[!(strand %in% c("+", "-")) | !(ltr_strand %in% c("+", "-")) | is.na(start) | is.na(stop) |
             start < 1L | stop < start | ltr_end < ltr_start | stop > ltr_end - ltr_start + 1L]
if (nrow(bad)) stop(nrow(bad), " STAT/IRF row(s) with an invalid strand or a motif outside its LTR; first: ", bad$sequence_name[1])
ltr_id <- unique(ifn[, .(sequence_name, subfamily, chrom, ltr_start, ltr_end, ltr_strand)])
if (anyDuplicated(ltr_id$sequence_name)) stop("an LTR sequence_name maps to more than one coordinate record")
enc <- ltr_id[sequence_name != paste0(subfamily, "_merged_pos_", chrom, "_", ltr_start, "_", ltr_end, "_strand_", ltr_strand)]
if (nrow(enc)) stop(nrow(enc), " LTR(s) whose sequence_name does not encode its subfamily/chrom/start/end/strand; first: ", enc$sequence_name[1])
n_exact_dup <- sum(duplicated(ifn))   # reported, never removed: identical calls collapse by the historical rule

# ------------------------------------------------------------------------------
# Historical helper (verbatim, apart from the RELEASE line recording group size and span)
# ------------------------------------------------------------------------------
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
  d[, `:=`(site_n_hits = .N, site_start = min(start), site_stop = max(stop)), by = site_id]   # RELEASE: site table only

  # pick ONE representative hit per site (best p-value; break ties by score)
  d[order(`p-value`, -score), .SD[1], by = site_id][, site_id := NULL][]
}

# densest-window scan of the historical scripts (identical loop in all three places)
densest_window <- function(p, window_bp) {
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
  list(max_in_win = max_in_win, max_start = max_start, max_end = max_end, n = n)
}

# ------------------------------------------------------------------------------
# One profile: the historical 01_analysis_isre.R pipeline (ISRE names; GAF renamed afterwards)
# ------------------------------------------------------------------------------
summarise_profile <- function(ifn, STAT_TF_ALTID, window_bp, min_frac) {
  stat_ltr_set <- unique(ifn[motif_alt_id == STAT_TF_ALTID]$sequence_name)
  filtered <- ifn[sequence_name %in% stat_ltr_set]

  stat_hits <- filtered[motif_alt_id == STAT_TF_ALTID]
  stat_hits_u <- stat_hits[, collapse_overlaps_minfrac_generic(.SD, min_frac = min_frac), by = sequence_name]
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
    w <- densest_window(pos, window_bp)
    .(
      max_stat1stat2_in_75bp = w$max_in_win,
      densest_stat1stat2_win_start = w$max_start,
      densest_stat1stat2_win_end   = w$max_end,
      prop_stat1stat2_in_densest_75bp = w$max_in_win / w$n
    )
  }, by = sequence_name]
  stat_counts <- merge(stat_counts, cluster_stats, by = "sequence_name", all.x = TRUE)

  irf_hits <- filtered[grepl(IRF_PATTERN, motif_alt_id, ignore.case = TRUE)]
  irf_sites <- irf_hits[, collapse_overlaps_minfrac_generic(.SD, min_frac = min_frac), by = sequence_name]
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
    w <- densest_window(irf_pos, window_bp)
    .(
      max_irf_in_75bp = w$max_in_win,
      densest_irf_win_start = w$max_start,
      densest_irf_win_end   = w$max_end,
      prop_irf_in_densest_75bp = if (w$n > 0) w$max_in_win / w$n else 0
    )
  }, by = sequence_name]
  stat_counts <- merge(stat_counts, irf_cluster, by = "sequence_name", all.x = TRUE)
  stat_counts[is.na(max_irf_in_75bp), `:=`(
    max_irf_in_75bp = 0L,
    prop_irf_in_densest_75bp = 0
  )]

  stat_sites_pos <- stat_hits_u[, .(sequence_name, pos = as.integer(floor((start + stop) / 2)), tf = "STAT")]
  irf_sites_pos <- irf_sites[, .(sequence_name, pos = as.integer(floor((start + stop) / 2)), tf = "IRF")]
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
        best_stat <- sum(tt[idx] == "STAT")
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
  stat_counts[is.na(max_stat1stat2_irf_in_75bp), `:=`(
    max_stat1stat2_irf_in_75bp = 0L,
    stat1stat2_in_best_ifn_window = 0L,
    irf_in_best_ifn_window = 0L
  )]

  stat_counts[, ifn_stimulation_potential := fcase(
    stat1stat2_in_best_ifn_window >= 1 & irf_in_best_ifn_window >= 1 & max_stat1stat2_irf_in_75bp >= 3,
    "1_Very High (Synergistic Cluster)",
    stat1stat2_in_best_ifn_window >= 1 & irf_in_best_ifn_window >= 1 & max_stat1stat2_irf_in_75bp < 3,
    "2_High (Synergistic Pair)",
    stat1stat2_in_best_ifn_window >= 3, paste0("3_High (", STAT_TF_ALTID, " Cluster)"),
    irf_in_best_ifn_window >= 3,       "3_High (IRF Cluster)",
    max_stat1stat2_irf_in_75bp == 2 & stat1stat2_in_best_ifn_window == 2, paste0("4_Medium (", STAT_TF_ALTID, " Pair)"),
    max_stat1stat2_irf_in_75bp == 2 & irf_in_best_ifn_window == 2,       "4_Medium (IRF Pair)",
    max_stat1stat2_irf_in_75bp == 1, "5_Low (Isolated Site)",
    default = "6_Unknown"
  )]
  list(summary = stat_counts, stat_sites = stat_hits_u, irf_sites = irf_sites, stat_raw = stat_hits, irf_raw = irf_hits)
}

# ------------------------------------------------------------------------------
# RELEASE: both profiles, site table, QC, output guards, atomic publication
# ------------------------------------------------------------------------------
res <- lapply(names(PROFILES), function(pn) {
  pr <- PROFILES[[pn]]
  r <- summarise_profile(ifn, pr$stat, WINDOW_BP, MIN_FRAC)
  if (!is.null(pr$rename)) setnames(r$summary, names(pr$rename), unname(pr$rename))
  r
})
names(res) <- names(PROFILES)

site_cols <- c("sequence_name", "site_id_in_class", "site_n_hits", "site_start", "site_stop", "motif_id", "motif_alt_id",
               "start", "stop", "strand", "score", "p-value", "pos")
sites <- rbindlist(lapply(names(res), function(pn) {
  one <- function(d, cls) {
    if (!nrow(d)) return(NULL)
    d <- copy(d)
    if (!"pos" %in% names(d)) d[, pos := as.integer(floor((start + stop) / 2))]
    setorder(d, sequence_name, site_start, site_stop)
    d[, site_id_in_class := seq_len(.N), by = sequence_name]
    cbind(data.table(profile = pn, site_class = cls), d[, ..site_cols])
  }
  rbind(one(res[[pn]]$stat_sites, PROFILES[[pn]]$stat), one(res[[pn]]$irf_sites, "IRF"))
}))
setorder(sites, profile, sequence_name, site_class, site_id_in_class)

qc <- list(
  data.table(scope = "source", metric = c("src09_rows", "src09_ltrs", "stat_irf_rows", "stat_irf_exact_duplicate_rows"),
             value = c(n_src_rows, n_src_ltrs, nrow(ifn), n_exact_dup)),
  data.table(scope = "parameters", metric = c("ifn_motif_overlap_collapse_fraction", "ifn_cluster_window_bp"), value = c(MIN_FRAC, WINDOW_BP)),
  data.table(scope = "source_motifs_excluded_by_the_historical_rules",
             metric = alt_counts[grepl("^(stat|irf)", motif_alt_id, ignore.case = TRUE) & !(motif_alt_id %in% stat_ids) &
                                   !grepl(IRF_PATTERN, motif_alt_id, ignore.case = TRUE)][order(motif_alt_id), paste0("rows:", motif_alt_id)],
             value = alt_counts[grepl("^(stat|irf)", motif_alt_id, ignore.case = TRUE) & !(motif_alt_id %in% stat_ids) &
                                  !grepl(IRF_PATTERN, motif_alt_id, ignore.case = TRUE)][order(motif_alt_id), N])
)
for (pn in names(res)) {
  r <- res[[pn]]; s <- r$summary; stat_n <- names(s)[2]; mx <- names(s)[25]
  raw <- rbind(r$stat_raw[, .(motif_alt_id, sequence_name)], r$irf_raw[, .(motif_alt_id, sequence_name)])
  qc[[length(qc) + 1]] <- rbind(
    data.table(scope = pn, metric = c("ltrs", "stat_raw_hits", "stat_sites", "irf_raw_hits", "irf_sites", "ltrs_with_irf_site",
                                      "ltrs_with_stat_window_ge2", "ltrs_with_joint_window_ge2", "ltrs_with_joint_window_ge3"),
               value = c(nrow(s), nrow(r$stat_raw), nrow(r$stat_sites), nrow(r$irf_raw), nrow(r$irf_sites), sum(s$has_irf),
                         sum(s[[names(s)[13]]] >= 2), sum(s[[mx]] >= 2), sum(s[[mx]] >= 3))),
    raw[, .(metric = paste0("raw_hits:", motif_alt_id), value = .N), by = motif_alt_id][order(motif_alt_id), .(metric, value)][, scope := pn],
    raw[, .(metric = paste0("ltrs_with_raw_hit:", motif_alt_id), value = uniqueN(sequence_name)), by = motif_alt_id][order(motif_alt_id), .(metric, value)][, scope := pn],
    rbind(r$stat_sites[, .(motif_alt_id)], r$irf_sites[, .(motif_alt_id)])[, .(metric = paste0("site_representatives:", motif_alt_id), value = .N), by = motif_alt_id][order(motif_alt_id), .(metric, value)][, scope := pn],
    s[, .(metric = paste0("ifn_stimulation_potential:", ifn_stimulation_potential), value = .N), by = ifn_stimulation_potential][order(ifn_stimulation_potential), .(metric, value)][, scope := pn],
    use.names = TRUE)
}
qc <- rbindlist(qc, use.names = TRUE)[, .(scope, metric, value)]

# output guards
for (pn in names(res)) {
  s <- res[[pn]]$summary
  if (anyDuplicated(s$sequence_name) || anyNA(s$sequence_name)) stop(pn, ": duplicate or missing LTR in the summary")
  if (!all(s$sequence_name %in% ltr_id$sequence_name)) stop(pn, ": LTR not present in SRC09")
  if (!all(s$ifn_stimulation_potential != "6_Unknown")) message("[note] ", pn, ": ", sum(s$ifn_stimulation_potential == "6_Unknown"), " 6_Unknown")
  if (any(s[[2]] < 1L)) stop(pn, ": an LTR of the STAT set without a STAT site")
}
if (anyDuplicated(sites[, .(profile, sequence_name, site_class, site_id_in_class)])) stop("duplicate site key")

parent <- dirname(normalizePath(FINAL_OUT, mustWork = FALSE))
dir.create(parent, recursive = TRUE, showWarnings = FALSE)
STAGE <- tempfile(pattern = ".stage_interferon_motifs.", tmpdir = parent)
dir.create(STAGE)
options(error = function() { unlink(STAGE, recursive = TRUE); quit(save = "no", status = 1) })
for (pn in names(res)) fwrite(res[[pn]]$summary, file.path(STAGE, PROFILES[[pn]]$file), sep = "\t")
fwrite(sites, file.path(STAGE, "ifn_motif_sites.tsv"), sep = "\t")
fwrite(qc, file.path(STAGE, "ifn_motif_qc.tsv"), sep = "\t")
fs <- sort(list.files(STAGE, pattern = "\\.tsv$"))
sha <- vapply(fs, function(f) digest::digest(file = file.path(STAGE, f), algo = "sha256"), "")
writeLines(paste0(sha, "  ", fs), file.path(STAGE, "checksums.sha256"))
Sys.chmod(list.files(STAGE, full.names = TRUE), "0644")
if (dir.exists(FINAL_OUT)) unlink(FINAL_OUT, recursive = TRUE)
if (!file.rename(STAGE, FINAL_OUT)) stop("could not publish ", FINAL_OUT)
Sys.chmod(FINAL_OUT, "0755")
for (pn in names(res)) message("[OK] ", pn, ": ", nrow(res[[pn]]$summary), " LTRs -> ", PROFILES[[pn]]$file)
message("[OK] sites ", nrow(sites), "; wrote ", FINAL_OUT)
