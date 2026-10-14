#!/usr/bin/env Rscript
# ==============================================================================
# 07_sparcs step 5 (frozen global order 46, effective order 46), LOCAL, offline:
# U3-R-U5 structural comparison of SPARCS-like LTRs against all other LTRs of the published catalogue.
#
# REFACTOR_SCRIPT of E_sparcs_and_dsrna_and_antisense_ltr/02_analysis_sparcs_ltrs.R (sha256 2581b523...). The
# analysis code is the historical one, verbatim apart from the RELEASE edits; 03_analysis_sparcs_ltrs_detail.R (a
# later variant whose outputs were never saved) is not migrated.
#
#   foreground   unique SPARCS-like LTRs = union of the coding_3UTR and lncRNA_last_exon branches (one row per LTR;
#                an LTR in both branches counts once; transcript rows collapsed by name)
#   structure    SRC08 HERV_LTR_U3_R_U5_catalogue.tsv (560,594 LTRs): status OK / LOW_CONF, role, u3_end_rel,
#                u5_start_rel, confidence. "reconstructed_both" = !is.na(u3_end_rel) & !is.na(u5_start_rel)
#   background   every other catalogue LTR (foreground excluded): group background_all; for the solo comparison,
#                role == "solo" LTRs (solo_background)
#   tests        Fisher (two-sided, conditional MLE OR, 95% CI) on status == "OK" and on reconstructed_both, and a
#                Wilcoxon rank-sum (normal approximation with continuity correction) on confidence, each for all
#                LTRs and for solo LTRs; 6 raw p-values, no multiple-testing correction; no subfamily adjustment
#
# Release changes: the foreground is the order-43 branch table (D033), guarded equal to the historical rule; SRC08 is
# content-verified; the contingency / Wilcoxon group order is fixed explicitly (historically it followed the locale's
# collation: en_US gives background_all first, which is kept); the test objects the script only printed are written
# (u3ru5_tests.tsv); an exact-count QC table is added; no plots (08_figures_supplement); staged atomic output.
# Usage: compare_u3ru5_integrity.R --candidates-dir 07_sparcs/results/candidates --u3r <SRC08> \
#          --source-status 01_sources/SOURCE_STATUS.tsv --outdir <dir>
# ==============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(digest)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("missing argument ", flag)
  args[[i + 1]]
}
CAND_DIR <- get_arg("--candidates-dir"); IN_U3R <- get_arg("--u3r"); SOURCE_STATUS <- get_arg("--source-status"); FINAL_OUT <- get_arg("--outdir")

# ------------------------------------------------------------------
# RELEASE: verified inputs
# ------------------------------------------------------------------
verify_manifest <- function(path) {
  m <- read.table(file.path(dirname(path), "checksums.sha256"), col.names = c("sha", "file"), colClasses = "character")
  want <- m$sha[m$file %in% c(basename(path), paste0("./", basename(path)))]
  if (length(want) != 1 || digest(file = path, algo = "sha256") != want) stop(path, " does not match its accepted manifest")
  message("[ok] ", basename(path), " sha256 ", want, " (manifest)")
}
IN_CANDIDATES <- file.path(CAND_DIR, "ltr_candidates_antisense_last_exon_or_utr3.tsv")
IN_BRANCHES   <- file.path(CAND_DIR, "sparcs_candidate_branches.tsv")
for (p in c(IN_CANDIDATES, IN_BRANCHES)) verify_manifest(p)
if (grepl("hervs_annotation4|ncRNA_layer|HERVs4_v2/analysis/(data|results)/ltr_candidates", CAND_DIR)) stop("order-43 inputs must be the accepted release tables")
src8 <- grep("^SRC08\t", readLines(SOURCE_STATUS), value = TRUE)
src8 <- sub(".*md5 ", "", src8)
want8 <- regmatches(src8, regexpr("^[0-9a-f]{32}", src8))
if (length(want8) != 1 || unname(tools::md5sum(IN_U3R)) != want8) stop(IN_U3R, " is not the registered SRC08 catalogue")
message("[ok] SRC08 md5 ", want8)

# ------------------------------------------------------------------
# Helper (historical)
# ------------------------------------------------------------------
infer_subfamily <- function(x) stringr::str_replace(x, "_merged.*$", "")

# ------------------------------------------------------------------
# Load accepted candidates (RELEASE: order-43 branches) and the structural catalogue
# ------------------------------------------------------------------
cand_all <- read_tsv(IN_CANDIDATES, show_col_types = FALSE) %>% as.data.table()
br  <- read_tsv(IN_BRANCHES, show_col_types = FALSE) %>% as.data.table()
u3r <- read_tsv(IN_U3R,      show_col_types = FALSE) %>% as.data.table()
key <- function(d) paste(d$a_name, d$tx_id, d$feature_class, sep = "\r")
if (!setequal(key(br[branch == "coding_3UTR"]), key(cand_all[gene_type == "protein_coding" & feature_class == "three_prime_UTR"])) ||
    !setequal(key(br[branch == "lncRNA_last_exon"]), key(cand_all[gene_type == "lncRNA" & feature_class == "last_exon"])))
  stop("order-43 branch labels differ from the historical 02 rule (protein_coding & 3' UTR; lncRNA & last exon)")
if (!all(c("name", "status", "role", "u3_end_rel", "u5_start_rel", "confidence", "start", "end") %in% names(u3r)) || anyDuplicated(u3r$name))
  stop("SRC08: unexpected schema or duplicate LTR names")
if (!all(u3r$status %in% c("OK", "LOW_CONF"))) stop("SRC08: unrecognised status value")

cand <- br[, !"branch"]
cand[, subfamily := infer_subfamily(a_name)]
u3r[, subfamily := infer_subfamily(name)]
cand[, name := a_name]

# ------------------------------------------------------------------
# Define the unique SPARCS-like LTR universe (historical)
# ------------------------------------------------------------------
cand_pc  <- cand[br$branch == "coding_3UTR"]
cand_lnc <- cand[br$branch == "lncRNA_last_exon"]
sparcs_like <- rbindlist(list(cand_pc, cand_lnc), use.names = TRUE, fill = TRUE)
sparcs_like_uniq <- sparcs_like[
  ,
  .(
    gene_types      = paste(sort(unique(gene_type)), collapse = ";"),
    feature_classes = paste(sort(unique(feature_class)), collapse = ";"),
    gene_names      = paste(sort(unique(gene_name)), collapse = ";"),
    n_tx            = uniqueN(tx_id),
    n_genes         = uniqueN(gene_id),
    subfamily       = subfamily[1]
  ),
  by = .(name)
]
if (anyDuplicated(sparcs_like_uniq$name)) stop("duplicate LTR in the unique foreground")
missing_ltrs <- setdiff(sparcs_like_uniq$name, u3r$name)

# ------------------------------------------------------------------
# Assign structural status and groups (historical)
# ------------------------------------------------------------------
u3r_dt <- copy(u3r)
u3r_dt[, ltr_len := end - start]
u3r_dt[, has_u5_boundary := !is.na(u5_start_rel)]
u3r_dt[, has_u3_boundary := !is.na(u3_end_rel)]
u3r_dt[, reconstructed_both := has_u3_boundary & has_u5_boundary]
u3r_dt[, R_len := ifelse(reconstructed_both, u5_start_rel - u3_end_rel - 1, NA_real_)]
u3r_dt[, status_simple := fifelse(status == "OK", "OK", "LOW_CONF")]
u3r_dt[, group := "background_all"]
u3r_dt[name %in% sparcs_like_uniq$name, group := "SPARCS_like"]
u3r_dt[, group_solo := fifelse(role == "solo", "solo_background", NA_character_)]
u3r_dt[role == "solo" & name %in% sparcs_like_uniq$name, group_solo := "SPARCS_like"]
sparcs_u3r <- merge(sparcs_like_uniq, u3r_dt, by = "name", all.x = TRUE)

summ <- function(d) d[, .(
  n = .N,
  n_OK = sum(status == "OK", na.rm = TRUE),
  frac_OK = mean(status == "OK", na.rm = TRUE),
  n_both = sum(reconstructed_both, na.rm = TRUE),
  frac_both = mean(reconstructed_both, na.rm = TRUE),
  median_conf = median(confidence, na.rm = TRUE),
  mean_conf = mean(confidence, na.rm = TRUE),
  median_u3 = median(u3_end_rel, na.rm = TRUE),
  median_u5 = median(u5_start_rel, na.rm = TRUE),
  median_R_len = median(R_len, na.rm = TRUE)
)]
summary_by_group <- u3r_dt[, summ(.SD), by = group]
summary_solo <- u3r_dt[role == "solo", summ(.SD), by = group_solo]

# ------------------------------------------------------------------
# Statistical comparison (historical tests; RELEASE: explicit group order = the historical en_US order)
# ------------------------------------------------------------------
LEV_ALL <- c("background_all", "SPARCS_like"); LEV_SOLO <- c("solo_background", "SPARCS_like")
solo <- u3r_dt[role == "solo"]
tab_ok_all   <- table(factor(u3r_dt$group, levels = LEV_ALL), u3r_dt$status == "OK", useNA = "no")
tab_ok_solo  <- table(factor(solo$group_solo, levels = LEV_SOLO), solo$status == "OK", useNA = "no")
tab_both_all <- table(factor(u3r_dt$group, levels = LEV_ALL), u3r_dt$reconstructed_both, useNA = "no")
tab_both_solo <- table(factor(solo$group_solo, levels = LEV_SOLO), solo$reconstructed_both, useNA = "no")
fisher_ok_all   <- fisher.test(tab_ok_all)
fisher_ok_solo  <- fisher.test(tab_ok_solo)
fisher_both_all <- fisher.test(tab_both_all)
fisher_both_solo <- fisher.test(tab_both_solo)
d_all  <- as.data.frame(u3r_dt[group %in% LEV_ALL]);   d_all$group <- factor(d_all$group, levels = LEV_ALL)
d_solo <- as.data.frame(solo[group_solo %in% LEV_SOLO]); d_solo$group_solo <- factor(d_solo$group_solo, levels = LEV_SOLO)
wilcox_conf_all  <- wilcox.test(confidence ~ group, data = d_all, exact = FALSE)
wilcox_conf_solo <- wilcox.test(confidence ~ group_solo, data = d_solo, exact = FALSE)

stats_report <- tibble(
  comparison = c(
    "OK proportion: SPARCS-like vs all",
    "OK proportion: SPARCS-like vs solo",
    "confidence: SPARCS-like vs all",
    "confidence: SPARCS-like vs solo",
    "reconstructed_both: SPARCS-like vs all",
    "reconstructed_both: SPARCS-like vs solo"
  ),
  p_value = c(
    fisher_ok_all$p.value,
    fisher_ok_solo$p.value,
    wilcox_conf_all$p.value,
    wilcox_conf_solo$p.value,
    fisher_both_all$p.value,
    fisher_both_solo$p.value
  )
)

# RELEASE: the printed test objects, in full
ft_row <- function(label, tab, ft) data.table(
  comparison = label, test = "fisher.test two-sided", background_group = rownames(tab)[1],
  background_FALSE = tab[1, "FALSE"], background_TRUE = tab[1, "TRUE"], sparcs_FALSE = tab[2, "FALSE"], sparcs_TRUE = tab[2, "TRUE"],
  sparcs_frac_TRUE = tab[2, "TRUE"] / sum(tab[2, ]), background_frac_TRUE = tab[1, "TRUE"] / sum(tab[1, ]),
  odds_ratio_sparcs_vs_background = unname(ft$estimate), ci95_low = ft$conf.int[1], ci95_high = ft$conf.int[2],
  statistic_W = NA_real_, p_value = ft$p.value)
wx_row <- function(label, wt, d, g) data.table(
  comparison = label, test = "wilcox.test normal approximation, continuity correction (exact = FALSE)", background_group = levels(d[[g]])[1],
  background_FALSE = NA_integer_, background_TRUE = sum(d[[g]] == levels(d[[g]])[1]), sparcs_FALSE = NA_integer_, sparcs_TRUE = sum(d[[g]] == "SPARCS_like"),
  sparcs_frac_TRUE = NA_real_, background_frac_TRUE = NA_real_, odds_ratio_sparcs_vs_background = NA_real_, ci95_low = NA_real_, ci95_high = NA_real_,
  statistic_W = unname(wt$statistic), p_value = wt$p.value)
tests <- rbindlist(list(
  ft_row("OK proportion: SPARCS-like vs all", tab_ok_all, fisher_ok_all),
  ft_row("OK proportion: SPARCS-like vs solo", tab_ok_solo, fisher_ok_solo),
  wx_row("confidence: SPARCS-like vs all", wilcox_conf_all, d_all, "group"),
  wx_row("confidence: SPARCS-like vs solo", wilcox_conf_solo, d_solo, "group_solo"),
  ft_row("reconstructed_both: SPARCS-like vs all", tab_both_all, fisher_both_all),
  ft_row("reconstructed_both: SPARCS-like vs solo", tab_both_solo, fisher_both_solo)
))

subfam_tab <- u3r_dt[
  ,
  .(
    n_total = .N,
    n_OK = sum(status == "OK", na.rm = TRUE),
    n_sparcs = sum(group == "SPARCS_like", na.rm = TRUE),
    n_sparcs_OK = sum(group == "SPARCS_like" & status == "OK", na.rm = TRUE)
  ),
  by = subfamily
][order(-n_sparcs)]

sparcs_ranked <- sparcs_u3r %>%
  as_tibble() %>%
  arrange(desc(status == "OK"), desc(confidence), desc(reconstructed_both), desc(n_tx))

# ------------------------------------------------------------------
# RELEASE: exact-count QC (branches are descriptive only: history tests the union)
# ------------------------------------------------------------------
qrow <- function(analysis, branch, population, names_) {
  d <- u3r_dt[name %in% names_]
  data.table(analysis = analysis, branch = branch, population = population, total_LTRs = length(unique(names_)),
             OK = sum(d$status == "OK"), LOW_CONF = sum(d$status == "LOW_CONF"), missing = length(setdiff(unique(names_), u3r_dt$name)),
             fraction_OK = sum(d$status == "OK") / nrow(d), reconstructed_both = sum(d$reconstructed_both), solo = sum(d$role == "solo"))
}
fg <- sparcs_like_uniq$name
qc <- rbindlist(list(
  qrow("tested (historical)", "union", "SPARCS-like foreground", fg),
  qrow("tested (historical)", "none", "background_all (catalogue minus foreground)", setdiff(u3r_dt$name, fg)),
  qrow("tested (historical)", "union", "SPARCS-like solo foreground", u3r_dt[group_solo == "SPARCS_like", name]),
  qrow("tested (historical)", "none", "solo_background", u3r_dt[group_solo == "solo_background", name]),
  qrow("descriptive only", "coding_3UTR", "coding_3UTR LTRs", unique(cand_pc$name)),
  qrow("descriptive only", "lncRNA_last_exon", "lncRNA_last_exon LTRs", unique(cand_lnc$name)),
  qrow("descriptive only", "both", "LTRs in both branches", intersect(cand_pc$name, cand_lnc$name)),
  qrow("descriptive only", "all", "SRC08 catalogue", u3r_dt$name)
))
qc[, rows_collapsed := c(nrow(sparcs_like), NA, NA, NA, nrow(cand_pc), nrow(cand_lnc), NA, NA)]
if (qc$total_LTRs[1] + qc$total_LTRs[2] != nrow(u3r_dt) || length(missing_ltrs)) stop("foreground / background do not partition the catalogue")

# ------------------------------------------------------------------
# Source-data export (staged, atomic)
# ------------------------------------------------------------------
parent <- dirname(normalizePath(FINAL_OUT, mustWork = FALSE))
dir.create(parent, recursive = TRUE, showWarnings = FALSE)
STAGE <- tempfile(pattern = ".stage_u3ru5.", tmpdir = parent)
dir.create(STAGE)
options(error = function() { unlink(STAGE, recursive = TRUE); quit(save = "no", status = 1) })
fwrite(sparcs_u3r, file.path(STAGE, "SPARCS_like_candidates_with_U3R_U5_catalogue.tsv"), sep = "\t")
fwrite(summary_by_group, file.path(STAGE, "summary_vs_all_LTRs.tsv"), sep = "\t")
fwrite(summary_solo,     file.path(STAGE, "summary_vs_soloLTRs.tsv"), sep = "\t")
write_tsv(stats_report, file.path(STAGE, "stats_report.tsv"))
fwrite(subfam_tab, file.path(STAGE, "subfamily_summary.tsv"), sep = "\t")
write_tsv(sparcs_ranked, file.path(STAGE, "SPARCS_like_candidates_ranked.tsv"))
fwrite(tests, file.path(STAGE, "u3ru5_tests.tsv"), sep = "\t")
fwrite(qc, file.path(STAGE, "u3ru5_qc_counts.tsv"), sep = "\t")
fs <- sort(list.files(STAGE, pattern = "\\.tsv$"))
writeLines(paste0(vapply(fs, function(f) digest(file = file.path(STAGE, f), algo = "sha256"), ""), "  ", fs), file.path(STAGE, "checksums.sha256"))
Sys.chmod(list.files(STAGE, full.names = TRUE), "0644")
if (dir.exists(FINAL_OUT)) unlink(FINAL_OUT, recursive = TRUE)
if (!file.rename(STAGE, FINAL_OUT)) stop("could not publish ", FINAL_OUT)
Sys.chmod(FINAL_OUT, "0755")
message("[OK] unique SPARCS-like LTRs ", length(fg), " (missing from SRC08: ", length(missing_ltrs), "); OK ", qc$OK[1], " / ", qc$total_LTRs[1],
        " vs background ", qc$OK[2], " / ", qc$total_LTRs[2], "; Fisher OR ", signif(unname(fisher_ok_all$estimate), 6), " p ", signif(fisher_ok_all$p.value, 6))
message("[OK] wrote ", FINAL_OUT)
