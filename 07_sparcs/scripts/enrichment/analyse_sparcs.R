#!/usr/bin/env Rscript
# ==============================================================================
# 07_sparcs step 4 (frozen global order 45, effective order 45), LOCAL, offline:
# SPARCS-like candidate analysis, IFN motif integration, GO / MSigDB C7 enrichment and TF-motif analysis.
#
# REFACTOR_SCRIPT of E_sparcs_and_dsrna_and_antisense_ltr/01_analysis_sparcs.R (sha256 a017bced...), split by
# concern. The analytical code of each block is the historical one, verbatim apart from the RELEASE edits:
#
#   A inputs     order-43 candidate table + branch table (D033), order-44 IFN summaries, SRC02 GTF, SRC09 FIMO,
#                SRC18 C7 (D034); every input is checksum-verified. Nothing upstream is recomputed: no antisense
#                geometry, no terminal-exon / 3' UTR assignment, no STAT/IRF filtering, collapse or window.
#   B IFN join   candidates x ISGF3 / GAF summaries on a_name = sequence_name; tiers; priority tables
#   C IFN tests  coding_3UTR LTRs vs every other SRC09 LTR, Fisher per tier mode x STAT profile, BH over the 10
#                (historically only printed; written here)
#   D GO         BP ORA, coding_3UTR gene symbols vs all GENCODE v48 protein_coding gene symbols (both mapped to
#                ENTREZ by org.Hs.eg.db); BH; IFN-high subset
#   E MSigDB C7  enricher() on the mapped symbols, TERM2GENE / TERM2NAME from SRC18 (MSigDB 2026.1.Hs; D034)
#   F TF motifs  per-motif Fisher: branch LTRs vs all other SRC09 LTRs (presence of >= 1 FIMO row), BH per branch;
#                GO BP over the Enriched / Depleted TF symbols
#   G Figure 5B  unique a_name and unique gene_id per branch (source data only; plotting is 08_figures_supplement)
#   H QC         units, branch union / intersection, D033 lineage
#
# Release changes: branches come from the order-43 branch table (guarded equal to the historical rule); C7 is read
# from SRC18 instead of a live msigdbr call (D034); thresholds from parameters.yaml; the tables the historical
# script only printed are written; no plots and no xlsx (08_figures_supplement); staged atomic output.
# Usage: analyse_sparcs.R --candidates-dir 07_sparcs/results/candidates --ifn-dir 07_sparcs/results/interferon_motifs
#          --gtf <SRC02> --source-status 01_sources/SOURCE_STATUS.tsv --fimo <SRC09> --archive-members 01_sources/manifests/archive_members.tsv
#          --msigdb-dir 01_sources/frozen/msigdb --snapshot-checksums 01_sources/manifests/manual_snapshot_checksums.tsv
#          --params config/parameters.yaml --outdir <dir>
# ==============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(AnnotationDbi)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(digest)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("missing argument ", flag)
  args[[i + 1]]
}
CAND_DIR <- get_arg("--candidates-dir"); IFN_DIR <- get_arg("--ifn-dir")
IN_GTF <- get_arg("--gtf"); SOURCE_STATUS <- get_arg("--source-status")
FIMO <- get_arg("--fimo"); MEMBERS <- get_arg("--archive-members")
MSIG_DIR <- get_arg("--msigdb-dir"); SNAP <- get_arg("--snapshot-checksums")
PARAMS <- get_arg("--params"); FINAL_OUT <- get_arg("--outdir")

# ------------------------------------------------------------------------------
# A. RELEASE: parameters and verified inputs
# ------------------------------------------------------------------------------
read_sparcs_param <- function(path, key) {
  txt <- readLines(path, warn = FALSE)
  s <- grep("^sparcs:", txt)
  blk <- txt[(s + 1):length(txt)]
  end <- grep("^\\S", blk)
  if (length(end)) blk <- blk[seq_len(end[1] - 1)]
  hit <- grep(paste0("^\\s+", key, ":\\s*[0-9.eE+-]+\\s*(#.*)?$"), blk, value = TRUE)
  if (length(s) != 1 || length(hit) != 1) stop("parameter sparcs.", key, " not found")
  as.numeric(sub(paste0("^\\s+", key, ":\\s*([0-9.eE+-]+).*$"), "\\1", hit))
}
P_CUT <- read_sparcs_param(PARAMS, "enrichment_pvalue_cutoff"); Q_CUT <- read_sparcs_param(PARAMS, "enrichment_qvalue_cutoff")
IFN_GO_MIN <- read_sparcs_param(PARAMS, "ifn_high_go_min_genes"); MOTIF_PADJ <- read_sparcs_param(PARAMS, "tf_motif_padj_cutoff")
TF_P_CUT <- read_sparcs_param(PARAMS, "tf_go_pvalue_cutoff"); TF_Q_CUT <- read_sparcs_param(PARAMS, "tf_go_qvalue_cutoff")
message("[ok] parameters: GO/C7 p ", P_CUT, " q ", Q_CUT, "; IFN-high GO min ", IFN_GO_MIN, "; motif p_adj ", MOTIF_PADJ,
        "; TF GO p ", TF_P_CUT, " q ", TF_Q_CUT)

verify_manifest <- function(path) {
  m <- read.table(file.path(dirname(path), "checksums.sha256"), col.names = c("sha", "file"), colClasses = "character")
  want <- m$sha[m$file %in% c(basename(path), paste0("./", basename(path)))]
  if (length(want) != 1 || digest(file = path, algo = "sha256") != want) stop(path, " does not match its accepted manifest")
  message("[ok] ", basename(path), " sha256 ", want, " (manifest)")
}
IN_CANDIDATES <- file.path(CAND_DIR, "ltr_candidates_antisense_last_exon_or_utr3.tsv")
IN_BRANCHES   <- file.path(CAND_DIR, "sparcs_candidate_branches.tsv")
IN_IFN_ISGF3  <- file.path(IFN_DIR, "LTR_IFN_STAT1STAT2_IRF_summary.tsv")   # alpha/beta proxy
IN_IFN_GAF    <- file.path(IFN_DIR, "LTR_IFN_STAT1_summary.tsv")            # gamma proxy
for (p in c(IN_CANDIDATES, IN_BRANCHES, IN_IFN_ISGF3, IN_IFN_GAF)) verify_manifest(p)
if (any(grepl("hervs_annotation4|ncRNA_layer|HERVs4_v2/analysis/(data|results)/ltr_candidates", c(CAND_DIR, IFN_DIR))))
  stop("order-43/44 inputs must be the accepted release tables, not historical copies")

src2 <- grep("^SRC02\t", readLines(SOURCE_STATUS), value = TRUE)
src2 <- sub(".*sha256 of the uncompressed GTF", "", src2)
want2 <- regmatches(src2, regexpr("[0-9a-f]{64}", src2))
if (length(want2) != 1 || digest(file = IN_GTF, algo = "sha256") != want2) stop(IN_GTF, " is not the registered SRC02 GTF")
message("[ok] SRC02 GTF sha256 ", want2)
m9 <- fread(MEMBERS, sep = "\t", colClasses = "character")
want9 <- m9[get(names(m9)[1]) == "SRC09"][[5]]
if (length(want9) != 1 || unname(tools::md5sum(FIMO)) != want9) stop(FIMO, " is not the registered SRC09 member")
message("[ok] SRC09 md5 ", want9)
snap <- fread(SNAP, sep = "\t", colClasses = "character")
T2G_FILE <- file.path(MSIG_DIR, "msigdb_2026.1.Hs_C7_term2gene.tsv"); T2N_FILE <- file.path(MSIG_DIR, "msigdb_2026.1.Hs_C7_term2name.tsv")
for (f in c(T2G_FILE, T2N_FILE)) {
  w <- snap[source_id == "SRC18" & filename == basename(f), sha256]
  if (length(w) != 1 || digest(file = f, algo = "sha256") != w) stop(f, " is not the registered SRC18 snapshot")
  message("[ok] SRC18 ", basename(f), " sha256 ", w)
}

# historical helpers (verbatim)
infer_subfamily <- function(x) str_replace(x, "_merged.*$", "")
tier_from_ifn_label <- function(x) {
  case_when(
    str_starts(x, "1_Very High") ~ "VeryHigh",
    str_starts(x, "2_High")      ~ "High",
    str_starts(x, "3_High")      ~ "High",
    str_starts(x, "4_Medium")    ~ "Medium",
    str_starts(x, "5_Low")       ~ "Low",
    TRUE                         ~ "None"
  )
}
keep_ifn_cols <- function(dt, prefix) {
  cols <- c(
    "sequence_name", "subfamily",
    "n_irf", "has_irf",
    "max_irf_in_75bp",
    "ifn_stimulation_potential",
    "best_ifn_win_start", "best_ifn_win_end"
  )
  stat_cols <- intersect(
    c("n_stat1stat2", "max_stat1stat2_in_75bp", "max_stat1stat2_irf_in_75bp",
      "stat1stat2_in_best_ifn_window", "irf_in_best_ifn_window",
      "best_p_stat1stat2", "best_score_stat1stat2", "stat1stat2_sites_per_kb"),
    names(dt)
  )
  stat1_cols <- intersect(
    c("n_stat1", "max_stat1_in_75bp", "max_stat1_irf_in_75bp",
      "stat1_in_best_ifn_window", "irf_in_best_ifn_window",
      "best_p_stat1", "best_score_stat1", "stat_sites_per_kb"),
    names(dt)
  )
  cols2 <- unique(c(cols, stat_cols, stat1_cols))
  out <- dt[, ..cols2]
  to_rename <- setdiff(names(out), "sequence_name")
  setnames(out, to_rename, paste0(prefix, "__", to_rename))
  out
}

cand <- read_tsv(IN_CANDIDATES, show_col_types = FALSE) %>% as.data.table()
branches <- fread(IN_BRANCHES, sep = "\t")
cand[, subfamily := infer_subfamily(a_name)]

# RELEASE (D033): the two branches are the order-43 labels; they must equal the historical rule on this table
key <- c("a_name", "tx_id", "feature_class")
bkey <- function(d) do.call(paste, c(d[, ..key], sep = "\r"))
coding_keys <- bkey(branches[branch == "coding_3UTR"]); lnc_keys <- bkey(branches[branch == "lncRNA_last_exon"])
if (!setequal(coding_keys, bkey(cand[gene_type == "protein_coding" & feature_class == "three_prime_UTR"])) ||
    !setequal(lnc_keys, bkey(cand[gene_type == "lncRNA"])) || length(intersect(coding_keys, lnc_keys)))
  stop("order-43 branch labels differ from the historical branch rule on the candidate table")
cand_pc  <- cand[bkey(cand) %chin% coding_keys]
cand_lnc <- cand[bkey(cand) %chin% lnc_keys]
message("[INFO] Candidate rows: ", nrow(cand), "; coding_3UTR rows ", nrow(cand_pc), "; lncRNA_last_exon rows ", nrow(cand_lnc))

# ------------------------------------------------------------------------------
# B. IFN integration (historical; the order-44 summaries are joined, never recomputed)
# ------------------------------------------------------------------------------
ifn_ab <- fread(IN_IFN_ISGF3)
ifn_g  <- fread(IN_IFN_GAF)
ifn_ab_small <- keep_ifn_cols(ifn_ab, "ifn_ab")
ifn_g_small  <- keep_ifn_cols(ifn_g,  "ifn_g")
joined <- copy(cand)
setnames(joined, "a_name", "sequence_name_tmp")
joined <- merge(joined, ifn_ab_small, by.x = "sequence_name_tmp", by.y = "sequence_name", all.x = TRUE)
joined <- merge(joined, ifn_g_small,  by.x = "sequence_name_tmp", by.y = "sequence_name", all.x = TRUE)
setnames(joined, "sequence_name_tmp", "a_name")

joined[, ifn_ab__tier := tier_from_ifn_label(ifn_ab__ifn_stimulation_potential)]
joined[, ifn_g__tier  := tier_from_ifn_label(ifn_g__ifn_stimulation_potential)]
joined[, ifn_ab__pos := ifn_ab__tier %in% c("VeryHigh","High","Medium")]
joined[, ifn_g__pos  := ifn_g__tier  %in% c("VeryHigh","High","Medium")]
joined[, ifn_any_tier := fcase(
  ifn_ab__pos & ifn_g__pos, "Both (Type I/III + Type II)",
  ifn_ab__pos,             "Type I/III-like (STAT1::STAT2)",
  ifn_g__pos,              "Type II-like (STAT1)",
  default = "Not IFN-high"
)]
joined[, ifn_ab__synergy := !is.na(ifn_ab__ifn_stimulation_potential) &
         ifn_ab__ifn_stimulation_potential %chin% c("1_Very High (Synergistic Cluster)", "2_High (Synergistic Pair)")]
joined[, ifn_g__synergy := !is.na(ifn_g__ifn_stimulation_potential) &
         ifn_g__ifn_stimulation_potential %chin% c("1_Very High (Synergistic Cluster)", "2_High (Synergistic Pair)")]
joined[, ifn_any_synergy := ifn_ab__synergy | ifn_g__synergy]

cols_keep <- c(
  "a_name","subfamily","gene_name","gene_type","feature_class","tx_id","tx_type",
  "overlap_bp_total","overlap_frac_of_ltr",
  "ifn_ab__tier","ifn_ab__ifn_stimulation_potential",
  "ifn_ab__n_stat1stat2","ifn_ab__stat1stat2_sites_per_kb","ifn_ab__stat1stat2_in_best_ifn_window",
  "ifn_ab__n_irf","ifn_ab__irf_in_best_ifn_window","ifn_ab__max_stat1stat2_irf_in_75bp",
  "ifn_g__tier","ifn_g__ifn_stimulation_potential",
  "ifn_g__n_stat1","ifn_g__stat1_in_best_ifn_window",
  "ifn_g__n_irf","ifn_g__irf_in_best_ifn_window","ifn_g__max_stat1_irf_in_75bp",
  "ifn_any_tier","ifn_any_synergy"
)
cols_keep <- intersect(cols_keep, names(joined))
priority <- function(rows) {
  pr <- rows[, ..cols_keep]
  pr[, sort_ab_density := suppressWarnings(as.numeric(ifn_ab__stat1stat2_sites_per_kb))]
  pr[, sort_g_nstat    := suppressWarnings(as.numeric(ifn_g__n_stat1))]
  pr[is.na(sort_ab_density), sort_ab_density := -Inf]
  pr[is.na(sort_g_nstat),    sort_g_nstat    := -Inf]
  setorder(pr, ifn_any_tier, -ifn_any_synergy, -sort_ab_density, -sort_g_nstat, gene_name)
  pr[, c("sort_ab_density","sort_g_nstat") := NULL]
  pr
}
pc_priority  <- priority(joined[bkey(joined) %chin% coding_keys])
lnc_priority <- priority(joined[bkey(joined) %chin% lnc_keys])

# ------------------------------------------------------------------------------
# D / E. GO BP and MSigDB C7 ORA of the coding_3UTR genes (historical)
# ------------------------------------------------------------------------------
gtf_genes <- fread(IN_GTF, sep = "\t", header = FALSE, quote = "",
                   col.names = c("seqname","source","feature","start","end","score","strand","frame","attribute"))[feature == "gene"]
get_attr <- function(x, key) {
  m <- str_match(x, paste0(key, ' "([^"]+)"'))
  m[,2]
}
gtf_genes[, gene_type := get_attr(attribute, "gene_type")]
gtf_genes[, gene_name := get_attr(attribute, "gene_name")]
universe_symbols <- unique(na.omit(gtf_genes[gene_type == "protein_coding", gene_name]))
rm(gtf_genes); invisible(gc())
message("[INFO] Universe protein_coding symbols: ", length(universe_symbols))

candidate_pc_genes <- unique(na.omit(joined[bkey(joined) %chin% coding_keys, gene_name]))
candidate_pc_genes <- intersect(candidate_pc_genes, universe_symbols)
candidate_pc_ifnhigh <- unique(na.omit(joined[bkey(joined) %chin% coding_keys & ifn_any_tier != "Not IFN-high", gene_name]))
candidate_pc_ifnhigh <- intersect(candidate_pc_ifnhigh, universe_symbols)
message("[INFO] PC SPARCS-like genes: ", length(candidate_pc_genes), "; IFN-high: ", length(candidate_pc_ifnhigh))

cand_map <- suppressMessages(suppressWarnings(bitr(candidate_pc_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)))
univ_map <- suppressMessages(suppressWarnings(bitr(universe_symbols,   fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)))
cand_entrez <- unique(cand_map$ENTREZID)
univ_entrez <- unique(univ_map$ENTREZID)

ego_bp_all <- enrichGO(gene = cand_entrez, universe = univ_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                       pAdjustMethod = "BH", pvalueCutoff = P_CUT, qvalueCutoff = Q_CUT, readable = TRUE)
bp_all <- as.data.frame(ego_bp_all)
bp_ifn <- NULL
if (length(candidate_pc_ifnhigh) >= IFN_GO_MIN) {
  cand_ifn_map <- suppressMessages(suppressWarnings(bitr(candidate_pc_ifnhigh, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)))
  ego_bp_ifn <- enrichGO(gene = unique(cand_ifn_map$ENTREZID), universe = univ_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                         ont = "BP", pAdjustMethod = "BH", pvalueCutoff = P_CUT, qvalueCutoff = Q_CUT, readable = TRUE)
  bp_ifn <- as.data.frame(ego_bp_ifn)
}

# RELEASE (D034): C7 from the frozen SRC18 snapshot (MSigDB 2026.1.Hs), identical to msigdbr's TERM2GENE / TERM2NAME
term2gene <- fread(T2G_FILE, sep = "\t", colClasses = "character", quote = "")
term2name <- fread(T2N_FILE, sep = "\t", colClasses = "character", quote = "")
ego_c7_all <- enricher(gene = cand_map$SYMBOL, universe = univ_map$SYMBOL, TERM2GENE = term2gene, TERM2NAME = term2name,
                       pAdjustMethod = "BH", pvalueCutoff = P_CUT, qvalueCutoff = Q_CUT)
c7_all <- as.data.frame(ego_c7_all)

# ------------------------------------------------------------------------------
# C / F. SRC09 LTR universe: IFN tier tests and per-motif enrichment (historical)
# ------------------------------------------------------------------------------
fimo_all <- fread(FIMO, select = c("sequence_name", "motif_alt_id"))
univ_ltrs <- unique(fimo_all$sequence_name)
fimo_hits <- unique(fimo_all[, .(sequence_name, motif_alt_id)])
rm(fimo_all); invisible(gc())

sparcs_ltrs <- unique(pc_priority$a_name)
tier_num <- function(x) {
  as.integer(str_match(x, "^\\s*([1-5])\\s*[_\\- ]")[,2])
}
is_synergy <- function(x) {
  str_detect(x, regex("^\\s*[12]\\s*[_\\- ]", ignore_case = TRUE))
}
pos_ltrs_from_summary <- function(dt, mode = c("non_low", "synergy_only", "tier1_only", "tier1_2_3", "tier1_2_3_4")) {
  mode <- match.arg(mode)
  t <- dt$ifn_stimulation_potential
  tn <- tier_num(t)
  keep <- switch(
    mode,
    non_low     = !is.na(tn) & tn %in% 1:4,
    synergy_only = is_synergy(t),
    tier1_only  = !is.na(tn) & tn == 1,
    tier1_2_3   = !is.na(tn) & tn %in% 1:3,
    tier1_2_3_4 = !is.na(tn) & tn %in% 1:4
  )
  unique(dt$sequence_name[keep])
}
fisher_enrich_ltrs <- function(fg_ltrs, pos_ltrs, universe_ltrs, label) {
  fg_ltrs <- intersect(unique(fg_ltrs), universe_ltrs)
  other_ltrs <- setdiff(universe_ltrs, fg_ltrs)
  pos_ltrs <- intersect(unique(pos_ltrs), universe_ltrs)
  a <- length(intersect(fg_ltrs, pos_ltrs))
  b <- length(fg_ltrs) - a
  c <- length(intersect(other_ltrs, pos_ltrs))
  d <- length(other_ltrs) - c
  mat <- matrix(c(a,b,c,d), nrow=2, byrow=TRUE, dimnames=list(c("FG","Other"), c("Pos","Neg")))
  ft <- fisher.test(mat)
  data.table(test = label, fg_total = length(fg_ltrs), fg_with = a, other_total = length(other_ltrs), other_with = c,
             odds_ratio = unname(ft$estimate), p_value = ft$p.value)
}
modes <- c("non_low","synergy_only","tier1_only","tier1_2_3","tier1_2_3_4")
res_all <- rbindlist(lapply(modes, function(mode) {
  pos_stat1    <- pos_ltrs_from_summary(ifn_g,  mode = mode)
  pos_stat1st2 <- pos_ltrs_from_summary(ifn_ab, mode = mode)
  rbindlist(list(
    fisher_enrich_ltrs(sparcs_ltrs, pos_stat1, univ_ltrs, paste0("STAT1 | mode=", mode)),
    fisher_enrich_ltrs(sparcs_ltrs, pos_stat1st2, univ_ltrs, paste0("STAT1::STAT2 | mode=", mode))
  ))
}))
res_all[, p_adj := p.adjust(p_value, method = "BH")]

fg_pc  <- unique(cand_pc$a_name)
fg_lnc <- unique(cand_lnc$a_name)
bg_pc  <- setdiff(univ_ltrs, fg_pc)
bg_lnc <- setdiff(univ_ltrs, fg_lnc)
motif_enrichment <- function(fg_ltrs, bg_ltrs, hits) {
  motifs <- sort(unique(hits$motif_alt_id))
  res <- lapply(motifs, function(m) {
    ltrs_with_m <- unique(hits[motif_alt_id == m, sequence_name])
    a <- length(intersect(fg_ltrs, ltrs_with_m))
    b <- length(fg_ltrs) - a
    c <- length(intersect(bg_ltrs, ltrs_with_m))
    d <- length(bg_ltrs) - c
    mat <- matrix(c(a,b,c,d), nrow=2, byrow=TRUE)
    ft <- suppressWarnings(fisher.test(mat))
    data.table(motif_alt_id = m, fg_with = a, fg_total = length(fg_ltrs), bg_with = c, bg_total = length(bg_ltrs),
               odds_ratio = unname(ft$estimate), p_value = ft$p.value)
  })
  res <- rbindlist(res)
  res[, p_adj := p.adjust(p_value, method="BH")]
  res[, log2_OR := log2(pmax(odds_ratio, .Machine$double.eps))]
  res[]
}
res_pc  <- motif_enrichment(fg_pc, bg_pc, fimo_hits)
res_lnc <- motif_enrichment(fg_lnc, bg_lnc, fimo_hits)
add_labels <- function(dt) {
  dt[, sig := fcase(
    p_adj < MOTIF_PADJ & odds_ratio > 1, "Enriched",
    p_adj < MOTIF_PADJ & odds_ratio < 1, "Depleted",
    default = "NS"
  )]
}
res_pc  <- add_labels(res_pc)
res_lnc <- add_labels(res_lnc)

explode_and_normalize <- function(x) {
  if (is.null(x) || length(x) == 0) return(character(0))
  x <- as.character(x)
  x <- x[!is.na(x)]
  if (!length(x)) return(character(0))
  parts <- unlist(strsplit(x, "::", fixed = TRUE), use.names = FALSE)
  parts <- trimws(parts)
  parts <- toupper(parts)
  parts <- sub("\\s+\\(.*\\)$", "", parts)
  parts <- gsub("[^A-Z0-9_\\-]", "", parts)
  parts <- parts[nchar(parts) > 0]
  unique(parts)
}
map_to_hgnc <- function(ids_upper) {
  ids_upper <- unique(ids_upper)
  if (length(ids_upper) == 0) return(character(0))
  sym_ok <- tryCatch(
    AnnotationDbi::mapIds(org.Hs.eg.db, keys = ids_upper, keytype = "SYMBOL", column = "SYMBOL", multiVals = "first"),
    error = function(e) setNames(rep(NA_character_, length(ids_upper)), ids_upper))
  need_alias <- names(sym_ok)[is.na(sym_ok)]
  alias_map <- if (length(need_alias)) {
    tryCatch(
      AnnotationDbi::mapIds(org.Hs.eg.db, keys = need_alias, keytype = "ALIAS", column = "SYMBOL", multiVals = "first"),
      error = function(e) setNames(rep(NA_character_, length(need_alias)), need_alias))
  } else {
    setNames(character(0), character(0))
  }
  out <- sym_ok
  out[is.na(out)] <- alias_map[names(out)[is.na(out)]]
  out <- unique(stats::na.omit(unname(out)))
  out
}
tf_go <- list()
run_tf_go <- function(res_dt, foreground_sig, out_prefix) {
  fg_raw <- unique(res_dt[sig == foreground_sig, motif_alt_id])
  bg_raw <- unique(res_dt$motif_alt_id)
  fg_norm <- explode_and_normalize(fg_raw)
  bg_norm <- explode_and_normalize(bg_raw)
  tf_symbols  <- suppressMessages(map_to_hgnc(fg_norm))
  tf_universe <- suppressMessages(map_to_hgnc(bg_norm))
  pre <- paste0(out_prefix, "_", tolower(foreground_sig))
  tf_go[[paste0(pre, "_TF_raw_ids.tsv")]] <<- data.table(raw_motif_id = unique(fg_raw))
  tf_go[[paste0(pre, "_TF_normalized_ids.tsv")]] <<- data.table(normalized_motif_id = unique(fg_norm))
  tf_go[[paste0(pre, "_TF_mapped_symbols.tsv")]] <<- data.table(mapped_symbol = unique(tf_symbols))
  unmapped_fg <- setdiff(unique(fg_norm), unique(tf_symbols))
  if (length(unmapped_fg)) tf_go[[paste0(pre, "_TF_unmapped_ids.tsv")]] <<- data.table(unmapped_motif_id = unmapped_fg)
  if (length(tf_symbols) == 0L || length(tf_universe) == 0L) return(invisible(NULL))
  ego_tf <- tryCatch(
    enrichGO(gene = tf_symbols, universe = tf_universe, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
             pAdjustMethod = "BH", pvalueCutoff = TF_P_CUT, qvalueCutoff = TF_Q_CUT, readable = TRUE),
    error = function(e) NULL)
  if (is.null(ego_tf) || NROW(as.data.frame(ego_tf)) == 0) return(invisible(NULL))
  tf_go[[paste0(out_prefix, "_GO_BP_over_", tolower(foreground_sig), "_TFsymbols.tsv")]] <<- as.data.table(as.data.frame(ego_tf))
  invisible(NULL)
}
run_tf_go(res_pc,  "Enriched", "sparcs_pc_enriched")
run_tf_go(res_pc,  "Depleted", "sparcs_pc_depleted")
run_tf_go(res_lnc, "Enriched", "sparcs_lnc_enriched")
run_tf_go(res_lnc, "Depleted", "sparcs_lnc_depleted")

# ------------------------------------------------------------------------------
# G. Figure 5B source data (historical counts; plot labels without the newline, plotting in module 08)
# ------------------------------------------------------------------------------
plot_df <- data.table(
  branch = c("coding_3UTR", "lncRNA_last_exon"),
  class = c("Protein-coding (3' UTR)", "lncRNA (last exon)"),
  n_insertions = c(uniqueN(cand_pc$a_name), uniqueN(cand_lnc$a_name)),
  n_affected = c(uniqueN(cand_pc$gene_id), uniqueN(cand_lnc$gene_id))
)
plot_df[, insertion_pct := n_insertions / sum(n_insertions)]   # historical: denominators are branch sums
plot_df[, affected_pct  := n_affected   / sum(n_affected)]

# ------------------------------------------------------------------------------
# H. RELEASE: QC (units never mixed) and gene lists behind every enrichment
# ------------------------------------------------------------------------------
u <- function(scope, d) data.table(scope = scope,
  metric = c("rows", "unique_ltrs", "unique_genes_gene_id", "unique_genes_gene_name", "unique_ltr_x_gene", "unique_transcripts"),
  value = c(nrow(d), uniqueN(d$a_name), uniqueN(d$gene_id), uniqueN(d$gene_name), nrow(unique(d[, .(a_name, gene_id)])), uniqueN(d$tx_id)))
qc <- rbindlist(list(
  u("all_terminal_candidates", cand), u("coding_3UTR", cand_pc), u("lncRNA_last_exon", cand_lnc),
  data.table(scope = "branches", metric = c("unique_ltrs_union", "unique_ltrs_branch_sum", "shared_ltrs", "shared_genes_gene_id"),
             value = c(length(union(fg_pc, fg_lnc)), length(fg_pc) + length(fg_lnc), length(intersect(fg_pc, fg_lnc)),
                       length(intersect(cand_pc$gene_id, cand_lnc$gene_id)))),
  data.table(scope = "all_terminal_candidates_genes_by_gene_type",
             metric = c("lncRNA", "protein_coding", "other"),
             value = c(uniqueN(cand[gene_type == "lncRNA", gene_id]), uniqueN(cand[gene_type == "protein_coding", gene_id]),
                       uniqueN(cand[!gene_type %in% c("lncRNA", "protein_coding"), gene_id]))),
  data.table(scope = "enrichment_universes",
             metric = c("go_universe_pc_symbols", "go_universe_entrez", "go_foreground_symbols", "go_foreground_entrez",
                        "go_foreground_ifn_high_symbols", "c7_foreground_symbols", "c7_universe_symbols", "src09_ltr_universe",
                        "ifn_test_foreground_ltrs"),
             value = c(length(universe_symbols), length(univ_entrez), length(candidate_pc_genes), length(cand_entrez),
                       length(candidate_pc_ifnhigh), length(cand_map$SYMBOL), length(univ_map$SYMBOL), length(univ_ltrs),
                       length(sparcs_ltrs)))
))
gene_lists <- rbindlist(list(
  data.table(list = "go_c7_foreground_symbol_entrez", symbol = cand_map$SYMBOL, entrez = cand_map$ENTREZID),
  data.table(list = "go_foreground_symbols_unmapped", symbol = setdiff(candidate_pc_genes, cand_map$SYMBOL), entrez = NA_character_),
  data.table(list = "go_c7_universe_symbol_entrez", symbol = univ_map$SYMBOL, entrez = univ_map$ENTREZID),
  data.table(list = "go_universe_symbols_unmapped", symbol = setdiff(universe_symbols, univ_map$SYMBOL), entrez = NA_character_),
  data.table(list = "go_foreground_ifn_high_symbols", symbol = sort(candidate_pc_ifnhigh), entrez = NA_character_)
))
setorder(gene_lists, list, symbol, entrez)

# ------------------------------------------------------------------------------
# Output (staged, atomic)
# ------------------------------------------------------------------------------
parent <- dirname(normalizePath(FINAL_OUT, mustWork = FALSE))
dir.create(parent, recursive = TRUE, showWarnings = FALSE)
STAGE <- tempfile(pattern = ".stage_enrichment.", tmpdir = parent)
dir.create(file.path(STAGE, "go"), recursive = TRUE)
options(error = function() { unlink(STAGE, recursive = TRUE); quit(save = "no", status = 1) })
w <- function(x, name) fwrite(x, file.path(STAGE, name), sep = "\t")
w(joined, "sparcs_candidates_with_ifn.tsv")
w(pc_priority, "sparcs_pc_ifn_priority.tsv")
w(lnc_priority, "sparcs_lnc_ifn_priority.tsv")
w(res_all[order(p_adj)], "ifn_motif_fisher_tests.tsv")
w(bp_all, "GO_ORA_BP_sparcs_pc_vs_proteinCodingUniverse.tsv")
if (!is.null(bp_ifn)) w(bp_ifn, "GO_ORA_BP_sparcs_pc_IFNhigh_vs_proteinCodingUniverse.tsv")
w(c7_all, "immune_msidb_candidates_vs_proteinCodingUniverse.tsv")
w(res_pc[order(p_adj)], "TF_motif_enrichment_pc.tsv")
w(res_lnc[order(p_adj)], "TF_motif_depleted_lncRNA.tsv")        # historical name: ALL lncRNA-branch motifs, by p_adj
for (nm in sort(names(tf_go))) w(tf_go[[nm]], file.path("go", nm))
w(plot_df, "figure5B_source_data.tsv")
w(qc, "sparcs_qc.tsv")
w(gene_lists, "enrichment_gene_lists.tsv")
fs <- sort(list.files(STAGE, pattern = "\\.tsv$", recursive = TRUE))
writeLines(paste0(vapply(fs, function(f) digest(file = file.path(STAGE, f), algo = "sha256"), ""), "  ", fs), file.path(STAGE, "checksums.sha256"))
Sys.chmod(list.files(STAGE, full.names = TRUE, recursive = TRUE), "0644")
Sys.chmod(file.path(STAGE, "go"), "0755")
if (dir.exists(FINAL_OUT)) unlink(FINAL_OUT, recursive = TRUE)
if (!file.rename(STAGE, FINAL_OUT)) stop("could not publish ", FINAL_OUT)
Sys.chmod(FINAL_OUT, "0755")
message("[OK] coding_3UTR ", plot_df$n_insertions[1], " LTRs / ", plot_df$n_affected[1], " genes; lncRNA_last_exon ",
        plot_df$n_insertions[2], " / ", plot_df$n_affected[2], "; GO terms ", nrow(bp_all), "; C7 terms ", nrow(c7_all),
        "; lncRNA depleted motifs ", sum(res_lnc$sig == "Depleted"))
message("[OK] wrote ", FINAL_OUT)
