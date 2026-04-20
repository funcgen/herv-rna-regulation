#!/usr/bin/env Rscript

# ============================================================
# SPARCS-like candidates + IFN motif architecture prioritization
# ============================================================
# Inputs:
#   - ../data/ltr_candidates_antisense_last_exon_or_utr3.tsv
#   - ../results/tables/LTR_IFN_STAT1STAT2_IRF_summary.tsv   (Type I/III proxy)
#   - ../results/tables/LTR_IFN_STAT1_summary.tsv           (Type II proxy)
#   - ../data/gencode.v48.primary_assembly.annotation.gtf   (for GO universe)
#
# Outputs (suggested):
#   - ../results/sparcs_like/sparcs_candidates_with_ifn.tsv
#   - ../results/sparcs_like/sparcs_pc_ifn_priority.tsv
#   - ../results/sparcs_like/sparcs_lnc_ifn_priority.tsv
#   - ../results/sparcs_like/GO_ORA_*.tsv
#
# Notes:
#   - Type I/III IFN proxy: STAT1::STAT2 (ISGF3-like) + IRF co-localization
#   - Type II IFN proxy   : STAT1 (GAF/GAS-like) + IRF co-localization
#   - We rely on the summary tables' ifn_stimulation_potential classes.

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(dplyr)
  library(stringr)
  
  # GO ORA (optional section)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# -----------------------
# Config
# -----------------------
IN_CANDIDATES   <- "../../../data/ltr_candidates_antisense_last_exon_or_utr3.tsv"
IN_IFN_ISGF3    <- "../../../results/interferon/LTR_IFN_STAT1STAT2_IRF_summary.tsv"  # alpha/beta proxy
IN_IFN_GAF      <- "../../../results/interferon/LTR_IFN_STAT1_summary.tsv"          # gamma proxy
IN_GTF          <- "../../../data/gencode.v48.primary_assembly.annotation.gtf"

OUTDIR          <- "../../../results/sparcs_like/"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

OUT_JOINED      <- file.path(OUTDIR, "sparcs_candidates_with_ifn.tsv")
OUT_PC_PRIOR    <- file.path(OUTDIR, "sparcs_pc_ifn_priority.tsv")
OUT_LNC_PRIOR   <- file.path(OUTDIR, "sparcs_lnc_ifn_priority.tsv")

OUT_GO_ALL_PC   <- file.path(OUTDIR, "GO_ORA_BP_sparcs_pc_vs_proteinCodingUniverse.tsv")
OUT_GO_IFN_PC   <- file.path(OUTDIR, "GO_ORA_BP_sparcs_pc_IFNhigh_vs_proteinCodingUniverse.tsv")
OUT_C7_ALL_PC   <- file.path(OUTDIR,"immune_msidb_candidates_vs_proteinCodingUniverse.tsv")

# -----------------------
# Helpers
# -----------------------
infer_subfamily <- function(x) str_replace(x, "_merged.*$", "")

# A simple, explicit tiering function for your ifn_stimulation_potential strings
# (works for both summary tables since they share labels)
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

# Keep only the "load-bearing" IFN columns to avoid huge joins
keep_ifn_cols <- function(dt, prefix) {
  # dt is a data.table read from summary
  # prefix is e.g. "ifn_ab" or "ifn_g"
  cols <- c(
    "sequence_name", "subfamily",
    "n_irf", "has_irf",
    "max_irf_in_75bp",
    "ifn_stimulation_potential",
    "best_ifn_win_start", "best_ifn_win_end"
  )
  
  # STAT columns differ between tables; detect them
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
  
  # rename (except sequence_name) with prefix for clarity after join
  to_rename <- setdiff(names(out), "sequence_name")
  setnames(out, to_rename, paste0(prefix, "__", to_rename))
  out
}

# -----------------------
# Load candidate LTRs (SPARCS-like context table)
# -----------------------
cand <- read_tsv(IN_CANDIDATES, show_col_types = FALSE) %>% as.data.table()

# Standardize
cand[, subfamily := infer_subfamily(a_name)]

# Candidate groups
cand_pc  <- cand[gene_type == "protein_coding" & feature_class == "three_prime_UTR"]
cand_lnc <- cand[gene_type == "lncRNA"]

# Basic counts (printed to console for notes)
message("[INFO] Candidate rows: ", nrow(cand))
message("[INFO] Candidate PC 3'UTR rows: ", nrow(cand_pc))
message("[INFO] Candidate lncRNA rows: ", nrow(cand_lnc))
message("[INFO] Unique PC genes:  ", length(unique(cand_pc$gene_name)))
message("[INFO] Unique lnc genes: ", length(unique(cand_lnc$gene_name)))
message("[INFO] Unique candidate LTRs (a_name): ", length(unique(cand$a_name)))

# -----------------------
# Load IFN motif architecture summaries
# -----------------------
ifn_ab <- fread(IN_IFN_ISGF3)  # STAT1::STAT2 + IRF
ifn_g  <- fread(IN_IFN_GAF)    # STAT1 + IRF

# Reduce & prefix columns to avoid collisions
ifn_ab_small <- keep_ifn_cols(ifn_ab, "ifn_ab")
ifn_g_small  <- keep_ifn_cols(ifn_g,  "ifn_g")

# -----------------------
# Join IFN summaries onto candidates by LTR id
# -----------------------
# Candidate LTR id is a_name, IFN tables use sequence_name
setnames(cand, "a_name", "sequence_name_tmp")
joined <- merge(
  cand,
  ifn_ab_small,
  by.x = "sequence_name_tmp", by.y = "sequence_name",
  all.x = TRUE
)
joined <- merge(
  joined,
  ifn_g_small,
  by.x = "sequence_name_tmp", by.y = "sequence_name",
  all.x = TRUE
)

# restore name
setnames(joined, "sequence_name_tmp", "a_name")

# -----------------------
# Derive IFN tier labels for alpha/beta and gamma
# -----------------------
joined[, ifn_ab__tier := tier_from_ifn_label(ifn_ab__ifn_stimulation_potential)]
joined[, ifn_g__tier  := tier_from_ifn_label(ifn_g__ifn_stimulation_potential)]

# IFN-positive = tiers 1–4 (VeryHigh/High/Medium). Excludes tier 5 (Low) and missing/None.
joined[, ifn_ab__pos := ifn_ab__tier %in% c("VeryHigh","High","Medium")]
joined[, ifn_g__pos  := ifn_g__tier  %in% c("VeryHigh","High","Medium")]

# Combined call (either/both), using the relaxed definition (exclude only tier 5)
joined[, ifn_any_tier := fcase(
  ifn_ab__pos & ifn_g__pos, "Both (Type I/III + Type II)",
  ifn_ab__pos,             "Type I/III-like (STAT1::STAT2)",
  ifn_g__pos,              "Type II-like (STAT1)",
  default = "Not IFN-high"
)]

# -----------------------
# Synergy flag (still strict: only tiers 1–2)
# -----------------------
joined[, ifn_ab__synergy := !is.na(ifn_ab__ifn_stimulation_potential) &
         ifn_ab__ifn_stimulation_potential %chin% c(
           "1_Very High (Synergistic Cluster)",
           "2_High (Synergistic Pair)"
         )]

joined[, ifn_g__synergy := !is.na(ifn_g__ifn_stimulation_potential) &
         ifn_g__ifn_stimulation_potential %chin% c(
           "1_Very High (Synergistic Cluster)",
           "2_High (Synergistic Pair)"
         )]

joined[, ifn_any_synergy := ifn_ab__synergy | ifn_g__synergy]

# -----------------------
# Export joined table
# -----------------------
fwrite(joined, OUT_JOINED, sep = "\t")
message("[OK] Wrote: ", OUT_JOINED)

# -----------------------
# Priority tables for notes
# -----------------------
# Choose a few useful columns to keep tables readable
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

# Protein-coding 3'UTR: prioritize IFN-high first
pc_priority <- joined[
  gene_type == "protein_coding" & feature_class == "three_prime_UTR"
][, ..cols_keep]

# create numeric sort keys safely (handles missing columns too)
pc_priority[, sort_ab_density := suppressWarnings(as.numeric(ifn_ab__stat1stat2_sites_per_kb))]
pc_priority[, sort_g_nstat    := suppressWarnings(as.numeric(ifn_g__n_stat1))]

# NAs -> -Inf so they go to the bottom when sorting decreasing
pc_priority[is.na(sort_ab_density), sort_ab_density := -Inf]
pc_priority[is.na(sort_g_nstat),    sort_g_nstat    := -Inf]

setorder(
  pc_priority,
  ifn_any_tier,
  -ifn_any_synergy,
  -sort_ab_density,
  -sort_g_nstat,
  gene_name
)

pc_priority[, c("sort_ab_density","sort_g_nstat") := NULL]

fwrite(pc_priority, OUT_PC_PRIOR, sep = "\t")
message("[OK] Wrote: ", OUT_PC_PRIOR)

# lncRNA: same idea
lnc_priority <- joined[
  gene_type == "lncRNA"
][, ..cols_keep]

lnc_priority[, sort_ab_density := suppressWarnings(as.numeric(ifn_ab__stat1stat2_sites_per_kb))]
lnc_priority[, sort_g_nstat    := suppressWarnings(as.numeric(ifn_g__n_stat1))]

lnc_priority[is.na(sort_ab_density), sort_ab_density := -Inf]
lnc_priority[is.na(sort_g_nstat),    sort_g_nstat    := -Inf]

setorder(
  lnc_priority,
  ifn_any_tier,
  -ifn_any_synergy,
  -sort_ab_density,
  -sort_g_nstat,
  gene_name
)

lnc_priority[, c("sort_ab_density","sort_g_nstat") := NULL]

fwrite(lnc_priority, OUT_LNC_PRIOR, sep = "\t")
message("[OK] Wrote: ", OUT_LNC_PRIOR)

# -----------------------
# Quick console summaries (good for notes)
# -----------------------
message("\n=== IFN tier counts among SPARCS-like candidates ===")
print(joined[, .N, by = ifn_any_tier][order(-N)])

message("\n=== Type I/III-like tiers (STAT1::STAT2) ===")
print(joined[, .N, by = ifn_ab__tier][order(-N)])

message("\n=== Type II-like tiers (STAT1) ===")
print(joined[, .N, by = ifn_g__tier][order(-N)])

message("\n=== Synergy (STAT+IRF same-window) ===")
print(joined[, .N, by = ifn_any_synergy][order(-N)])

# -----------------------
# GO ORA (Protein-coding SPARCS-like genes)
# -----------------------
# Universe: all protein-coding genes in the GTF (gene_name)
# Candidates: protein-coding SPARCS-like candidate genes (unique gene_name)
# Also run an ORA on the IFN-high subset (optional but informative)

# ---- read 'gene' rows from GTF
gtf_genes <- fread(
  IN_GTF,
  sep = "\t",
  header = FALSE,
  quote = "",
  col.names = c("seqname","source","feature","start","end","score","strand","frame","attribute")
)[feature == "gene"]

get_attr <- function(x, key) {
  m <- str_match(x, paste0(key, ' "([^"]+)"'))
  m[,2]
}

gtf_genes[, gene_type := get_attr(attribute, "gene_type")]
gtf_genes[, gene_name := get_attr(attribute, "gene_name")]

universe_symbols <- unique(na.omit(gtf_genes[gene_type == "protein_coding", gene_name]))
message("\n[INFO] Universe protein_coding symbols: ", length(universe_symbols))

candidate_pc_genes <- unique(na.omit(joined[
  gene_type == "protein_coding" & feature_class == "three_prime_UTR",
  gene_name
]))
candidate_pc_genes <- intersect(candidate_pc_genes, universe_symbols)

# IFN-high PC subset (either tier VeryHigh/High)
candidate_pc_ifnhigh <- unique(na.omit(joined[
  gene_type == "protein_coding" & feature_class == "three_prime_UTR" &
    ifn_any_tier != "Not IFN-high",
  gene_name
]))
candidate_pc_ifnhigh <- intersect(candidate_pc_ifnhigh, universe_symbols)

message("[INFO] PC SPARCS-like genes: ", length(candidate_pc_genes))
message("[INFO] PC SPARCS-like IFN-high genes: ", length(candidate_pc_ifnhigh))

# ---- map SYMBOL -> ENTREZ
cand_map <- bitr(candidate_pc_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
univ_map <- bitr(universe_symbols,   fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

cand_entrez <- unique(cand_map$ENTREZID)
univ_entrez <- unique(univ_map$ENTREZID)

# ORA: all PC candidates
ego_bp_all <- enrichGO(
  gene          = cand_entrez,
  universe      = univ_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20,
  readable      = TRUE
)

bp_all <- as.data.frame(ego_bp_all)
fwrite(bp_all, OUT_GO_ALL_PC, sep="\t")
message("[OK] Wrote: ", OUT_GO_ALL_PC)

# ORA: IFN-high subset (only if enough genes)
if (length(candidate_pc_ifnhigh) >= 10) {
  cand_ifn_map <- bitr(candidate_pc_ifnhigh, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  cand_ifn_entrez <- unique(cand_ifn_map$ENTREZID)
  
  ego_bp_ifn <- enrichGO(
    gene          = cand_ifn_entrez,
    universe      = univ_entrez,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20,
    readable      = TRUE
  )
  
  bp_ifn <- as.data.frame(ego_bp_ifn)
  fwrite(bp_ifn, OUT_GO_IFN_PC, sep="\t")
  message("[OK] Wrote: ", OUT_GO_IFN_PC)
} else {
  message("[INFO] Skipping GO ORA for IFN-high subset (too few genes: ", length(candidate_pc_ifnhigh), ")")
}




# -----------------------
# MSigDB C7 ORA (Protein-coding SPARCS-like genes)
# -----------------------
# Universe: all protein-coding genes in the GTF (gene_name)
# Candidates: protein-coding SPARCS-like candidate genes (unique gene_name)
library(msigdbr)
# ---- MSigDB C7 gene sets
# C7 contains IMMUNESIGDB and VAX subcollections.
# Usually IMMUNESIGDB is the main one people use for immune signatures.
msig_c7 <- msigdbr(species = "Homo sapiens", category = "C7")

# TERM2GENE and TERM2NAME for enricher()
term2gene <- msig_c7 %>%
  dplyr::select(gs_name, gene_symbol) %>%
  distinct()

term2name <- msig_c7 %>%
  dplyr::select(gs_name, gs_description) %>%
  distinct()

# ---- restrict universe/candidates to genes represented in C7
msig_symbols <- unique(term2gene$gene_symbol)

# ---- ORA: all PC candidates against MSigDB C7
ego_c7_all <- enricher(
  gene          = cand_map$SYMBOL,
  universe      = univ_map$SYMBOL,
  TERM2GENE     = term2gene,
  TERM2NAME     = term2name,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.20
)

c7_all <- as.data.frame(ego_c7_all)
fwrite(c7_all, OUT_C7_ALL_PC, sep = "\t")
message("[OK] Wrote: ", OUT_C7_ALL_PC)









suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

# -----------------------
# Inputs
# -----------------------
fimo_all_file <- "../../../data/merged_fimo_parsed_v4.tsv"
fimo_all      <- fread(fimo_all_file)

univ_ltrs <- unique(fimo_all$sequence_name)

# Foreground: SPARCS-like PC candidates (LTR IDs)
sparcs_ltrs <- unique(pc_priority$a_name)

# -----------------------
# Helpers: tier parsing
# -----------------------

# extract tier number 1..5 from your label strings
tier_num <- function(x) {
  # e.g. "1_Very High (...)" -> 1
  # returns NA if missing/doesn't match
  as.integer(str_match(x, "^\\s*([1-5])\\s*[_\\- ]")[,2])
}

is_synergy <- function(x) {
  # tiers 1-2 in your scheme correspond to STAT+IRF in same window
  str_detect(x, regex("^\\s*[12]\\s*[_\\- ]", ignore_case = TRUE))
}

# Define which LTRs are "positive" under different rules
pos_ltrs_from_summary <- function(dt, mode = c("non_low", "synergy_only", "tier1_only", "tier1_2_3", "tier1_2_3_4")) {
  mode <- match.arg(mode)
  t <- dt$ifn_stimulation_potential
  tn <- tier_num(t)
  
  keep <- switch(
    mode,
    non_low     = !is.na(tn) & tn %in% 1:4,          # your current "exclude only tier 5"
    synergy_only = is_synergy(t),                    # tiers 1-2 only
    tier1_only  = !is.na(tn) & tn == 1,
    tier1_2_3   = !is.na(tn) & tn %in% 1:3,
    tier1_2_3_4 = !is.na(tn) & tn %in% 1:4           # same as non_low (kept for readability)
  )
  
  unique(dt$sequence_name[keep])
}

# Fisher enrichment: FG vs (Universe \ FG)
fisher_enrich_ltrs <- function(fg_ltrs, pos_ltrs, universe_ltrs, label) {
  fg_ltrs <- intersect(unique(fg_ltrs), universe_ltrs)
  other_ltrs <- setdiff(universe_ltrs, fg_ltrs)
  
  pos_ltrs <- intersect(unique(pos_ltrs), universe_ltrs)
  
  a <- length(intersect(fg_ltrs, pos_ltrs))   # FG with
  b <- length(fg_ltrs) - a                    # FG without
  c <- length(intersect(other_ltrs, pos_ltrs))# Other with
  d <- length(other_ltrs) - c                 # Other without
  
  mat <- matrix(c(a,b,c,d), nrow=2, byrow=TRUE,
                dimnames=list(c("FG","Other"), c("Pos","Neg")))
  
  ft <- fisher.test(mat)
  
  data.table(
    test = label,
    fg_total = length(fg_ltrs),
    fg_with  = a,
    other_total = length(other_ltrs),
    other_with  = c,
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value
  )
}

# -----------------------
# Run multiple tier modes (keeps your old + adds new)
# -----------------------
modes <- c("non_low","synergy_only","tier1_only","tier1_2_3","tier1_2_3_4")  

res_all <- rbindlist(lapply(modes, function(mode) {
  
  pos_stat1    <- pos_ltrs_from_summary(ifn_g,  mode = mode)
  pos_stat1st2 <- pos_ltrs_from_summary(ifn_ab, mode = mode)
  
  r1 <- fisher_enrich_ltrs(
    fg_ltrs = sparcs_ltrs,
    pos_ltrs = pos_stat1,
    universe_ltrs = univ_ltrs,
    label = paste0("STAT1 | mode=", mode)
  )
  
  r2 <- fisher_enrich_ltrs(
    fg_ltrs = sparcs_ltrs,
    pos_ltrs = pos_stat1st2,
    universe_ltrs = univ_ltrs,
    label = paste0("STAT1::STAT2 | mode=", mode)
  )
  
  rbindlist(list(r1, r2))
}))

res_all[, p_adj := p.adjust(p_value, method = "BH")]
print(res_all[order(p_adj)])




#-------------
# Motif enrichment analysis in SPARCS-like LTRs
# -----------

# Define foregrounds
fg_pc  <- unique(cand_pc$a_name)
fg_lnc <- unique(cand_lnc$a_name)

univ_ltrs <- unique(fimo_all$sequence_name)

# Background
bg_pc  <- setdiff(univ_ltrs, fg_pc)
bg_lnc <- setdiff(univ_ltrs, fg_lnc)

# Deduplicate motif presence
fimo_hits <- unique(
  fimo_all[, .(sequence_name, motif_alt_id)]
)



# Fisher enrichment function
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
    
    data.table(
      motif_alt_id = m,
      fg_with = a,
      fg_total = length(fg_ltrs),
      bg_with = c,
      bg_total = length(bg_ltrs),
      odds_ratio = unname(ft$estimate),
      p_value = ft$p.value
    )
    
  })
  
  res <- rbindlist(res)
  res[, p_adj := p.adjust(p_value, method="BH")]
  res[, log2_OR := log2(pmax(odds_ratio, .Machine$double.eps))]
  
  res[]
}

# Run enrichment
res_pc  <- motif_enrichment(fg_pc, bg_pc, fimo_hits)
res_lnc <- motif_enrichment(fg_lnc, bg_lnc, fimo_hits)

# Add interpretation labels
add_labels <- function(dt) {
  dt[, sig := fcase(
    p_adj < 0.05 & odds_ratio > 1, "Enriched",
    p_adj < 0.05 & odds_ratio < 1, "Depleted",
    default = "NS"
  )]
}

res_pc  <- add_labels(res_pc)
res_lnc <- add_labels(res_lnc)

res_pc


# ------------------------------------------------------------------------------
# GO over enriched/depleted TF motif IDs from SPARCS-like motif ORA
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
})

go_dir <- file.path(OUTDIR, "go")
dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)

# ---- helpers: split complex IDs, normalize, and map to human SYMBOLs ----
explode_and_normalize <- function(x) {
  if (is.null(x) || length(x) == 0) return(character(0))
  x <- as.character(x)
  x <- x[!is.na(x)]
  if (!length(x)) return(character(0))
  
  # Split complexes like "MAX::MYC" -> c("MAX","MYC")
  parts <- unlist(strsplit(x, "::", fixed = TRUE), use.names = FALSE)
  parts <- trimws(parts)
  parts <- toupper(parts)
  
  # Remove trailing notes if any
  parts <- sub("\\s+\\(.*\\)$", "", parts)
  
  # Keep only letters, numbers, underscore, hyphen
  parts <- gsub("[^A-Z0-9_\\-]", "", parts)
  parts <- parts[nchar(parts) > 0]
  
  unique(parts)
}

map_to_hgnc <- function(ids_upper) {
  ids_upper <- unique(ids_upper)
  if (length(ids_upper) == 0) return(character(0))
  
  # 1) Direct SYMBOL mapping
  sym_ok <- tryCatch(
    AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = ids_upper,
      keytype = "SYMBOL",
      column = "SYMBOL",
      multiVals = "first"
    ),
    error = function(e) setNames(rep(NA_character_, length(ids_upper)), ids_upper)
  )
  
  # 2) Fallback through ALIAS
  need_alias <- names(sym_ok)[is.na(sym_ok)]
  alias_map <- if (length(need_alias)) {
    tryCatch(
      AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = need_alias,
        keytype = "ALIAS",
        column = "SYMBOL",
        multiVals = "first"
      ),
      error = function(e) setNames(rep(NA_character_, length(need_alias)), need_alias)
    )
  } else {
    setNames(character(0), character(0))
  }
  
  out <- sym_ok
  out[is.na(out)] <- alias_map[names(out)[is.na(out)]]
  out <- unique(stats::na.omit(unname(out)))
  out
}

# ---- reusable GO function ----
run_tf_go <- function(res_dt,
                      foreground_sig = c("Enriched", "Depleted"),
                      label = "SPARCS_like",
                      out_prefix = "sparcs") {
  
  foreground_sig <- match.arg(foreground_sig)
  
  stopifnot(is.data.table(res_dt) || is.data.frame(res_dt))
  res_dt <- as.data.table(res_dt)
  
  if (!all(c("motif_alt_id", "sig") %in% names(res_dt))) {
    stop("res_dt must contain columns: motif_alt_id, sig")
  }
  
  # Foreground = enriched or depleted motifs
  fg_raw <- unique(res_dt[sig == foreground_sig, motif_alt_id])
  
  # Universe = all tested motifs
  bg_raw <- unique(res_dt$motif_alt_id)
  
  fg_norm <- explode_and_normalize(fg_raw)
  bg_norm <- explode_and_normalize(bg_raw)
  
  tf_symbols  <- map_to_hgnc(fg_norm)
  tf_universe <- map_to_hgnc(bg_norm)
  
  message(sprintf(
    "[INFO] %s | %s TF GO: foreground raw=%d -> norm=%d -> mapped=%d; universe raw=%d -> norm=%d -> mapped=%d",
    label, foreground_sig,
    length(unique(fg_raw)), length(fg_norm), length(tf_symbols),
    length(unique(bg_raw)), length(bg_norm), length(tf_universe)
  ))
  
  # Save audit tables
  fwrite(
    data.table(raw_motif_id = unique(fg_raw)),
    file.path(go_dir, paste0(out_prefix, "_", tolower(foreground_sig), "_TF_raw_ids.tsv")),
    sep = "\t"
  )
  
  fwrite(
    data.table(normalized_motif_id = unique(fg_norm)),
    file.path(go_dir, paste0(out_prefix, "_", tolower(foreground_sig), "_TF_normalized_ids.tsv")),
    sep = "\t"
  )
  
  fwrite(
    data.table(mapped_symbol = unique(tf_symbols)),
    file.path(go_dir, paste0(out_prefix, "_", tolower(foreground_sig), "_TF_mapped_symbols.tsv")),
    sep = "\t"
  )
  
  # Unmapped audit
  unmapped_fg <- setdiff(unique(fg_norm), unique(tf_symbols))
  if (length(unmapped_fg)) {
    fwrite(
      data.table(unmapped_motif_id = unmapped_fg),
      file.path(go_dir, paste0(out_prefix, "_", tolower(foreground_sig), "_TF_unmapped_ids.tsv")),
      sep = "\t"
    )
  }
  
  # Skip gracefully if empty
  if (length(tf_symbols) == 0L) {
    message(sprintf("[INFO] %s | %s TF GO: nothing to test — no mapped TF symbols.", label, foreground_sig))
    return(NULL)
  }
  
  if (length(tf_universe) == 0L) {
    message(sprintf("[INFO] %s | %s TF GO: empty universe after mapping.", label, foreground_sig))
    return(NULL)
  }
  
  ego_tf <- tryCatch(
    enrichGO(
      gene           = tf_symbols,
      universe       = tf_universe,
      OrgDb          = org.Hs.eg.db,
      keyType        = "SYMBOL",
      ont            = "BP",
      pAdjustMethod  = "BH",
      pvalueCutoff   = 0.1,
      qvalueCutoff   = 0.2,
      readable       = TRUE
    ),
    error = function(e) {
      message(sprintf("[INFO] %s | %s TF GO error: %s", label, foreground_sig, conditionMessage(e)))
      NULL
    }
  )
  
  if (is.null(ego_tf) || NROW(as.data.frame(ego_tf)) == 0) {
    message(sprintf("[INFO] %s | %s TF GO: no significant BP terms.", label, foreground_sig))
    return(NULL)
  }
  
  ego_df <- as.data.table(as.data.frame(ego_tf))
  out_tsv <- file.path(go_dir, paste0(out_prefix, "_GO_BP_over_", tolower(foreground_sig), "_TFsymbols.tsv"))
  fwrite(ego_df, out_tsv, sep = "\t")
  
  p_tf <- enrichplot::dotplot(ego_tf, showCategory = 20, font.size = 12) +
    ggtitle(sprintf("GO: BP — %s TF motifs (%s)", foreground_sig, label))
  
  out_png <- file.path(go_dir, paste0(out_prefix, "_GO_BP_over_", tolower(foreground_sig), "_TFsymbols.png"))
  ggsave(out_png, p_tf, width = 10, height = 7, dpi = 300, bg = "white")
  
  message("[OK] Wrote: ", out_tsv)
  message("[OK] Wrote: ", out_png)
  
  invisible(list(
    ego = ego_tf,
    table = ego_df,
    plot = p_tf,
    tf_symbols = tf_symbols,
    tf_universe = tf_universe
  ))
}


go_pc_enriched <- run_tf_go(
  res_dt = res_pc,
  foreground_sig = "Enriched",
  label = "SPARCS-like protein-coding 3UTR LTRs",
  out_prefix = "sparcs_pc_enriched"
)

go_pc_depleted <- run_tf_go(
  res_dt = res_pc,
  foreground_sig = "Depleted",
  label = "SPARCS-like protein-coding 3UTR LTRs",
  out_prefix = "sparcs_pc_depleted"
)

go_lnc_enriched <- run_tf_go(
  res_dt = res_lnc,
  foreground_sig = "Enriched",
  label = "SPARCS-like lncRNA LTRs",
  out_prefix = "sparcs_lnc_enriched"
)

go_lnc_depleted <- run_tf_go(
  res_dt = res_lnc,
  foreground_sig = "Depleted",
  label = "SPARCS-like lncRNA LTRs",
  out_prefix = "sparcs_lnc_depleted"
)





# ============================================================
# Panel: Global counts / composition of SPARCS-like candidates
# Built from final candidate subsets:
#   - cand_pc  : protein-coding, 3'UTR
#   - cand_lnc : lncRNA
# Counts use:
#   - unique a_name  = unique LTR insertions
#   - unique gene_id = affected genes / loci
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# -----------------------
# Output paths
# -----------------------
OUT_PANEL_COUNTS_PDF <- file.path(OUTDIR, "panel_global_counts_composition.pdf")
OUT_PANEL_COUNTS_PNG <- file.path(OUTDIR, "panel_global_counts_composition.png")
OUT_COUNTS_TSV       <- file.path(OUTDIR, "panel_global_counts_composition_source_data.tsv")

# -----------------------
# Build source table
# -----------------------
plot_df <- data.table(
  class = factor(
    c("Protein-coding\n(3' UTR)", "lncRNA\n(last exon)"),
    levels = c("Protein-coding\n(3' UTR)", "lncRNA\n(last exon)")
  ),
  n_insertions = c(
    uniqueN(cand_pc$a_name),   # 1912
    uniqueN(cand_lnc$a_name)   # 4682
  ),
  n_affected = c(
    uniqueN(cand_pc$gene_id),  # 1398
    uniqueN(cand_lnc$gene_id)  # 3400
  )
)

# percentages for optional labels / notes
plot_df[, insertion_pct := n_insertions / sum(n_insertions)]
plot_df[, affected_pct  := n_affected   / sum(n_affected)]

# save source data
fwrite(plot_df, OUT_COUNTS_TSV, sep = "\t")

# -----------------------
# Colors
# -----------------------
class_cols <- c(
  "Protein-coding\n(3' UTR)" = "#4C78A8",
  "lncRNA\n(last exon)"      = "#54A24B"
)

# -----------------------
# Left: unique insertions
# -----------------------
p_insertions <- ggplot(
  plot_df,
  aes(x = class, y = n_insertions, fill = class)
) +
  geom_col(width = 0.68, color = "black", linewidth = 0.25, show.legend = FALSE) +
  geom_text(
    aes(label = comma(n_insertions)),
    vjust = -0.45,
    size = 4.4
  ) +
  scale_fill_manual(values = class_cols) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Unique LTR insertions",
    x = NULL,
    y = "Count"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 11)
  )

# -----------------------
# Right: affected genes / loci
# -----------------------
p_affected <- ggplot(
  plot_df,
  aes(x = class, y = n_affected, fill = class)
) +
  geom_col(width = 0.68, color = "black", linewidth = 0.25, show.legend = FALSE) +
  geom_text(
    aes(label = comma(n_affected)),
    vjust = -0.45,
    size = 4.4
  ) +
  scale_fill_manual(values = class_cols) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Affected genes / loci",
    x = NULL,
    y = "Count"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 11)
  )

# -----------------------
# Combine
# -----------------------
p_counts_panel <- p_insertions + p_affected +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "Global composition of SPARCS-like transcript architectures",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
    )
  )

# save
ggsave(OUT_PANEL_COUNTS_PDF, p_counts_panel, width = 8.5, height = 4.8, bg = "white")
ggsave(OUT_PANEL_COUNTS_PNG, p_counts_panel, width = 8.5, height = 4.8, dpi = 300, bg = "white")

message("[OK] Wrote: ", OUT_PANEL_COUNTS_PDF)
message("[OK] Wrote: ", OUT_PANEL_COUNTS_PNG)
message("[OK] Wrote: ", OUT_COUNTS_TSV)

# quick console check
print(plot_df)
