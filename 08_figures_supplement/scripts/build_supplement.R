#!/usr/bin/env Rscript
# 08_figures_supplement step 8 (frozen and effective order 54; extended by effective orders 54.2, erratum D038, 54.3, erratum
# D041, and 54.4, erratum D042): build the supplementary workbooks and tables.
#
# PUBLICATION / EXPORT ONLY. Every value comes from an accepted upstream result table of 02_rbp, 05_lncrna, 06_lncrna_function or
# 07_sparcs, read as text. This script selects accepted columns, selects rows by accepted labels, adds constant label columns,
# orders sheets, sets Excel cell types and display formats and writes XLSX / TSV / one tar.gz. It runs no statistical test,
# multiple-testing adjustment, enrichment, motif scan, classification or coverage calculation, and it does not change any
# threshold or universe. The one assembled sheet (lncRNA functional annotation, order 54.3) joins accepted 06 tables by gene and
# only re-presents accepted values (see derive_functional_annotation).
#
# Workbooks are named by scientific content (content_id), never by number. The publication labels and the labels cited by the
# 2026-09-22 manuscript snapshot are metadata only (validation/supplement/supplement_numbering_decision.tsv); the manuscript text
# is updated separately. Supplemental Tables 1-6 were APPROVED_BY_AUTHOR on 2026-10-01 (order-54 adjudication); on 2026-10-06
# effective order 54.2 (D038) added lncrna_herv_coverage as Supplemental Table 4 and renumbered the protein-coding, C7 and TF
# workbooks 4-6 -> 5-7 (APPROVED_BY_AUTHOR); on 2026-10-11 effective order 54.3 (D041) added lncrna_functional_annotation as
# Supplemental Table 6 and renumbered the C7 and TF workbooks 6-7 -> 7-8 (APPROVED_BY_AUTHOR); on 2026-10-11 effective order 54.4
# (D042) replaced the C7 workbook by sparcs_like_catalogue (Supplemental Table 7: the SPARCS-like catalogue, Figure 5C GO, the
# unchanged C7 sheet, interferon-motif and U3-R-U5 tests). Numbering and release dates are per workbook.
#
# Published workbooks hold only their data sheets and a final Column_definitions sheet (Sheet / Column / Definition, from
# scripts/supplement/column_dictionary.tsv): no title, number, date, note or path (author decision 2026-10-08, order 54.2). The
# publication metadata (label, dates, title, notes, display rule) is written to the release record
# panel_sources/supplement/supplement_workbook_notes.tsv instead.
#
# Display rules (recorded in supplement_workbook_notes.tsv):
#   - numeric cells hold the accepted source token verbatim (e.g. 3.2787860589288393e-187), so each Excel cell is the accepted
#     double exactly; openxlsx's own 15-significant-digit fixed-point rendering is replaced in the sheet XML (normalise_xlsx)
#   - p-value-like cells are numeric with an "0.00E+00" display format; blank = NA in the accepted source (no valid test)
#   - an exact 0 in a p-value-like source column STOPS the build unless that cell is classified in
#     scripts/supplement/zero_p_value_adjudications.tsv (empty: the accepted canonical sources contain no such zero; the order-54
#     audit is validation/supplement/p_value_zero_audit.md). An adjudicated zero is shown with its adjudicated text and listed in
#     panel_sources/supplement/underflow_display_cells.tsv; the TSV support file keeps the source "0"
#   - identifiers (GO ids, Ensembl ids, gene / motif / gene-set names, ratios such as "32/1315") are written as text
#
# Usage: build_supplement.R --release <ncRNA root> --out <staging root> --work <scratch dir>
suppressPackageStartupMessages({
  library(openxlsx)
  library(data.table)
})
Sys.setenv(TZ = "UTC")  # ZIP member times are local-time DOS fields: fix the zone so the bytes do not depend on the host

args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("missing ", flag)
  args[[i + 1]]
}
REL <- normalizePath(arg("--release"), mustWork = TRUE)
OUT <- arg("--out")
WORK <- arg("--work")
M <- file.path(REL, "08_figures_supplement")
RELEASE_DATE <- "2026-10-01"                                  # container date (core.xml created) of every workbook
FIXED_TIME <- as.POSIXct("2026-10-01 00:00:00", tz = "UTC")    # every ZIP member time; a build constant, not a release claim
NUMERIC_TOKEN <- "^-?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?$"  # the only source spellings accepted as numeric cells (no Inf/NaN)
MISSING_TOKEN <- "NA"                                          # the upstream missing-value spelling; blank cell in the XLSX
# an empty field in a numeric column is also missing (blank cell): the order-46 U3-R-U5 test table leaves the Fisher-only columns
# empty in its Wilcoxon rows (order 54.4). No source consumed before order 54.4 has an empty numeric field (the build stopped on it)
NUMBERING_STATUS <- "APPROVED_BY_AUTHOR"                       # order-54 and order-54.2 (D038) adjudications of the labels

dir.create(file.path(OUT, "supplement/tables/by_content"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "supplement/tables/by_content_tsv"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "panel_sources/supplement"), recursive = TRUE, showWarnings = FALSE)
dir.create(WORK, recursive = TRUE, showWarnings = FALSE)

dict <- fread(file.path(M, "scripts/supplement/column_dictionary.tsv"), sep = "\t", colClasses = "character", na.strings = NULL)
zero_adj <- fread(file.path(M, "scripts/supplement/zero_p_value_adjudications.tsv"), sep = "\t", colClasses = "character",
                  na.strings = NULL)

# ------------------------------------------------------------------------------------------------------------- source registry
R2 <- "02_rbp/results"
SRC <- list(
  sub_all     = c("02", "23", file.path(R2, "enrichment/subfamily/rbp_subfamily_composition_controlled_enrichment.tsv")),
  sub_hk      = c("02", "23", file.path(R2, "enrichment/subfamily/rbp_HERVH_HERVK_composition_controlled_enrichment.tsv")),
  sub_sig     = c("02", "23", file.path(R2, "enrichment/subfamily/rbp_HERVH_HERVK_significant.tsv")),
  gsea        = c("02", "24", file.path(R2, "enrichment/hervh_hervk_followup/GSEA_HERVH_vs_HERVK_GO_BP_complete.tsv")),
  gsea_rank   = c("02", "24", file.path(R2, "enrichment/hervh_hervk_followup/GSEA_HERVH_vs_HERVK_RBP_ranking.tsv")),
  t_global    = c("02", "27", file.path(R2, "sensitivity_strict/analysis/threshold_comparison/global_threshold_comparison.tsv")),
  t_family    = c("02", "27", file.path(R2, "sensitivity_strict/analysis/threshold_comparison/HERVH_HERVK_global_threshold_comparison.tsv")),
  t_sets      = c("02", "27", file.path(R2, "sensitivity_strict/analysis/threshold_comparison/HERVH_HERVK_significant_set_summary.tsv")),
  t_sparsity  = c("02", "27", file.path(R2, "sensitivity_strict/analysis/threshold_comparison/strict_threshold_sparsity.tsv")),
  t_per_rbp   = c("02", "27", file.path(R2, "sensitivity_strict/analysis/threshold_comparison/HERVH_HERVK_per_RBP_threshold_comparison.tsv")),
  t_func      = c("02", "27", file.path(R2, "sensitivity_strict/analysis/sensitivity_summary/fimo_threshold_functional_summary.tsv")),
  go_h_1e4    = c("02", "24", file.path(R2, "enrichment/hervh_hervk_followup/GO_HERVH_enriched_RBPs_BP_significant_FDR05.tsv")),
  go_k_1e4    = c("02", "24", file.path(R2, "enrichment/hervh_hervk_followup/GO_HERVK_enriched_RBPs_BP_significant_FDR05.tsv")),
  go_h_1e5    = c("02", "26", file.path(R2, "sensitivity_strict/analysis/hervh_hervk_followup/GO_HERVH_enriched_RBPs_BP_significant_FDR05.tsv")),
  go_k_1e5    = c("02", "26", file.path(R2, "sensitivity_strict/analysis/hervh_hervk_followup/GO_HERVK_enriched_RBPs_BP_significant_FDR05.tsv")),
  gsea_1e4    = c("02", "24", file.path(R2, "enrichment/hervh_hervk_followup/GSEA_HERVH_vs_HERVK_GO_BP_significant_FDR05.tsv")),
  gsea_1e5    = c("02", "26", file.path(R2, "sensitivity_strict/analysis/hervh_hervk_followup/GSEA_HERVH_vs_HERVK_GO_BP_significant_FDR05.tsv")),
  lnc_cov     = c("05", "35", "05_lncrna/results/herv_coverage/lncRNA_herv_coverage.tsv"),
  pc          = c("05", "36", "05_lncrna/results/pc_herv_coverage/protein_coding_overlapping_herv_coverage.tsv"),
  fun_int     = c("06", "42", "06_lncrna_function/results/integrated_annotation/herv_lncrna_with_literature.tsv"),
  fun_mem     = c("06", "42", "06_lncrna_function/results/integrated_annotation/annotation_membership.tsv"),
  fun_cat     = c("06", "42", "06_lncrna_function/results/integrated_annotation/figure4_panelB_gene_categories.tsv"),
  fun_catn    = c("06", "42", "06_lncrna_function/results/integrated_annotation/figure4_panelB_functional_categories.tsv"),
  c7          = c("07", "45", "07_sparcs/results/enrichment/immune_msidb_candidates_vs_proteinCodingUniverse.tsv"),
  sp_branch   = c("07", "43", "07_sparcs/results/candidates/sparcs_candidate_branches.tsv"),
  sp_go       = c("07", "45", "07_sparcs/results/enrichment/GO_ORA_BP_sparcs_pc_vs_proteinCodingUniverse.tsv"),
  sp_ifn      = c("07", "45", "07_sparcs/results/enrichment/ifn_motif_fisher_tests.tsv"),
  sp_u3       = c("07", "46", "07_sparcs/results/u3ru5_integrity/u3ru5_tests.tsv"),
  tf_lnc      = c("07", "45", "07_sparcs/results/enrichment/TF_motif_depleted_lncRNA.tsv"),
  tf_go       = c("07", "45", "07_sparcs/results/enrichment/go/sparcs_lnc_depleted_GO_BP_over_depleted_TFsymbols.tsv")
)

manifest_sha <- function(rel) {
  # The accepted checksum of the source, read from the upstream module's acceptance records (never recomputed into the release).
  # Modules keep these records at different levels (02_rbp: module checksums.sha256; 05_lncrna: results/<dir>/checksums.sha256
  # and RESULTS_MANIFEST.tsv; 07_sparcs: all three), so every checksums.sha256 from the file's directory up to the module root is
  # consulted, plus <module>/RESULTS_MANIFEST.tsv (paths relative to results/). At least one exact-path entry must exist and all
  # entries must agree.
  parts <- strsplit(rel, "/", fixed = TRUE)[[1]]
  mod <- parts[1]
  hits <- character(0)
  for (depth in seq(length(parts) - 1, 1)) {
    man <- file.path(REL, paste(parts[seq_len(depth)], collapse = "/"), "checksums.sha256")
    if (!file.exists(man)) next
    lines <- readLines(man)
    sha <- sub("  .*$", "", lines); path <- sub("^\\./", "", sub("^[0-9a-f]{64}  ", "", lines))
    hits <- c(hits, sha[path == paste(parts[-seq_len(depth)], collapse = "/")])
  }
  rm_path <- file.path(REL, mod, "RESULTS_MANIFEST.tsv")
  if (file.exists(rm_path) && parts[2] == "results") {
    rm <- fread(rm_path, sep = "\t", colClasses = "character", skip = "path\t")
    sha_col <- intersect(c("sha256", "decompressed_sha256"), names(rm))[1]
    hits <- c(hits, rm[[sha_col]][sub("^\\./", "", rm$path) == paste(parts[-(1:2)], collapse = "/") & rm$status == "ACCEPTED"])
  }
  hits <- unique(hits)
  if (length(hits) != 1) stop("source must have exactly one accepted checksum in its module records (found ", length(hits), "): ", rel)
  hits
}
file_sha <- function(path) digest::digest(file = path, algo = "sha256")

read_src <- function(key) {
  rel <- SRC[[key]][3]
  path <- file.path(REL, rel)
  if (file_sha(path) != manifest_sha(rel)) stop("accepted source changed since acceptance: ", rel)
  # line terminator: every line CRLF (05_lncrna protein-coding table, historical byte-compatible) or every line LF; anything else
  # (mixed terminators, a CR inside a field) stops the build. The terminator is not part of any field; support TSVs use LF.
  b <- readBin(path, "raw", file.size(path))
  n_cr <- sum(b == as.raw(13)); n_lf <- sum(b == as.raw(10))
  n_crlf <- if (length(b) > 1) sum(b[-length(b)] == as.raw(13) & b[-1] == as.raw(10)) else 0
  if (!(n_cr == 0 || (n_cr == n_crlf && n_crlf == n_lf))) stop("inconsistent line terminators in accepted source: ", rel)
  # raw strings: na.strings = NULL keeps "NA" as the literal upstream token, so the TSV support files are the source fields exactly
  dt <- fread(path, sep = "\t", colClasses = "character", na.strings = NULL, encoding = "UTF-8", quote = "")
  attr(dt, "src_key") <- key
  attr(dt, "terminator") <- if (n_cr == 0) "LF" else "CRLF"
  dt
}
source_terminator <- function(key) attr(read_src(key), "terminator")

# ------------------------------------------------------------------------------------------------------------- products
#   sheet: name, src (one key, or several keys stacked with constant label columns), keep (columns; NULL = all minus drop),
#          drop, where (column, accepted label), add (named constant columns, per stacked key),
#          derive (the name of an assembly function; its src keys are the accepted tables it joins, recorded in the inventory),
#          scope (dictionary scope of this sheet when it differs from the product's; order 54.4),
#          same_test (c(mode shown twice, mode kept)): rows of the first interferon-tier mode are dropped after checking that each
#          equals the row of the second mode with the same STAT profile in every other column (order 54.4, author decision)
#   numbering_date / release_date: when the publication label was approved / the workbook content was first released
#          (default 2026-10-01, the order-54 adjudication; set only on the workbooks added or renumbered by orders 54.2 and 54.3)
# Columns of accepted sources omitted from the publication tables (author decision 2026-10-08, order 54.2): redundant per-kb
# densities, descriptive-only medians / within-subfamily adjustment, and a locus count with no shuffled counterpart. The accepted
# 02_rbp tables keep them.
DISPLAY_COLS <- c("p_value_display", "p_adj_global_display")
TRIM_SUBFAMILY <- c("observed_hits_per_kb", "mean_shuffled_hits_per_kb", "density_difference_per_kb", "loci_with_observed_hit",
                    "median_observed_hits_per_locus", "median_shuffled_hits_per_locus", "p_adj_within_subfamily")
TRIM_THRESHOLD <- c("observed_hits_per_kb_1e4", "mean_shuffled_hits_per_kb_1e4", "observed_hits_per_kb_1e5", "mean_shuffled_hits_per_kb_1e5")
TRIM_GSEA <- "log2err"                                           # fgsea p-value accuracy estimate (Table 2, author decision 2026-10-09)
TRIM_RANKING <- c("rank_statistic_used", "rank_position")       # tie-broken copy of the contrast; row order (Table 2, 2026-10-09)
# Table 3 (author decision 2026-10-09): 5 of 13 sheets kept and renamed; descriptive columns omitted; source archive no longer published
TRIM_FAMILY <- c("total_length_kb", "observed_hits_per_kb", "mean_shuffled_hits_per_kb", "median_observed_locus_density",
                 "median_shuffled_locus_density", "median_distinct_rbps_per_locus")
TRIM_COMPARISON <- c("observed_total_hits_1e4", "mean_shuffled_total_hits_1e4", "direction_1e4", "direction_1e5", "effect_difference_1e5_minus_1e4")
TRIM_NOTE_T3 <- "omitted from the publication table (author decision 2026-10-09, effective order 54.2): derivable, in Table 1, or descriptive only; kept in the accepted source"
# Tables 4 and 5 (author decision 2026-10-09): gene_type (constant) omitted; percentages first; rows by exact gene-body coverage
COVERAGE_LEAD <- c("gene_stable_id", "gene_name", "gene_body_herv_percent", "exonic_herv_percent")
COVERAGE_DROP <- c("gene_type", "gene_body_herv_fraction", "exonic_herv_fraction")  # constant column; fractions duplicate the percentages
TRIM_NOTE_T7 <- "omitted from the publication table (author decision 2026-10-09, effective order 54.2): constant (totals stated in the column definitions) or redundant (log2 of the odds ratio; direction given by the sheet); kept in the accepted source"
TRIM_NOTE_T45 <- "omitted from the publication table (author decision 2026-10-09, effective order 54.2): constant column (the table holds one gene type) or fraction duplicating a percentage column; kept in the accepted source"
# Table 7 (author decision 2026-10-09): constant totals (stated in the definitions), log2_OR (= log2 odds_ratio) and sig (= sheet) omitted
TRIM_TF <- c("fg_total", "bg_total", "log2_OR", "sig")
TRIM_NOTE <- "omitted from the publication table (author decision 2026-10-08, effective order 54.2): redundant or descriptive only; kept in the accepted source"
# Table 7 (author decision 2026-10-11, order 54.4): columns constant within a catalogue sheet, derivable (ltr_len_bp = a_end - a_start;
# overlap_frac_of_ltr = overlap / length), or not cited; the IFN and U3-R-U5 totals are stated in the column definitions
SPARCS_LEAD <- c("a_name", "a_chrom", "a_start", "a_end", "a_strand", "gene_id", "gene_name")
TRIM_BRANCH <- c("branch", "gene_type", "tx_type", "feature_class", "ltr_len_bp", "overlap_frac_of_ltr")
TRIM_IFN <- c("fg_total", "other_total")
TRIM_U3 <- c("background_group", "ci95_low", "ci95_high", "statistic_W")
U3_LEAD <- c("comparison", "test", "sparcs_TRUE", "sparcs_FALSE", "background_TRUE", "background_FALSE")
TRIM_NOTE_SPC <- "omitted from the publication table (author decision 2026-10-11, effective order 54.4): constant within the sheet, derivable, a constant total stated in the column definitions, or not cited; kept in the accepted source"
TRIM_NOTE_T2 <- "omitted from the publication table (author decision 2026-10-09, effective order 54.2): technical or redundant; kept in the accepted source"
P <- list(
  list(id = "rbp_subfamily_enrichment", scope = "rbp",
       title = "Composition-controlled RBP motif enrichment across HERV subfamilies (primary FIMO threshold P <= 1e-4)",
       current = "Supplemental Table 1", publication = "Supplemental Table 1",
       notes = c("Primary analysis: FIMO P <= 1e-4; 151 motif matrices for 93 RBPs; observed HERV internal regions vs five dinucleotide-preserving shuffled versions of each locus.",
                 "HERVH: 20 enriched / 40 depleted RBPs; HERVK: 15 enriched / 37 depleted (global BH adjusted P < 0.05).",
                 "Predicted motif enrichment does not demonstrate physical RBP binding."),
       sheets = list(
         list(name = "HERVH_HERVK", src = "sub_hk", drop = c(DISPLAY_COLS, TRIM_SUBFAMILY)),
         list(name = "HERVH_HERVK_significant", src = "sub_sig", drop = c(DISPLAY_COLS, TRIM_SUBFAMILY)),
         list(name = "All_subfamilies", src = "sub_all", drop = c(DISPLAY_COLS, TRIM_SUBFAMILY)))),
  list(id = "rbp_gsea", scope = "gsea",
       title = "Pre-ranked GO Biological Process analysis of the HERVH-versus-HERVK composition-controlled RBP contrast",
       current = "Supplemental Table 2 (RBP section; number also used for the protein-coding HERV-coverage table)", publication = "Supplemental Table 2",
       notes = c("All 93 tested RBPs ranked by HERVH - HERVK composition-controlled log2(O/E); fgseaMultilevel (fgsea 1.28.0), GO BP gene sets of 5-50 ranked RBPs (org.Hs.eg.db / GO.db 3.18.0), seed 260730.",
                 "Significant at BH adjusted P < 0.05: GO:0000398 mRNA splicing, via spliceosome, toward the HERVK-biased end.",
                 "Primary FIMO threshold P <= 1e-4."),
       sheets = list(
         list(name = "GSEA_GO_BP", src = "gsea", drop = TRIM_GSEA),
         list(name = "RBP_ranking", src = "gsea_rank", drop = TRIM_RANKING))),
  list(id = "rbp_threshold_sensitivity", scope = "sens",
       title = "Sensitivity of the composition-controlled RBP motif analyses to the FIMO threshold (P <= 1e-5 sensitivity vs P <= 1e-4 primary)",
       current = "Supplemental Table X (placeholder)", publication = "Supplemental Table 3",
       notes = c("PRIMARY threshold = FIMO P <= 1e-4 (the main analysis). SENSITIVITY threshold = FIMO P <= 1e-5, obtained by filtering the complete primary FIMO output; it is not the primary result.",
                 "The same five dinucleotide-preserving shuffled controls, Wilcoxon tests and global BH correction were applied at both thresholds.",
                 "Published sheets (author decision 2026-10-09): Global_burden, Family_burden (accepted HERVH_HERVK global threshold comparison), Significant_RBP_counts (accepted significant-set summary), RBP_threshold_comparison, Functional_summary. Not published (kept in the accepted 02_rbp sources): significant-set members, strict sparsity, status transitions, significant reversals, effect correlation, significant GO / GSEA rows, interpretation.",
                 "Source files of the historical Supplemental Table X archive: panel_sources/supplement/rbp_threshold_sensitivity_source_files.tar.gz (11 members; release record, not published)."),
       sheets = list(
         list(name = "Global_burden", src = "t_global", drop = "observed_minus_null"),
         list(name = "Family_burden", src = "t_family", drop = TRIM_FAMILY),
         list(name = "Significant_RBP_counts", src = "t_sets"),
         list(name = "RBP_threshold_comparison", src = "t_per_rbp", drop = c(TRIM_THRESHOLD, TRIM_COMPARISON)),
         list(name = "Functional_summary", src = "t_func"))),
  list(id = "lncrna_herv_coverage", scope = "lnc",
       title = "HERV internal-region coverage of lncRNA genes overlapping same-strand HERV internal regions",
       current = "not cited (workbook added at effective order 54.2, erratum D038)", publication = "Supplemental Table 4",
       numbering_date = "2026-10-06", release_date = "2026-10-06",
       notes = c("4,074 lncRNA genes (GENCODE v48, GRCh38) overlapping at least one same-strand HERV internal region; rows sorted by exact gene-body HERV coverage (descending), then exact exonic coverage (descending), then gene id.",
                 "Coverage uses the D027-corrected coordinate semantics (HERV BED intervals read as 0-based half-open).",
                 "Gene-body HERV coverage: <= 10% in 3,043 genes, > 10% to 50% in 750, > 50% to 80% in 137, > 80% in 144, with each interval evaluated on the exact ratio herv_bp_in_gene_body / gene_body_bp. The rounded fraction and percentage columns are for display: ENSG00000225964 shows 10.0% but its exact gene-body coverage is 10.004% (2,101 / 21,002 bp), so counting on the rounded column gives 3,044 and 749 for the first two intervals.",
                 "Gene-body coverage includes intronic sequence; exonic HERV coverage is reported separately (herv_bp_in_exons, exonic_herv_fraction, exonic_herv_percent). Of the 144 genes with > 80% gene-body coverage, 39 have > 80% exonic HERV coverage and 70 have <= 50%.",
                 "Locus shown in Figure 3C: ENSG00000286892."),
       sheets = list(list(name = "lncRNA_genes", src = "lnc_cov", drop = COVERAGE_DROP, lead = COVERAGE_LEAD, sort = "coverage"))),
  list(id = "pc_herv_coverage", scope = "pc",
       title = "HERV internal-region coverage of protein-coding genes overlapping same-strand HERV internal regions",
       current = "Supplemental Table 2 (lncRNA section; number also used for the RBP GSEA table)", publication = "Supplemental Table 5",
       numbering_date = "2026-10-06",
       notes = c("1,521 protein-coding genes (GENCODE v48, GRCh38) overlapping at least one same-strand HERV internal region; rows sorted by exact gene-body HERV coverage (descending), then exact exonic coverage (descending), then gene id.",
                 "Coverage uses the D028-corrected coordinate semantics (HERV BED intervals read as 0-based half-open).",
                 "Named examples in the text: ERVW-1, ERVH48-1, ERVFRD-1, ENSG00000293569, ENSG00000293570."),
       sheets = list(list(name = "Protein_coding_genes", src = "pc", drop = COVERAGE_DROP, lead = COVERAGE_LEAD, sort = "coverage"))),
  list(id = "lncrna_functional_annotation", scope = "fun",
       title = "Functional annotation of lncRNA genes overlapping same-strand HERV internal regions: LncRNAWiki curation, GO terms and PubMed literature mining",
       current = "not cited (workbook added at effective order 54.3, erratum D041)", publication = "Supplemental Table 6",
       numbering_date = "2026-10-11", release_date = "2026-10-11",
       notes = c("The 4,074 lncRNA genes of Supplemental Table 4, in the same order (exact gene-body HERV coverage, descending; then exact exonic coverage; then gene id).",
                 "A gene is annotated when it has a LncRNAWiki match, at least one PubMed literature match or a GO annotation: 593 genes (LncRNAWiki only 120, literature only 111, GO only 20, more than one source 342); 3,481 genes have none (annotation_sources = None).",
                 "functional_categories: the Figure 4B categories assigned to the gene by keyword patterns over the literature-derived terms and the LncRNAWiki fields (GO terms not used); categories are not mutually exclusive and are listed in the Figure 4B order. 381 annotated genes have at least one category; 173 are ceRNA-related.",
                 "Gene-body HERV coverage > 80%: 5 of 144 genes are annotated (CERNA2 and LINC02617 in LncRNAWiki; LINC02701, LINC02159 and LINC02009 by literature only).",
                 "lncrnawiki_mechanism and go_terms are shown with repeated items removed (first occurrence kept) and '; ' as separator; this changes 4 LncRNAWiki mechanism cells and no category, count or ceRNA assignment. LncRNAWiki portal export of 10 March 2026 (frozen).",
                 "Annotations describe functions reported for the lncRNA genes; they do not indicate that the HERV-derived sequence mediates them."),
       sheets = list(list(name = "lncRNA_annotation", src = c("fun_int", "fun_mem", "fun_cat", "fun_catn", "lnc_cov"), derive = "functional_annotation"))),
  list(id = "sparcs_like_catalogue", scope = "spc",
       title = "SPARCS-like candidate configurations: antisense LTRs in protein-coding 3' UTRs and lncRNA last exons, with GO Biological Process and MSigDB C7 over-representation, interferon-motif and U3-R-U5 tests",
       current = "Supplemental Table 3", publication = "Supplemental Table 7",
       numbering_date = "2026-10-11", release_date = "2026-10-11",
       notes = c("Catalogue sheets: one row per LTR x transcript where a HERV LTR overlaps (>= 1 bp) a GENCODE v48 transcript feature on the opposite strand. Coding_3UTR (protein-coding gene, 3' UTR): 3,084 rows, 1,913 unique LTRs, 1,399 genes. lncRNA_last_exon (lncRNA gene, last exon): 9,077 rows, 4,683 unique LTRs, 3,401 genes. 13 LTRs occur in both; 6,583 unique LTRs in total. LTR coordinates are 0-based half-open (BED); the LTR identifier gives 1-based coordinates.",
                 "GO_BP_coding_genes (Figure 5C): over-representation among the protein-coding genes of Coding_3UTR against all GENCODE v48 protein-coding genes (clusterProfiler enrichGO, BH; reported at p.adjust <= 0.05, q <= 0.20).",
                 "Immune_MSigDB_C7 (unchanged from the earlier C7 workbook). Gene-set collection: MSigDB 2026.1.Hs C7 (IMMUNESIGDB, VAX), frozen as release source SRC18 (msigdbr 26.1.0, offline). The MSigDB release of the historical analysis is UNKNOWN (erratum D034); these statistics are a regeneration under SRC18, not a reproduction of the historical table.",
                 "Immune_MSigDB_C7 foreground: protein-coding genes of the SPARCS-like coding branch; universe: all protein-coding gene names of the GTF (clusterProfiler enricher, BH; reported at p and p.adjust < 0.05, q < 0.20). Three signatures are significant.",
                 "IFN_motif_tests: the 1,903 Coding_3UTR LTRs present in the LTR motif scans against the 554,958 other scanned LTRs; 8 two-sided Fisher tests (STAT1 or STAT1::STAT2 by four interferon-tier definitions); all adjusted P = 1 (IRF class IRF1-IRF9, erratum D040). The accepted table also runs the tiers 1-4 comparison a second time under the mode name non_low, with identical counts and P; those 2 rows are not shown (author decision 2026-10-11). Its BH adjustment was over those 10 rows; over the 8 distinct tests it gives the same adjusted P (all 1).",
                 "U3RU5_comparison: the 6,583 SPARCS-like LTRs against all other LTRs of the U3-R-U5 catalogue (554,011; solo LTRs: 6,382 against 531,529); Fisher tests of the proportion of plausibly reconstructed (OK) loci and Wilcoxon tests of the confidence score. The reconstructed_both rows repeat the OK-proportion rows because both classifications select the same loci."),
       sheets = list(
         list(name = "Coding_3UTR", src = "sp_branch", where = c("branch", "coding_3UTR"), drop = TRIM_BRANCH, lead = SPARCS_LEAD),
         list(name = "lncRNA_last_exon", src = "sp_branch", where = c("branch", "lncRNA_last_exon"), drop = TRIM_BRANCH, lead = SPARCS_LEAD),
         list(name = "GO_BP_coding_genes", src = "sp_go", scope = "spgo"),
         list(name = "Immune_MSigDB_C7", src = "c7", scope = "c7"),
         list(name = "IFN_motif_tests", src = "sp_ifn", scope = "spifn", drop = TRIM_IFN, same_test = c("non_low", "tier1_2_3_4")),
         list(name = "U3RU5_comparison", src = "sp_u3", scope = "spu3", drop = TRIM_U3, lead = U3_LEAD))),
  list(id = "sparcs_tf_analysis", scope = "tf",
       title = "TF motifs depleted and enriched in lncRNA SPARCS-like antisense terminal LTRs, and GO analysis of the depleted TFs",
       current = "Supplemental Table 4", publication = "Supplemental Table 8",
       numbering_date = "2026-10-11",
       notes = c("Foreground: the 4,683 unique lncRNA-branch SPARCS-like LTRs (D033); background: SRC09 LTRs minus the foreground; Fisher's exact test per motif, BH over 730 motifs.",
                 "Depleted / Enriched sheets: motifs with sig = Depleted / Enriched (p_adj < 0.05), ordered by p_adj as upstream.",
                 "GO sheet: over-representation of the depleted TF symbols (clusterProfiler enrichGO, BH); terms are reported at p and p.adjust < 0.10, q < 0.20 (07_sparcs tf_go cutoffs); none reaches p.adjust < 0.05.",
                 "Supplemental Figure 2 is drawn from the same GO table; it is not built in this step."),
       sheets = list(
         list(name = "Depleted_TF_motifs", src = "tf_lnc", where = c("sig", "Depleted"), drop = TRIM_TF),
         list(name = "Enriched_TF_motifs", src = "tf_lnc", where = c("sig", "Enriched"), drop = TRIM_TF),
         list(name = "GO_BP_depleted_TFs", src = "tf_go")))
)

# ------------------------------------------------------------------------------------------------------------- helpers
kind_of <- function(product_scope, column_name) {
  # exactly one dictionary row for (product scope x column). Base-R logical indexing on the data.frame, so the arguments can never
  # be captured by the dictionary's own same-named columns (inside dict[...] data.table would read `scope` as dict$scope)
  hit <- which(dict[["scope"]] == product_scope & dict[["column"]] == column_name)
  if (length(hit) != 1) stop("expected exactly one dictionary row for ", product_scope, " / ", column_name, "; found ", length(hit))
  dict[hit]
}
as_kind <- function(x, kind, label) {
  # the typed XLSX value of one source column (the source strings themselves are never altered; MISSING_TOKEN -> blank cell)
  miss <- x == MISSING_TOKEN
  if (kind %in% c("id_text", "label_text")) return(ifelse(miss, NA_character_, x))
  if (kind == "logical") {
    if (!all(miss | x %in% c("TRUE", "FALSE"))) stop("non-logical value in ", label)
    return(ifelse(miss, NA, x == "TRUE"))
  }
  miss <- miss | x == ""  # numeric kinds only: an empty field is missing too (order 54.4)
  bad <- !miss & !grepl(NUMERIC_TOKEN, x)
  if (any(bad)) stop("non-numeric token in ", label, ": ", paste(head(x[bad], 3), collapse = ", "))
  v <- ifelse(miss, NA_real_, suppressWarnings(as.numeric(x)))
  if (kind == "integer" && any(!is.na(v) & v != round(v))) stop("non-integer value in ", label)
  if (kind == "p_value" && any(!is.na(v) & (v < 0 | v > 1))) stop("p-value-like value outside [0, 1] in ", label)
  v
}
zero_adjudication <- function(content_id, sheet_name, column_name, excel_row) {
  # an exact-zero p-value is displayed only with an adjudicated classification; there is no default text (never assume underflow)
  hit <- which(zero_adj[["content_id"]] == content_id & zero_adj[["sheet"]] == sheet_name & zero_adj[["column"]] == column_name &
               zero_adj[["excel_row"]] == as.character(excel_row))
  if (length(hit) != 1) stop("exact-zero p-value without an adjudicated classification: ", content_id, " / ", sheet_name, " / ",
                             column_name, " / Excel row ", excel_row, ". Audit it and record it in scripts/supplement/zero_p_value_adjudications.tsv")
  zero_adj[hit]
}

provenance <- list(); inventory <- list(); underflow <- list(); sheet_summary <- list(); p_audit <- list(); notes <- list()
hs <- createStyle(textDecoration = "bold", fgFill = "#E8EEF4", wrapText = TRUE, valign = "top")
st_p <- createStyle(numFmt = "0.00E+00"); st_int <- createStyle(numFmt = "0"); st_num <- createStyle(numFmt = "GENERAL")
st_txt <- createStyle(numFmt = "@")
NUMERIC_KINDS <- c("integer", "numeric", "p_value")

build_sheet <- function(prod, sh) {
  # select rows / columns of accepted tables; stacked sources get their constant label columns first
  if (!is.null(sh$derive)) {
    if (sh$derive != "functional_annotation") stop("unknown assembly function: ", sh$derive)
    return(derive_functional_annotation(sh))
  }
  parts <- lapply(seq_along(sh$src), function(i) {
    dt <- read_src(sh$src[i])
    if (!is.null(sh$add)) for (nm in rev(names(sh$add))) dt <- cbind(setNames(data.table(rep(sh$add[[nm]][i], nrow(dt))), nm), dt)
    dt
  })
  rows_per_source <- vapply(parts, nrow, 1L)
  dt <- rbindlist(parts, use.names = TRUE)
  if (!is.null(sh$where)) dt <- dt[get(sh$where[1]) == sh$where[2]]
  if (!is.null(sh$same_test)) {  # the accepted IFN table runs the tiers 1-4 comparison twice, under two mode names
    dup <- grepl(paste0("\\| mode=", sh$same_test[1], "$"), dt$test)
    twin <- sub(paste0("mode=", sh$same_test[1], "$"), paste0("mode=", sh$same_test[2]), dt$test[dup])
    other <- setdiff(names(dt), "test")
    same <- vapply(seq_along(twin), function(i) {
      a <- dt[which(dup)[i], ..other]; b <- dt[test == twin[i], ..other]
      nrow(b) == 1 && identical(unlist(a), unlist(b))
    }, TRUE)
    if (!any(dup) || !all(same)) stop("same_test: a '", sh$same_test[1], "' row is not identical to its '", sh$same_test[2], "' row")
    dt <- dt[!dup]
  }
  missing_drop <- setdiff(sh$drop, names(dt))
  if (length(missing_drop)) stop("columns to drop are absent from the source: ", paste(missing_drop, collapse = ", "))
  keep <- setdiff(names(dt), sh$drop)
  if (!is.null(sh$lead)) {  # presentation: these columns first, the rest in source order (author decision 2026-10-09)
    if (!all(sh$lead %in% keep)) stop("lead columns absent: ", paste(setdiff(sh$lead, keep), collapse = ", "))
    keep <- c(sh$lead, setdiff(keep, sh$lead))
  }
  dt <- dt[, ..keep]
  if (identical(sh$sort, "coverage")) dt <- sort_by_coverage(dt)
  list(data = dt, rows_per_source = rows_per_source)
}

sort_by_coverage <- function(dt) {
  # Tables 4 and 5 (author decision 2026-10-09): exact gene-body HERV coverage herv_bp_in_gene_body / gene_body_bp descending,
  # ties by exact exonic coverage herv_bp_in_exons / exonic_bp descending, then gene_stable_id ascending (C locale). Values are
  # untouched; only the row order changes. The double ratios order exactly here (integer operands below 2^26, so distinct
  # fractions never round to the same double); the result is re-checked with exact integer cross-multiplication.
  gb <- as.numeric(dt$herv_bp_in_gene_body); gl <- as.numeric(dt$gene_body_bp)
  ex <- as.numeric(dt$herv_bp_in_exons); el <- as.numeric(dt$exonic_bp)
  if (any(gl <= 0 | el <= 0) || max(c(gb, gl, ex, el)) >= 2^26) stop("coverage sort: lengths must be positive and below 2^26")
  o <- order(-(gb / gl), -(ex / el), dt$gene_stable_id, method = "radix")
  dt <- dt[o]; gb <- gb[o]; gl <- gl[o]; ex <- ex[o]; el <- el[o]
  if (nrow(dt) > 1) {
    i <- seq_len(nrow(dt) - 1); j <- i + 1
    c1 <- gb[i] * gl[j] - gb[j] * gl[i]; c2 <- ex[i] * el[j] - ex[j] * el[i]  # exact: products below 2^52
    bad <- c1 < 0 | (c1 == 0 & (c2 < 0 | (c2 == 0 & dt$gene_stable_id[i] >= dt$gene_stable_id[j])))
    if (any(bad)) stop("coverage sort: order not confirmed by exact integer comparison")
  }
  dt
}

dedup_items <- function(x, sep) {
  # display only: repeated items of a multi-value cell removed (first occurrence kept), re-joined with "; " (author decision
  # 2026-10-11, order 54.3); the item set is unchanged and MISSING_TOKEN stays missing
  vapply(x, function(v) if (v == MISSING_TOKEN) v else paste(unique(strsplit(v, sep, fixed = TRUE)[[1]]), collapse = "; "), "",
         USE.NAMES = FALSE)
}

derive_functional_annotation <- function(sh) {
  # Supplemental Table 6 (effective order 54.3, D041; author decisions 2026-10-11): one row per HERV-associated lncRNA gene, in the
  # Table 4 order, assembled from the accepted order-42 tables. Nothing is classified here: the source flags, categories, mechanism
  # and GO fields and literature PMIDs are the accepted values, re-presented as follows
  #   annotation_sources    the names of the accepted has_lncRNAwiki / has_lit / has_go flags equal to 1, "; "-joined; None if no flag
  #   functional_categories the gene's accepted Figure 4B categories (long table), "; "-joined in the accepted Figure 4B order
  #   lncrnawiki_mechanism  LncRNAWiki functional_mechanisms, dedup_items; the 06 placeholder for unmatched genes -> missing
  #   go_terms              GO term names, dedup_items
  #   gene_body_herv_percent from the Table 4 source (lnc_cov, order 35), so both tables show the same token (06 writes 34 whole
  #                         numbers without the ".0" of 05, e.g. 100 for 100.0); the two must agree numerically for every gene
  # Every flag is cross-checked against the integrated table's own evidence columns; any disagreement stops the build.
  int <- read_src("fun_int"); mem <- read_src("fun_mem"); gcat <- read_src("fun_cat"); catn <- read_src("fun_catn"); cov <- read_src("lnc_cov")
  if (anyDuplicated(int$gene_stable_id) || anyDuplicated(mem$gene_stable_id) || !setequal(int$gene_stable_id, mem$gene_stable_id))
    stop("functional annotation: the integrated and membership tables must hold the same unique genes")
  if (anyDuplicated(catn$category) || !all(gcat$category %in% catn$category) || !all(gcat$gene_stable_id %in% int$gene_stable_id) ||
      anyDuplicated(paste(gcat$gene_stable_id, gcat$category)))
    stop("functional annotation: the Figure 4B gene-category table does not match its category list or the gene set")
  int <- sort_by_coverage(int)
  m <- mem[match(int$gene_stable_id, mem$gene_stable_id)]
  if (any(m$gene_name != int$gene_name)) stop("functional annotation: gene names differ between the integrated and membership tables")
  if (anyDuplicated(cov$gene_stable_id) || !setequal(cov$gene_stable_id, int$gene_stable_id))
    stop("functional annotation: the gene set differs from the Table 4 coverage source")
  pct <- cov$gene_body_herv_percent[match(int$gene_stable_id, cov$gene_stable_id)]
  if (any(as.numeric(pct) != as.numeric(int$gene_body_herv_percent))) stop("functional annotation: gene-body percentages differ from the Table 4 source")
  flags <- c("has_lncRNAwiki", "has_lit", "has_go", "has_any_ann")
  if (!all(unlist(m[, ..flags]) %in% c("0", "1"))) stop("functional annotation: membership flags must be 0/1")
  wiki <- int$lncrnawiki_annotation_status == "matched"; lit <- int$n_literature_rows != MISSING_TOKEN; go <- int$has_go_annotation == "TRUE"
  if (any((m$has_lncRNAwiki == "1") != wiki | (m$has_lit == "1") != lit | (m$has_go == "1") != go | (m$has_any_ann == "1") != (wiki | lit | go)))
    stop("functional annotation: membership flags disagree with the integrated table's evidence columns")
  if (any(int$functional_mechanisms[!wiki] != "no lncRNAwiki annotation") || any(int$functional_mechanisms[wiki] == "no lncRNAwiki annotation"))
    stop("functional annotation: unexpected LncRNAWiki placeholder pattern")
  if (any(int$go_terms[!go] != MISSING_TOKEN) || any(int$go_terms[go] == MISSING_TOKEN) ||
      any((int$n_pmids != MISSING_TOKEN) != lit) || any((int$pmids_literature != MISSING_TOKEN) != lit))
    stop("functional annotation: GO or literature fields present for a gene without that evidence (or absent with it)")
  n_listed <- lengths(strsplit(int$pmids_literature[lit], "; ", fixed = TRUE))
  if (any(n_listed != as.integer(int$n_pmids[lit]))) stop("functional annotation: n_pmids differs from the number of listed PMIDs")
  src_names <- c("LncRNAWiki", "Literature", "GO")
  sources <- vapply(seq_len(nrow(m)), function(i) {
    s <- src_names[c(m$has_lncRNAwiki[i], m$has_lit[i], m$has_go[i]) == "1"]
    if (length(s)) paste(s, collapse = "; ") else "None"
  }, "")
  by_gene <- split(gcat$category, gcat$gene_stable_id)
  cats <- vapply(int$gene_stable_id, function(g) {
    x <- by_gene[[g]]
    if (is.null(x)) MISSING_TOKEN else paste(catn$category[catn$category %in% x], collapse = "; ")
  }, "", USE.NAMES = FALSE)
  if (any(cats != MISSING_TOKEN & m$has_any_ann != "1")) stop("functional annotation: a category assigned to an unannotated gene")
  mech <- ifelse(wiki, int$functional_mechanisms, MISSING_TOKEN)
  dt <- data.table(gene_stable_id = int$gene_stable_id, gene_name = int$gene_name, gene_body_herv_percent = pct,
                   annotation_sources = sources, functional_categories = cats,
                   lncrnawiki_mechanism = dedup_items(mech, ";"), go_terms = dedup_items(int$go_terms, ";"),
                   n_literature_pmids = int$n_pmids, literature_pmids = int$pmids_literature)
  INT <- SRC$fun_int[3]; MEM <- SRC$fun_mem[3]
  prov <- list(
    gene_stable_id = c(paste(INT, MEM, sep = ";"), "gene_stable_id", "join key (1:1, identical gene sets); rows in the Table 4 order (exact gene-body coverage, exonic coverage, gene id)"),
    gene_name = c(paste(INT, MEM, sep = ";"), "gene_name", "passthrough (identical in both tables)"),
    gene_body_herv_percent = c(SRC$lnc_cov[3], "gene_body_herv_percent", "passthrough from the Table 4 source (numerically equal to the 06 value, checked)"),
    annotation_sources = c(MEM, "has_lncRNAwiki;has_lit;has_go", "names of the flags equal to 1 (LncRNAWiki; Literature; GO), '; '-joined; None when all are 0; flags cross-checked against the integrated table"),
    functional_categories = c(paste(SRC$fun_cat[3], SRC$fun_catn[3], sep = ";"), "category", "the gene's categories from the long table, '; '-joined in the order of the category table; missing when none"),
    lncrnawiki_mechanism = c(INT, "functional_mechanisms", "repeated items removed (first kept), '; '-joined; the placeholder 'no lncRNAwiki annotation' -> missing"),
    go_terms = c(INT, "go_terms", "repeated items removed (first kept), '; '-joined"),
    n_literature_pmids = c(INT, "n_pmids", "passthrough"),
    literature_pmids = c(INT, "pmids_literature", "passthrough"))
  list(data = dt, rows_per_source = nrow(dt), prov = prov)
}

write_product <- function(prod) {
  wb <- createWorkbook(creator = "ncRNA release")
  modifyBaseFont(wb, fontSize = 10, fontName = "Arial")
  tsv_dir <- file.path(OUT, "supplement/tables/by_content_tsv", prod$id)
  dir.create(tsv_dir, recursive = TRUE, showWarnings = FALSE)
  readme_cols <- list(); exact <- list(); p_min <- NULL
  for (sh in prod$sheets) {
    b <- build_sheet(prod, sh); dt <- b$data
    srcs <- vapply(sh$src, function(k) SRC[[k]][3], "")
    # raw TSV support file: the accepted source strings, unparsed
    fwrite(dt, file.path(tsv_dir, paste0(sh$name, ".tsv")), sep = "\t", quote = FALSE, eol = "\n")
    out <- as.data.frame(dt, stringsAsFactors = FALSE)
    kinds <- character(0)
    for (cl in names(out)) {
      k <- kind_of(if (is.null(sh$scope)) prod$scope else sh$scope, cl)
      kinds[cl] <- k$kind
      out[[cl]] <- as_kind(out[[cl]], k$kind, paste(prod$id, sh$name, cl))
      src_col <- if (!is.null(sh$add) && cl %in% names(sh$add)) "(constant label added in module 08)" else cl
      pv <- b$prov[[cl]]  # assembled sheet: per-column source, source column(s) and presentation step
      provenance[[length(provenance) + 1]] <<- data.table(
        content_id = prod$id, sheet = sh$name, column = cl, source_column = if (is.null(pv)) src_col else pv[2],
        source_file = if (is.null(pv)) paste(srcs, collapse = ";") else pv[1],
        transformation = if (!is.null(pv)) pv[3] else if (src_col == cl) "passthrough" else paste0("constant per stacked source: ", paste(sh$add[[cl]], collapse = " | ")),
        cell_type = k$kind,
        display_format = switch(k$kind, p_value = "0.00E+00 (exact zero -> adjudicated text only)", integer = "0", numeric = "General",
                                logical = "TRUE/FALSE", "text (@)"))
      readme_cols[[length(readme_cols) + 1]] <- c(sh$name, cl, k$definition)
    }
    for (cl in sh$drop) {
      provenance[[length(provenance) + 1]] <<- data.table(
        content_id = prod$id, sheet = sh$name, column = "(not exported)", source_column = cl, source_file = paste(srcs, collapse = ";"),
        transformation = if (cl %in% DISPLAY_COLS) paste0("dropped: upstream 4-significant-digit display string; the full-precision numeric column ",
                                                          sub("_display$", "", cl), " is exported") else if (prod$id == "rbp_gsea") TRIM_NOTE_T2
                         else if (prod$id %in% c("lncrna_herv_coverage", "pc_herv_coverage")) TRIM_NOTE_T45
                         else if (prod$id == "sparcs_tf_analysis") TRIM_NOTE_T7
                         else if (prod$id == "sparcs_like_catalogue") TRIM_NOTE_SPC
                         else if (prod$id == "rbp_threshold_sensitivity" && !(cl %in% TRIM_THRESHOLD)) TRIM_NOTE_T3 else TRIM_NOTE,
        cell_type = "none", display_format = "none")
    }
    addWorksheet(wb, sh$name)
    # p-value-like columns: audit (rows, NA, zero, smallest positive) and exact-zero cells
    zero_cells <- list()
    for (cl in names(out)[kinds == "p_value"]) {
      v <- out[[cl]]; pos <- which(!is.na(v) & v > 0); z <- which(!is.na(v) & v == 0)
      min_tok <- if (length(pos)) dt[[cl]][pos[which.min(v[pos])]] else "NA"
      p_audit[[length(p_audit) + 1]] <<- data.table(content_id = prod$id, sheet = sh$name, column = cl, n_rows = length(v),
                                                    n_na = sum(is.na(v)), n_zero = length(z), min_positive = min_tok,
                                                    max_value = if (any(!is.na(v))) dt[[cl]][which.max(v)] else "NA")
      if (length(pos) && (is.null(p_min) || v[pos[which.min(v[pos])]] < p_min$value))
        p_min <- list(value = min(v[pos]), token = min_tok, where = paste0(sh$name, " / ", cl))
      if (length(z)) zero_cells[[cl]] <- z
    }
    writeData(wb, sh$name, out, headerStyle = hs, keepNA = FALSE)
    for (cl in names(zero_cells)) {
      j <- match(cl, names(out))
      for (i in zero_cells[[cl]]) {
        a <- zero_adjudication(prod$id, sh$name, cl, i + 1)
        writeData(wb, sh$name, a$display_text, startCol = j, startRow = i + 1)
        underflow[[length(underflow) + 1]] <<- data.table(content_id = prod$id, sheet = sh$name, column = cl, excel_row = i + 1,
                                                          source_value = dt[[cl]][i], classification = a$classification,
                                                          display = a$display_text)
      }
    }
    n <- nrow(out)
    for (j in seq_along(out)) {
      st <- switch(kinds[[j]], p_value = st_p, integer = st_int, numeric = st_num, logical = NULL, st_txt)
      if (!is.null(st) && n > 0) addStyle(wb, sh$name, st, rows = 2:(n + 1), cols = j, gridExpand = TRUE)
    }
    setColWidths(wb, sh$name, cols = seq_along(out), widths = pmin(60, pmax(10, nchar(names(out)) + 2)))
    freezePane(wb, sh$name, firstRow = TRUE)
    # exact source tokens of the numeric cells, written into the sheet XML by normalise_xlsx (adjudicated zeros are text cells)
    num_cols <- which(kinds %in% NUMERIC_KINDS)
    tokens <- lapply(num_cols, function(j) {
      x <- dt[[j]]; x[x == MISSING_TOKEN | x == ""] <- NA
      if (names(out)[j] %in% names(zero_cells)) x[zero_cells[[names(out)[j]]]] <- NA
      x
    })
    exact[[sh$name]] <- list(cols = num_cols, tokens = tokens)
    for (s in unique(sh$src)) {
      inventory[[length(inventory) + 1]] <<- data.table(
        content_id = prod$id, canonical_file = paste0("supplement/tables/by_content/", prod$id, ".xlsx"), canonical_sheet = sh$name,
        support_tsv = paste0("supplement/tables/by_content_tsv/", prod$id, "/", sh$name, ".tsv"),
        upstream_module = SRC[[s]][1], upstream_order = SRC[[s]][2], upstream_path = SRC[[s]][3], upstream_sha256 = manifest_sha(SRC[[s]][3]),
        source_line_terminator = source_terminator(s), support_tsv_line_terminator = "LF",
        rows_in_sheet = n, snapshot_manuscript_label = prod$current, publication_label = prod$publication,
        numbering_status = NUMBERING_STATUS)
    }
    sheet_summary[[length(sheet_summary) + 1]] <<- data.table(content_id = prod$id, sheet = sh$name, n_rows = n, n_cols = ncol(out),
                                                              header = paste(names(out), collapse = "|"))
  }
  # publication metadata -> release record (not published in the workbook)
  p_rule <- if (is.null(p_min)) "This workbook has no p-value-like column." else if (length(underflow) == 0 || !any(vapply(underflow, function(u) u$content_id == prod$id, TRUE)))
    paste0("p-value-like columns are shown in scientific format; no p-value-like cell in this workbook is zero (smallest positive value ",
           p_min$token, ", ", p_min$where, "). A blank p-value cell means NA in the accepted source (no valid test).")
  else paste0("p-value-like columns are shown in scientific format; exact-zero source values are shown with their adjudicated text ",
              "(panel_sources/supplement/underflow_display_cells.tsv). A blank p-value cell means NA in the accepted source (no valid test).")
  r <- list(c("title", prod$title), c("reference", "GRCh38; GENCODE v48 where gene annotation is used"))
  for (i in seq_along(prod$notes)) r[[length(r) + 1]] <- c(paste0("note_", i), prod$notes[i])
  r[[length(r) + 1]] <- c("display_rule", paste0("Each numeric cell holds the accepted value exactly as written in the release ",
                                                 "(supplement/tables/by_content_tsv/", prod$id, "/ has the same values as text). ", p_rule,
                                                 " Identifiers are stored as text."))
  for (x in r) notes[[length(notes) + 1]] <<- data.table(
    content_id = prod$id, publication_label = prod$publication,
    numbering_date = if (is.null(prod$numbering_date)) "2026-10-01" else prod$numbering_date,
    release_date = if (is.null(prod$release_date)) "2026-10-01" else prod$release_date,
    snapshot_manuscript_label = prod$current, field = x[1], value = x[2])
  # the only non-data sheet: column definitions
  addWorksheet(wb, "Column_definitions")
  cols <- rbindlist(lapply(readme_cols, function(x) data.table(Sheet = x[1], Column = x[2], Definition = x[3])))
  writeData(wb, "Column_definitions", cols, headerStyle = hs)
  setColWidths(wb, "Column_definitions", cols = 1:3, widths = c(26, 34, 100))
  addStyle(wb, "Column_definitions", createStyle(wrapText = TRUE, valign = "top"), rows = 1:(nrow(cols) + 1), cols = 1:3, gridExpand = TRUE, stack = TRUE)
  freezePane(wb, "Column_definitions", firstRow = TRUE)
  raw <- file.path(WORK, paste0(prod$id, ".raw.xlsx"))
  saveWorkbook(wb, raw, overwrite = TRUE)
  normalise_xlsx(raw, file.path(OUT, "supplement/tables/by_content", paste0(prod$id, ".xlsx")), exact)
}

col_letters_to_index <- function(s) vapply(strsplit(s, ""), function(ch) sum((match(ch, LETTERS)) * 26^(rev(seq_along(ch)) - 1)), 1)
index_to_col_letters <- function(j) { s <- ""; while (j > 0) { r <- (j - 1) %% 26; s <- paste0(LETTERS[r + 1], s); j <- (j - 1) %/% 26 }; s }

sheet_xml_paths <- function(d) {
  # sheet name -> worksheet part, from xl/workbook.xml and its relationships (no assumption about openxlsx part numbering)
  wbx <- paste(readLines(file.path(d, "xl/workbook.xml"), warn = FALSE, encoding = "UTF-8"), collapse = "")
  rel <- paste(readLines(file.path(d, "xl/_rels/workbook.xml.rels"), warn = FALSE), collapse = "")
  sh <- regmatches(wbx, gregexpr("<sheet [^>]*/>", wbx))[[1]]
  nm <- sub(".* name=\"([^\"]*)\".*", "\\1", sh); rid <- sub(".* r:id=\"([^\"]*)\".*", "\\1", sh)
  rl <- regmatches(rel, gregexpr("<Relationship [^>]*/>", rel))[[1]]
  rl_id <- sub(".* Id=\"([^\"]*)\".*", "\\1", rl); rl_tg <- sub(".* Target=\"([^\"]*)\".*", "\\1", rl)
  setNames(file.path(d, "xl", sub("^/?xl/", "", rl_tg[match(rid, rl_id)])), nm)
}

resolve_part <- function(rels_member, target) {
  # package path of a relationship target, relative to the part that owns the .rels file (OPC rules)
  owner_dir <- dirname(dirname(rels_member))
  segs <- strsplit(if (startsWith(target, "/")) sub("^/", "", target) else file.path(owner_dir, target), "/", fixed = TRUE)[[1]]
  out <- character(0)
  for (sg in segs[segs != "." & segs != ""]) out <- if (sg == "..") head(out, -1) else c(out, sg)
  paste(out, collapse = "/")
}

drop_dangling_parts <- function(d, members) {
  # openxlsx 4.2.5 declares drawing / vmlDrawing parts it never writes (worksheet .rels and [Content_Types].xml overrides);
  # a relationship to a missing part makes the package invalid (openpyxl refuses it; Excel may offer a repair). Remove them.
  for (rm in members[grepl("\\.rels$", members)]) {
    f <- file.path(d, rm); x <- readChar(f, file.size(f), useBytes = TRUE)
    rel <- regmatches(x, gregexpr("<Relationship [^>]*/>", x))[[1]]
    tgt <- sub(".* Target=\"([^\"]*)\".*", "\\1", rel); ext <- grepl("TargetMode=\"External\"", rel)
    dangling <- !ext & !vapply(tgt, function(t) resolve_part(rm, t) %in% members, TRUE)
    for (r in rel[dangling]) x <- sub(r, "", x, fixed = TRUE)
    writeChar(x, f, eos = NULL, useBytes = TRUE)
  }
  f <- file.path(d, "[Content_Types].xml"); x <- readChar(f, file.size(f), useBytes = TRUE)
  ov <- regmatches(x, gregexpr("<Override [^>]*/>", x))[[1]]
  part <- sub("^/", "", sub(".* PartName=\"([^\"]*)\".*", "\\1", ov))
  for (o in ov[!part %in% members]) x <- sub(o, "", x, fixed = TRUE)
  writeChar(x, f, eos = NULL, useBytes = TRUE)
}

normalise_xlsx <- function(raw, final, exact) {
  # (1) numeric data cells: openxlsx's 15-significant-digit fixed-point text -> the accepted source token (bit-exact, scientific)
  # (2) <dimension> set to the used range (openxlsx writes "A1")
  # (3) relationships / content-type overrides to parts absent from the package are removed (openxlsx 4.2.5 defect)
  # (4) deterministic container: fixed creation date in docProps/core.xml, fixed member times (UTC), fixed member order
  d <- file.path(WORK, paste0(basename(final), ".d"))
  unlink(d, recursive = TRUE); dir.create(d)
  members <- unzip(raw, list = TRUE)$Name
  unzip(raw, exdir = d)
  drop_dangling_parts(d, members)
  parts <- sheet_xml_paths(d)
  cell_pat <- "<c r=\"([A-Z]+)([0-9]+)\"( s=\"[0-9]+\")? t=\"n\"><v>([^<]*)</v></c>"
  for (sn in names(parts)) {
    f <- parts[[sn]]
    x <- readChar(f, file.size(f), useBytes = TRUE)
    if (!is.null(exact[[sn]])) {
      e <- exact[[sn]]
      m <- gregexpr(cell_pat, x, perl = TRUE)
      cells <- regmatches(x, m)[[1]]
      cl <- sub(cell_pat, "\\1", cells, perl = TRUE); rw <- as.integer(sub(cell_pat, "\\2", cells, perl = TRUE))
      sty <- sub(cell_pat, "\\3", cells, perl = TRUE); old <- sub(cell_pat, "\\4", cells, perl = TRUE)
      if (any(rw < 2)) stop("numeric cell in the header row of ", sn)
      jj <- match(col_letters_to_index(cl), e$cols)
      if (anyNA(jj)) stop("numeric cell outside the numeric columns of ", sn)
      tok <- mapply(function(j, i) e$tokens[[j]][i], jj, rw - 1)
      if (anyNA(tok)) stop("numeric cell without a source token in ", sn)
      expected <- sum(vapply(e$tokens, function(t) sum(!is.na(t)), 1))
      if (length(cells) != expected) stop("numeric cell count ", length(cells), " != source token count ", expected, " in ", sn)
      ov <- as.numeric(old); tv <- as.numeric(tok)
      if (any(abs(ov - tv) > 1e-14 * abs(tv))) stop("openxlsx value and source token disagree in ", sn)
      regmatches(x, m) <- list(paste0("<c r=\"", cl, rw, "\"", sty, " t=\"n\"><v>", tok, "</v></c>"))
    }
    refs <- regmatches(x, gregexpr("<c r=\"[A-Z]+[0-9]+\"", x))[[1]]
    if (length(refs)) {
      max_row <- max(as.integer(sub("<c r=\"[A-Z]+([0-9]+)\"", "\\1", refs)))
      max_col <- max(col_letters_to_index(sub("<c r=\"([A-Z]+)[0-9]+\"", "\\1", refs)))
      x <- sub("<dimension ref=\"A1\"/>", paste0("<dimension ref=\"A1:", index_to_col_letters(max_col), max_row, "\"/>"), x, fixed = TRUE)
    }
    writeChar(x, f, eos = NULL, useBytes = TRUE)
  }
  core <- file.path(d, "docProps/core.xml")
  x <- readLines(core, warn = FALSE, encoding = "UTF-8")
  x <- gsub("<dcterms:created xsi:type=\"dcterms:W3CDTF\">[^<]*</dcterms:created>",
            paste0("<dcterms:created xsi:type=\"dcterms:W3CDTF\">", RELEASE_DATE, "T00:00:00Z</dcterms:created>"), x)
  writeLines(x, core, sep = "", useBytes = TRUE)
  for (m in members) Sys.setFileTime(file.path(d, m), FIXED_TIME)
  for (dd in unique(dirname(members))) if (dd != ".") Sys.setFileTime(file.path(d, dd), FIXED_TIME)
  unlink(final)
  zip::zip(normalizePath(final, mustWork = FALSE), files = members, root = d, mode = "mirror", include_directories = FALSE)
  unlink(d, recursive = TRUE)
}

for (prod in P) write_product(prod)

# ------------------------------------------------------------------------------------------------------------- source archive
# Supplemental Table X source files: the 11 historical members, from their canonical release paths (scripted, deterministic)
arch_keys <- c("t_global", "t_family", "t_per_rbp", "t_sets", "t_sparsity", "go_h_1e4", "go_k_1e4", "gsea_1e4", "go_h_1e5", "go_k_1e5", "gsea_1e5")
arch_files <- sort(vapply(arch_keys, function(k) { read_src(k); SRC[[k]][3] }, ""), method = "radix")  # C-locale byte order
# a release record since order 54.2 (author decision 2026-10-09): written to panel_sources/supplement, not to the published supplement/
tarball <- normalizePath(file.path(OUT, "panel_sources/supplement/rbp_threshold_sensitivity_source_files.tar"), mustWork = FALSE)
st <- system2("tar", c("--sort=name", "--mtime=@1790812800", "--owner=0", "--group=0", "--numeric-owner", "--mode=go-w,a+r",
                       "-C", shQuote(REL), "-cf", shQuote(tarball), shQuote(arch_files)))
if (st != 0) stop("tar failed")
st <- system2("gzip", c("-n", "-9", "-f", shQuote(tarball)))
if (st != 0) stop("gzip failed")
archive_rows <- data.table(member = arch_files, member_sha256 = vapply(arch_files, manifest_sha, ""))

# ------------------------------------------------------------------------------------------------------------- tables
PS <- file.path(OUT, "panel_sources/supplement")
w <- function(x, f) fwrite(x, file.path(PS, f), sep = "\t", quote = FALSE, na = "NA", eol = "\n")
w(rbindlist(provenance), "supplement_column_provenance.tsv")
w(rbindlist(inventory), "supplement_table_inventory.tsv")
w(rbindlist(sheet_summary), "supplement_sheet_summary.tsv")
w(rbindlist(p_audit), "supplement_p_value_audit.tsv")
w(rbindlist(notes), "supplement_workbook_notes.tsv")
uf <- if (length(underflow)) rbindlist(underflow) else
  data.table(content_id = character(), sheet = character(), column = character(), excel_row = integer(), source_value = character(),
             classification = character(), display = character())
w(uf, "underflow_display_cells.tsv")
w(archive_rows, "rbp_threshold_sensitivity_archive_members.tsv")

# ------------------------------------------------------------------------------------------------------------- checksums
write_checksums <- function(dir) {
  files <- sort(setdiff(list.files(dir, recursive = TRUE), "checksums.sha256"), method = "radix")
  writeLines(paste0(vapply(file.path(dir, files), file_sha, ""), "  ./", files), file.path(dir, "checksums.sha256"))
}
write_checksums(file.path(OUT, "supplement")); write_checksums(PS)
cat("[OK] ", length(P), " workbooks, ", sum(vapply(P, function(p) length(p$sheets), 1L)), " data sheets, archive of ",
    nrow(archive_rows), " members; exact-zero p-value cells displayed: ", nrow(uf), "\n", sep = "")
