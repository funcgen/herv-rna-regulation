#!/usr/bin/env Rscript
# ==============================================================================
# 06_lncrna_function step 4 (frozen global order 41, effective order 42), LOCAL, offline:
# integrate curated (LncRNAWiki), GO and literature annotation with the canonical HERV coverage
# categories, and write the Figure 4 source tables.
#
# REFACTOR_SCRIPT of G_literature_integration/exploration_lit_mining.R (sha256 393780db...), sections
# 7, 8, 9, 11, 12 and the Panel A/B/C data-frame code, copied verbatim apart from the RELEASE edits:
#   D032  the literature matching (sections 1, 3-6, 10) is effective order 41: lit_summary_per_gene is
#         read from results/literature_classification/, never recomputed here
#   D029  herv_contribution is the order-37 herv_coverage_category (the historical rounded-percent
#         case_when is not re-applied)
#   the HERV lncRNA table is the canonical order-38 full table (TSV) instead of res_cov.RDS
#   no plotting (08_figures_supplement), no exploratory prints; herv_lncrna_simplified is a TSV
# Rules (historical): has_lncRNAwiki = status != placeholder; has_lit = n_literature_rows > 0;
# has_go = has_go_annotation; has_any_ann = has_lncRNAwiki | has_lit | has_go. Panel A classes by the
# three flags; Panel B categories by regex over the lower-cased union of functions_detected and five
# LncRNAWiki fields, among annotated genes (ceRNA = 'cerna|competing endogenous rna|mirna sponge|
# microrna sponge|mirna'); Panel C = annotated / total per coverage category.
# ==============================================================================
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(forcats)
  library(scales)
  library(data.table)
  library(digest)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("missing argument ", flag)
  args[[i + 1]]
}
ANNOT <- get_arg("--annotation")
LITSUM <- get_arg("--literature-summary")
CATEG <- get_arg("--categories")
FINAL_OUT <- get_arg("--outdir")

verify_manifest <- function(path) {
  m <- read.table(file.path(dirname(path), "checksums.sha256"), col.names = c("sha", "file"), colClasses = "character")
  want <- m$sha[m$file %in% c(basename(path), paste0("./", basename(path)))]
  if (length(want) != 1 || digest(file = path, algo = "sha256") != want) stop(path, " does not match its accepted manifest")
  message("[ok] ", basename(path), " sha256 ", want, " (manifest)")
}
for (p in c(ANNOT, LITSUM, CATEG)) verify_manifest(p)

# ------------------------------------------------------------------------------
# RELEASE: canonical inputs (order 38 table, order 41 per-gene literature summary, order 37 categories)
# ------------------------------------------------------------------------------
herv_lncrna <- fread(ANNOT, sep = "\t", quote = "", na.strings = "NA", data.table = FALSE, encoding = "UTF-8")
lit_summary_per_gene <- read_tsv(LITSUM, col_types = cols(.default = col_character(), n_literature_rows = col_integer(),
                                                         n_pmids = col_integer(), min_candidates_in_sentence = col_double(),
                                                         max_candidates_in_sentence = col_double()),
                                 quote = "", na = "NA", progress = FALSE)
categories <- fread(CATEG, sep = "\t", quote = "", data.table = FALSE, encoding = "UTF-8")
if (nrow(herv_lncrna) != 4074 || anyDuplicated(herv_lncrna$gene_stable_id)) stop("order-38 universe must be 4,074 unique genes")
if (!all(lit_summary_per_gene$gene_stable_id %in% herv_lncrna$gene_stable_id)) stop("literature genes outside the universe")
cc <- categories$herv_coverage_category[match(herv_lncrna$gene_stable_id, categories$gene_stable_id)]
if (anyNA(cc) || !identical(cc, herv_lncrna$herv_coverage_category)) stop("order-38 categories differ from order 37")
if (!all(cc %in% c("minor", "substantial", "dominant", "predominant"))) stop("invalid coverage category")

parent <- dirname(normalizePath(FINAL_OUT, mustWork = FALSE))
dir.create(parent, recursive = TRUE, showWarnings = FALSE)
STAGE <- tempfile(pattern = ".stage_integrated_annotation.", tmpdir = parent)
dir.create(STAGE)
options(error = function() { unlink(STAGE, recursive = TRUE); quit(save = "no", status = 1) })   # no partial output on failure

# ============================================================
# 7. Merge literature summary back into HERV table
# ============================================================

merge_pmids <- function(x, sep = "; ") {
  x <- x[!is.na(x) & x != "" & x != "no lncRNAwiki annotation"]
  if (length(x) == 0) return(NA_character_)
  
  x %>%
    paste(collapse = ";") %>%
    str_split(";") %>%
    unlist() %>%
    str_trim() %>%
    discard(~ .x == "" || .x == "NA" || .x == "no lncRNAwiki annotation") %>%
    unique() %>%
    sort() %>%
    paste(collapse = sep)
}


herv_with_literature <- herv_lncrna %>%
  left_join(lit_summary_per_gene, by = c("gene_stable_id", "gene_name")) %>%
  mutate(
    pmids_lncrnawiki = pmids,
    pmids_merged = purrr::map2_chr(pmids_lncrnawiki, pmids_literature, ~ merge_pmids(c(.x, .y)))
  )

# ============================================================
# 8. Build comparison table: lncRNAWiki vs literature
# ============================================================

check_tbl <- herv_with_literature %>%
  mutate(
    has_lncrnawiki = n_lncrnawiki_rows > 0,
    has_literature = !is.na(n_pmids) & n_pmids > 0
  )


# ============================================================
# 9. Inspect literature-only genes
# ============================================================

lit_only <- check_tbl %>%
  filter(!has_lncrnawiki, has_literature)

lit_only_inspection <- lit_only %>%
  dplyr::select(
    gene_stable_id,
    gene_name,
    matched_symbols,
    matched_synonyms,
    n_pmids,
    n_literature_rows,
    evidence_tiers,
    sentence_types,
    mechanistic_terms,
    phenotypic_terms,
    associative_terms,
    pathway_terms,
    regulatory_verbs,
    functions_detected,
    evidence_sentences
  )

write.table(
  lit_only_inspection,
  file = file.path(STAGE, "lit_only_inspection.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

# ============================================================
# 11. Integrated final table
# ============================================================

write.table(
  herv_with_literature,
  file = file.path(STAGE, "herv_lncrna_with_literature.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)



# ============================================================
# 12. Simplification: generate smaller table
# ============================================================

herv_simplified <- herv_with_literature %>% 
  mutate(
    has_lncRNAwiki = as.integer(
      !is.na(lncrnawiki_annotation_status) &
        lncrnawiki_annotation_status != "no lncRNAwiki annotation"
    ),
    has_lit = as.integer(
      !is.na(n_literature_rows) & n_literature_rows > 0
    ),
    has_go = as.integer(
      !is.na(has_go_annotation) & has_go_annotation
    )
  ) %>% 
  mutate(
    has_any_ann = as.integer(has_lncRNAwiki == 1 | has_lit == 1 | has_go == 1)
  ) %>% 
  dplyr::select(
    gene_stable_id, gene_name, matched_symbols, matched_synonyms,
    chrom, gene_start, gene_end, strand, 
    n_overlapping_hervs, overlapping_subfamilies, 
    exonic_herv_percent, gene_body_herv_percent,
    diseases, go_bp_terms, functions_detected,
    functional_mechanisms, biological_processes, 
    molecular_functions, pathways, 
    has_lncRNAwiki, has_lit, has_go, has_any_ann,
    herv_coverage_category   # RELEASE (D029)
  ) %>% 
  # RELEASE (D029): the canonical order-37 category is carried; no cut point is re-applied here
  mutate(
    herv_contribution = factor(
      herv_coverage_category,
      levels = c("minor", "substantial", "dominant", "predominant"),
      ordered = TRUE
    )
  ) %>% 
  mutate(
    across(
      where(is.character),
      ~ na_if(.x, "no lncRNAwiki annotation")
    )
  ) %>% 
  dplyr::select(
    gene_stable_id, gene_name, matched_symbols, matched_synonyms,
    chrom, gene_start, gene_end, strand, 
    n_overlapping_hervs, overlapping_subfamilies, 
    exonic_herv_percent, gene_body_herv_percent, herv_contribution,
    has_lncRNAwiki, has_lit, has_go, has_any_ann,
    diseases, go_bp_terms, functions_detected,
    functional_mechanisms, biological_processes, 
    molecular_functions, pathways
  ) %>% 
  arrange(desc(gene_body_herv_percent))


# ============================================================
# Panel A — Functional annotation coverage of HERV-associated lncRNAs
# ============================================================


panelA_df <- herv_simplified %>%
  distinct(gene_stable_id, .keep_all = TRUE) %>%
  mutate(
    evidence_class = case_when(
      has_lncRNAwiki == 0 & has_lit == 0 & has_go == 0 ~ "No annotation",
      has_lncRNAwiki == 1 & has_lit == 0 & has_go == 0 ~ "LncRNAWiki only",
      has_lncRNAwiki == 0 & has_lit == 1 & has_go == 0 ~ "Literature only",
      has_lncRNAwiki == 0 & has_lit == 0 & has_go == 1 ~ "GO only",
      TRUE ~ "Multiple sources"
    ),
    evidence_class = factor(
      evidence_class,
      levels = c(
        "No annotation",
        "GO only",
        "LncRNAWiki only",
        "Literature only",
        "Multiple sources"
      )
    )
  ) %>%
  count(evidence_class, name = "n") %>%
  mutate(
    fraction = n / sum(n),
    label = paste0(comma(n), "\n(", percent(fraction, accuracy = 0.1), ")")
  )

# ============================================================
# Panel B: broad functional categories among annotated HERV-lncRNAs
# ============================================================

panelB_input <- herv_simplified %>%
  filter(has_any_ann == 1) %>%
  mutate(
    functions_combined = pmap_chr(
      list(functions_detected,
           functional_mechanisms,
           biological_processes,
           molecular_functions,
           pathways,
           diseases),
      ~ paste(na.omit(c(...)), collapse = "; ")
    ),
    functions_combined = str_to_lower(functions_combined)
  ) %>%
  distinct(gene_stable_id, .keep_all = TRUE)

panelB_long <- panelB_input %>%
  mutate(
    category_ceRNA = str_detect(functions_combined, "cerna|competing endogenous rna|mirna sponge|microrna sponge|mirna"),
    category_proliferation = str_detect(functions_combined, "proliferat|growth"),
    category_apoptosis = str_detect(functions_combined, "apoptosis|apoptotic|cell death|anti-apoptotic|pro-apoptotic"),
    category_migration_invasion = str_detect(functions_combined, "migration|invasion|migratory|invasive"),
    category_metastasis = str_detect(functions_combined, "metastasis|metastatic"),
    category_emt = str_detect(functions_combined, "emt|epithelial[- ]mesenchymal transition"),
    category_cell_cycle = str_detect(functions_combined, "cell cycle"),
    category_differentiation = str_detect(functions_combined, "differentiation|stemness|stem cell"),
    category_immune_inflammation = str_detect(functions_combined, "immune|immunity|inflamm|interferon|cytokine"),
    category_signaling = str_detect(functions_combined, "signaling|signal transduction|pathway|akt|erk|mapk|wnt|nf-kb|tgf|pi3k|stat"),
    category_cancer = str_detect(functions_combined, "cancer|tumou?r|oncogen|carcinoma|malignan|metastasis")
  ) %>%
  dplyr::select(
    gene_stable_id,
    starts_with("category_")
  ) %>%
  pivot_longer(
    cols = starts_with("category_"),
    names_to = "category",
    values_to = "present"
  ) %>%
  filter(present) %>%
  mutate(
    category = recode(
      category,
      category_ceRNA = "ceRNA / miRNA sponge",
      category_proliferation = "Proliferation",
      category_apoptosis = "Apoptosis / cell death",
      category_migration_invasion = "Migration / invasion",
      category_metastasis = "Metastasis",
      category_emt = "EMT",
      category_cell_cycle = "Cell cycle",
      category_differentiation = "Differentiation",
      category_immune_inflammation = "Immune / inflammation",
      category_signaling = "Signaling / pathways",
      category_cancer = "Cancer-related"
    )
  ) %>%
  distinct(gene_stable_id, category)


panelB_df <- panelB_long %>%
  count(category, name = "n_genes", sort = TRUE) %>%
  mutate(
    fraction = n_genes / n_distinct(panelB_input$gene_stable_id),
    label = paste0(percent(fraction, accuracy = 0.1))
  )

panelC_df <- herv_simplified %>%
  distinct(gene_stable_id, .keep_all = TRUE) %>%
  group_by(herv_contribution) %>%
  summarise(
    n_total = n(),
    n_annotated = sum(has_any_ann == 1, na.rm = TRUE),
    frac_annotated = n_annotated / n_total,
    .groups = "drop"
  ) %>%
  mutate(
    herv_contribution = factor(
      herv_contribution,
      levels = c("minor", "substantial", "dominant", "predominant")
    ),
    label = paste0(
      percent(frac_annotated, accuracy = 0.1),
      "\n(",
      n_annotated, "/", n_total,
      ")"
    )
  )


# ------------------------------------------------------------------------------
# RELEASE: canonical integrated table, Figure 4 source tables and gene-level membership
# ------------------------------------------------------------------------------
stopifnot(nrow(herv_with_literature) == 4074, nrow(herv_simplified) == 4074)
w <- function(df, name) write.table(df, file = file.path(STAGE, name), quote = FALSE, sep = "\t", row.names = FALSE, fileEncoding = "UTF-8")
w(herv_simplified, "herv_lncrna_simplified.tsv")
# the plot-only 'label' strings (with embedded newlines) are left to 08_figures_supplement
w(dplyr::select(panelA_df, -label), "figure4_panelA_annotation_sources.tsv")
w(dplyr::select(panelB_df, -label), "figure4_panelB_functional_categories.tsv")
w(dplyr::select(panelC_df, -label), "figure4_panelC_annotation_by_coverage_category.tsv")
w(panelB_long, "figure4_panelB_gene_categories.tsv")      # gene x category membership (source of the counts)
membership <- herv_simplified %>%
  dplyr::select(gene_stable_id, gene_name, herv_coverage_category = herv_contribution, has_lncRNAwiki, has_lit, has_go, has_any_ann) %>%
  mutate(herv_coverage_category = as.character(herv_coverage_category)) %>%
  left_join(panelA_class <- herv_simplified %>%
              distinct(gene_stable_id, .keep_all = TRUE) %>%
              transmute(gene_stable_id,
                        evidence_class = case_when(
                          has_lncRNAwiki == 0 & has_lit == 0 & has_go == 0 ~ "No annotation",
                          has_lncRNAwiki == 1 & has_lit == 0 & has_go == 0 ~ "LncRNAWiki only",
                          has_lncRNAwiki == 0 & has_lit == 1 & has_go == 0 ~ "Literature only",
                          has_lncRNAwiki == 0 & has_lit == 0 & has_go == 1 ~ "GO only",
                          TRUE ~ "Multiple sources")),
            by = "gene_stable_id") %>%
  left_join(panelB_long %>% filter(category == "ceRNA / miRNA sponge") %>% transmute(gene_stable_id, cerna_related = 1L),
            by = "gene_stable_id") %>%
  mutate(cerna_related = ifelse(is.na(cerna_related), 0L, cerna_related)) %>%
  arrange(gene_stable_id)
w(membership, "annotation_membership.tsv")
message("[ok] annotated ", sum(membership$has_any_ann), " / ", nrow(membership), "; ceRNA-related ", sum(membership$cerna_related),
        "; predominant annotated ", sum(membership$has_any_ann[membership$herv_coverage_category == "predominant"]), " / ",
        sum(membership$herv_coverage_category == "predominant"))
Sys.chmod(list.files(STAGE, full.names = TRUE), "0644")
if (dir.exists(FINAL_OUT)) unlink(FINAL_OUT, recursive = TRUE)
if (!file.rename(STAGE, FINAL_OUT)) stop("could not publish ", FINAL_OUT)
Sys.chmod(FINAL_OUT, "0755")
message("[OK] wrote ", FINAL_OUT)
