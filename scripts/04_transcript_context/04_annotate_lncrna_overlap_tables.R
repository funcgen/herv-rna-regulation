#!/usr/bin/env Rscript

# ==============================================================================
# 03_annotate_lncrna_overlap_tables.R
# ==============================================================================
# Purpose
#   Annotate lncRNA/HERV overlap tables with external functional information.
#
# Main questions addressed
#   1. Which lncRNAs overlapping HERVs have existing annotations in LncRNAWiki?
#   2. Can these lncRNAs be matched by gene ID, transcript ID, symbol, or synonym?
#   3. Which GO terms are associated with the overlapped lncRNAs?
#   4. How many high-HERV-content lncRNAs already have literature/database support?
#
# Inputs
#   - ../../../results/lncRNA/lncRNA_herv_coverage.tsv
#   - ../../../../lncRNAs/lncrnawiki.json
#
# Outputs
#   Written to:
#     ../../../results/lncRNA/
#
#   Main annotation outputs:
#     - lncRNA_herv_coverage.annotated_with_lncrnawiki.full.tsv
#     - lncRNA_herv_coverage.annotated_with_lncrnawiki.entry_level.tsv
#     - lncRNA_herv_coverage.annotated_with_lncrnawiki.gene_level.tsv
#     - lncRNA_herv_coverage.no_lncrnawiki_match.tsv
#
#   Summary / QC outputs:
#     - lncRNA_herv_coverage.class_counts.tsv
#     - lncRNA_herv_coverage.annotation_counts.tsv
#     - lncRNA_herv_coverage.high_herv_annotation_counts.tsv
#
# Dependencies
#   - utils_transcript_context.R
#
# Notes
#   - This is the third script to run in the transcript-context module.
#   - It expects the lncRNA/HERV overlap or coverage table generated upstream.
#   - Matching to LncRNAWiki is intentionally permissive and uses multiple
#     identifier strategies:
#       * gene ID
#       * stripped gene ID
#       * transcript ID
#       * stripped transcript ID
#       * symbol
#       * synonym
#   - GO annotation is queried live from Ensembl BioMart.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(jsonlite)
  library(biomaRt)
})

# ------------------------------------------------------------------------------
# Load shared helper functions
# ------------------------------------------------------------------------------
# Provides:
#   - norm_symbol()
#   - strip_version()
#   - split_synonyms()
#   - collapse_chr()
#   - collapse_unique_keep_na()
source("utils_transcript_context.R")

# ------------------------------------------------------------------------------
# Input / output paths
# ------------------------------------------------------------------------------
IN_COVERAGE  <- "../../../results/lncRNA/lncRNA_herv_coverage.tsv"
IN_WIKI_JSON <- "../../../../lncRNAs/lncrnawiki.json"
OUT_DIR      <- "../../../results/lncRNA"

# Ensure output directory exists
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1) Query GO annotations from Ensembl BioMart
# ------------------------------------------------------------------------------
# Given a set of Ensembl gene IDs (without version suffixes), retrieve:
#   - GO IDs
#   - GO term names
#   - GO namespaces (BP / MF / CC)
#
# Notes:
#   - Queries are chunked to avoid overly large requests.
#   - Empty input returns an empty but well-formed table.
query_go_annotations_biomart <- function(gene_ids_stripped, dataset = "hsapiens_gene_ensembl") {
  gene_ids_stripped <- unique(as.character(gene_ids_stripped))
  gene_ids_stripped <- gene_ids_stripped[!is.na(gene_ids_stripped) & gene_ids_stripped != ""]

  if (length(gene_ids_stripped) == 0) {
    return(data.frame(
      ensembl_gene_id = character(),
      go_id = character(),
      name_1006 = character(),
      namespace_1003 = character(),
      stringsAsFactors = FALSE
    ))
  }

  message("Connecting to Ensembl BioMart...")
  mart <- useEnsembl(biomart = "genes", dataset = dataset)
  print(mart)
  
  chunk_size <- 500
  chunks <- split(gene_ids_stripped, ceiling(seq_along(gene_ids_stripped) / chunk_size))

  res_list <- lapply(seq_along(chunks), function(i) {
    message("Querying GO chunk ", i, "/", length(chunks), " ...")
    getBM(
      attributes = c("ensembl_gene_id", "go_id", "name_1006", "namespace_1003"),
      filters = "ensembl_gene_id",
      values = chunks[[i]],
      mart = mart
    )
  })
  
  res <- bind_rows(res_list)
  
  if (!all(c("ensembl_gene_id", "go_id", "name_1006", "namespace_1003") %in% colnames(res))) {
    stop("Missing expected columns. Found: ", paste(colnames(res), collapse = ", "))
  }
  print(head(res))
  res <- res %>%
    # mutate(
    #   ensembl_gene_id = as.character(.data[["ensembl_gene_id"]]),
    #   go_id = na_if(as.character(.data[["go_id"]]), ""),
    #   name_1006 = na_if(as.character(.data[["name_1006"]]), ""),
    #   namespace_1003 = na_if(as.character(.data[["namespace_1003"]]), "")
    # ) %>%
    distinct()
  
  print(head(res))
  
  return(res)
  # bind_rows(res_list) %>%
  #   mutate(
  #     ensembl_gene_id = as.character(ensembl_gene_id),
  #     go_id = na_if(as.character(go_id), ""),
  #     name_1006 = na_if(as.character(name_1006), ""),
  #     namespace_1003 = na_if(as.character(namespace_1003), "")
  #   ) %>%
  #   distinct()
}

# ------------------------------------------------------------------------------
# 2) Collapse raw GO annotations to one row per gene
# ------------------------------------------------------------------------------
# Input:
#   GO annotation table returned by query_go_annotations_biomart()
#
# Output:
#   One row per stripped Ensembl gene ID, with collapsed GO summaries:
#     - all GO IDs
#     - all term names
#     - all namespaces
#     - BP-only / MF-only / CC-only term summaries
#     - number of GO terms
#     - Boolean flag for presence/absence of GO annotations
summarise_go_annotations <- function(go_df) {
  if (nrow(go_df) == 0) {
    return(data.frame(
      gene_stable_id_stripped = character(),
      go_ids = character(),
      go_terms = character(),
      go_namespaces = character(),
      go_bp_terms = character(),
      go_mf_terms = character(),
      go_cc_terms = character(),
      n_go_terms = integer(),
      has_go_annotation = logical(),
      stringsAsFactors = FALSE
    ))
  }
  
  go_df %>%
    dplyr::rename(gene_stable_id_stripped = ensembl_gene_id) %>%
    group_by(gene_stable_id_stripped) %>%
    summarise(
      go_ids = collapse_unique_keep_na(go_id),
      go_terms = collapse_unique_keep_na(name_1006),
      go_namespaces = collapse_unique_keep_na(namespace_1003),
      go_bp_terms = collapse_unique_keep_na(name_1006[namespace_1003 == "biological_process"]),
      go_mf_terms = collapse_unique_keep_na(name_1006[namespace_1003 == "molecular_function"]),
      go_cc_terms = collapse_unique_keep_na(name_1006[namespace_1003 == "cellular_component"]),
      n_go_terms = dplyr::n_distinct(go_id[!is.na(go_id) & go_id != ""]),
      has_go_annotation = n_go_terms > 0,
      .groups = "drop"
    )
}

# ------------------------------------------------------------------------------
# 3) Main annotation function
# ------------------------------------------------------------------------------
# Annotate a lncRNA table by matching entries to LncRNAWiki and optionally
# enriching them with GO annotations.
#
# Matching strategies:
#   - full gene ID
#   - stripped gene ID
#   - full transcript ID
#   - stripped transcript ID
#   - gene symbol
#   - gene synonym
#
# Outputs written by this function:
#   - full appended table
#   - entry-level matches
#   - gene-level annotated subset
#   - no-match subset
annotate_lncrna_table <- function(in_file,
                                  out_prefix,
                                  keep_input_cols = NULL,
                                  placeholder = "no lncRNAwiki annotation",
                                  add_go_annotations = TRUE) {
  cat("\n========================================\n")
  cat("Processing:", in_file, "\n")
  cat("========================================\n")
  
  # --------------------------------------------------------------------------
  # Read input table
  # --------------------------------------------------------------------------
  df <- fread(in_file, sep = "\t", header = TRUE, quote = "", data.table = FALSE)
  cat("Input rows:", nrow(df), "\n")
  
  required_cols <- c("gene_stable_id", "gene_name", "gene_type")
  missing_required <- setdiff(required_cols, colnames(df))
  if (length(missing_required) > 0) {
    stop("Missing required columns: ", paste(missing_required, collapse = ", "))
  }
  
  # Detect transcript identifier column if present
  transcript_col <- intersect(c("transcript_stable_id", "transcript_id"), colnames(df))
  transcript_col <- if (length(transcript_col) > 0) transcript_col[1] else NULL
  
  # Standardize identifiers for matching
  df2 <- df %>%
    mutate(
      gene_name = as.character(gene_name),
      gene_name_norm = norm_symbol(gene_name),
      gene_stable_id = as.character(gene_stable_id),
      gene_stable_id_stripped = strip_version(gene_stable_id)
    )
  
  if (!is.null(transcript_col)) {
    df2 <- df2 %>%
      mutate(
        transcript_input = as.character(.data[[transcript_col]]),
        transcript_input_stripped = strip_version(.data[[transcript_col]])
      )
  }
  
  # --------------------------------------------------------------------------
  # Read and prepare LncRNAWiki
  # --------------------------------------------------------------------------
  wiki_json <- fromJSON(IN_WIKI_JSON)
  
  wiki <- as.data.frame(wiki_json$geneInfo, stringsAsFactors = FALSE) %>%
    dplyr::select(-any_of(c("seq", "descriptions")))
  
  wiki2 <- wiki %>%
    mutate(
      wiki_symbol = as.character(symbol),
      wiki_symbol_norm = norm_symbol(symbol),
      synonyms = as.character(synonyms),
      wiki_gene_id = if ("gene_id" %in% colnames(.)) as.character(gene_id) else NA_character_,
      wiki_gene_id_stripped = if ("gene_id" %in% colnames(.)) strip_version(gene_id) else NA_character_,
      wiki_transcript_id = if ("transcript_id" %in% colnames(.)) as.character(transcript_id) else NA_character_,
      wiki_transcript_id_stripped = if ("transcript_id" %in% colnames(.)) strip_version(transcript_id) else NA_character_
    )
  
  # Expand synonyms to long format for synonym-based matching
  wiki_syn_long <- wiki2 %>%
    dplyr::select(-wiki_symbol_norm) %>%
    mutate(syn_list = split_synonyms(synonyms)) %>%
    unnest_longer(syn_list, values_to = "synonym_raw") %>%
    mutate(synonym_norm = norm_symbol(synonym_raw)) %>%
    filter(!is.na(synonym_norm), synonym_norm != "")
  
  # Small helper to safely skip empty join results in bind_rows()
  join_if_any <- function(tbl) {
    if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
    tbl
  }
  
  # --------------------------------------------------------------------------
  # Perform all matching strategies
  # --------------------------------------------------------------------------
  geneid_full_matches <- df2 %>%
    inner_join(
      wiki2 %>% filter(!is.na(wiki_gene_id), wiki_gene_id != ""),
      by = c("gene_stable_id" = "wiki_gene_id"),
      relationship = "many-to-many"
    ) %>%
    mutate(match_type = "gene_id")
  
  geneid_stripped_matches <- df2 %>%
    inner_join(
      wiki2 %>% filter(!is.na(wiki_gene_id_stripped), wiki_gene_id_stripped != ""),
      by = c("gene_stable_id_stripped" = "wiki_gene_id_stripped"),
      relationship = "many-to-many"
    ) %>%
    mutate(match_type = "gene_id_stripped")
  
  direct_matches <- df2 %>%
    inner_join(
      wiki2 %>% filter(!is.na(wiki_symbol_norm), wiki_symbol_norm != ""),
      by = c("gene_name_norm" = "wiki_symbol_norm"),
      relationship = "many-to-many"
    ) %>%
    mutate(match_type = "symbol")
  
  synonym_matches <- df2 %>%
    inner_join(
      wiki_syn_long %>% filter(!is.na(synonym_norm), synonym_norm != ""),
      by = c("gene_name_norm" = "synonym_norm"),
      relationship = "many-to-many"
    ) %>%
    mutate(match_type = "synonym")
  
  transcript_full_matches <- NULL
  transcript_stripped_matches <- NULL
  
  if (!is.null(transcript_col)) {
    transcript_full_matches <- df2 %>%
      inner_join(
        wiki2 %>% filter(!is.na(wiki_transcript_id), wiki_transcript_id != ""),
        by = c("transcript_input" = "wiki_transcript_id"),
        relationship = "many-to-many"
      ) %>%
      mutate(match_type = "transcript_id")
    
    transcript_stripped_matches <- df2 %>%
      inner_join(
        wiki2 %>% filter(!is.na(wiki_transcript_id_stripped), wiki_transcript_id_stripped != ""),
        by = c("transcript_input_stripped" = "wiki_transcript_id_stripped"),
        relationship = "many-to-many"
      ) %>%
      mutate(match_type = "transcript_id_stripped")
  }
  
  # --------------------------------------------------------------------------
  # Combine all matches into one table
  # --------------------------------------------------------------------------
  all_matches <- bind_rows(
    join_if_any(geneid_full_matches),
    join_if_any(geneid_stripped_matches),
    join_if_any(transcript_full_matches),
    join_if_any(transcript_stripped_matches),
    join_if_any(direct_matches),
    join_if_any(synonym_matches)
  ) %>%
    distinct(
      gene_stable_id, gene_name, entry_id, disease,
      biological_process, functional_mechanism, match_type,
      .keep_all = TRUE
    )
  
  cat("Matched rows:", nrow(all_matches), "\n")
  cat("Matched genes:", all_matches %>% distinct(gene_stable_id) %>% nrow(), "\n")
  cat("Matches by type:\n")
  print(all_matches %>% count(match_type, sort = TRUE))
  
  # Keep only user-requested input columns in downstream outputs
  if (is.null(keep_input_cols)) keep_input_cols <- colnames(df)
  keep_input_cols <- intersect(keep_input_cols, colnames(df2))
  
  # --------------------------------------------------------------------------
  # Collapse matches to one summary row per gene
  # --------------------------------------------------------------------------
  gene_summary <- all_matches %>%
    group_by(gene_stable_id) %>%
    summarise(
      n_lncrnawiki_rows = n(),
      match_types = collapse_chr(match_type),
      matched_symbols = collapse_chr(symbol),
      matched_gene_ids = collapse_chr(gene_id),
      matched_transcript_ids = collapse_chr(transcript_id),
      matched_synonyms = collapse_chr(synonyms),
      functional_mechanisms = collapse_chr(functional_mechanism),
      biological_processes = collapse_chr(biological_process),
      molecular_functions = collapse_chr(molecular_function),
      pathways = collapse_chr(pathway),
      diseases = collapse_chr(disease),
      biological_contexts = collapse_chr(biological_context),
      context_details = collapse_chr(context_detail),
      targets = collapse_chr(target),
      target_types = collapse_chr(target_type),
      target_interactions = collapse_chr(target_interaction),
      target_effects = collapse_chr(target_effect),
      regulators = collapse_chr(regulator),
      regulator_types = collapse_chr(regulator_type),
      regulator_effects = collapse_chr(regulator_effect),
      regulator_interactions = collapse_chr(regulator_interaction),
      expression_patterns = collapse_chr(expression),
      expression_details = collapse_chr(expression_detail),
      experiment_methods = collapse_chr(experiment_method),
      target_experiments = collapse_chr(target_experiment),
      pmids = collapse_chr(pmid),
      entry_ids = collapse_chr(entry_id),
      .groups = "drop"
    )
  
  # Append annotation summary back to the original input table
  full_annotated <- df2 %>%
    left_join(gene_summary, by = "gene_stable_id") %>%
    dplyr::select(-any_of(c("gene_name_norm", "transcript_input", "transcript_input_stripped")))
  
  # --------------------------------------------------------------------------
  # Optional GO annotation
  # --------------------------------------------------------------------------
  if (add_go_annotations) {
    cat("Querying GO annotations from Ensembl BioMart...\n")
    
    go_raw <- query_go_annotations_biomart(unique(full_annotated$gene_stable_id_stripped))
    go_summary <- summarise_go_annotations(go_raw)
    
    full_annotated <- full_annotated %>%
      left_join(go_summary, by = "gene_stable_id_stripped") %>%
      mutate(
        n_go_terms = ifelse(is.na(n_go_terms), 0L, as.integer(n_go_terms)),
        has_go_annotation = ifelse(is.na(has_go_annotation), FALSE, has_go_annotation)
      )
  }
  
  # Clean temporary columns and add annotation-status label
  full_annotated <- full_annotated %>%
    dplyr::select(-any_of("gene_stable_id_stripped")) %>%
    mutate(
      n_lncrnawiki_rows = ifelse(is.na(n_lncrnawiki_rows), 0L, as.integer(n_lncrnawiki_rows)),
      lncrnawiki_annotation_status = ifelse(
        n_lncrnawiki_rows > 0,
        "matched",
        "no lncRNAwiki annotation"
      )
    )
  
  # For genes with no match, fill selected annotation columns with a placeholder
  annotation_cols <- c(
    "match_types", "matched_symbols", "matched_gene_ids", "matched_transcript_ids", "matched_synonyms",
    "functional_mechanisms", "biological_processes", "molecular_functions", "pathways", "diseases",
    "biological_contexts", "context_details", "targets", "target_types", "target_interactions",
    "target_effects", "regulators", "regulator_types", "regulator_effects", "regulator_interactions",
    "expression_patterns", "expression_details", "experiment_methods", "target_experiments", "pmids", "entry_ids"
  )
  
  for (cc in intersect(annotation_cols, colnames(full_annotated))) {
    idx_nomatch <- full_annotated$n_lncrnawiki_rows == 0 &
      (is.na(full_annotated[[cc]]) | full_annotated[[cc]] == "")
    full_annotated[[cc]][idx_nomatch] <- placeholder
  }
  
  # Split full table into annotated and non-annotated subsets
  gene_level_annotated <- full_annotated %>% filter(n_lncrnawiki_rows > 0)
  no_match <- full_annotated %>% filter(n_lncrnawiki_rows == 0)
  
  # --------------------------------------------------------------------------
  # Prepare output file paths
  # --------------------------------------------------------------------------
  out_full    <- file.path(OUT_DIR, paste0(out_prefix, ".annotated_with_lncrnawiki.full.tsv"))
  out_entry   <- file.path(OUT_DIR, paste0(out_prefix, ".annotated_with_lncrnawiki.entry_level.tsv"))
  out_gene    <- file.path(OUT_DIR, paste0(out_prefix, ".annotated_with_lncrnawiki.gene_level.tsv"))
  out_nomatch <- file.path(OUT_DIR, paste0(out_prefix, ".no_lncrnawiki_match.tsv"))
  
  # Reorder entry-level output to prioritize important columns first
  preferred_entry_cols <- c(
    keep_input_cols, "match_type", "symbol", "gene_id", "transcript_id", "synonyms", "entry_id",
    "functional_mechanism", "biological_process", "disease", "biological_context", "context_detail",
    "molecular_function", "pathway", "target", "target_type", "target_interaction", "target_effect",
    "regulator", "regulator_type", "regulator_effect", "regulator_interaction", "expression",
    "expression_detail", "experiment_method", "target_experiment", "pmid", "record_created", "last_modified"
  )
  
  all_matches_out <- all_matches[, c(
    intersect(preferred_entry_cols, colnames(all_matches)),
    setdiff(colnames(all_matches), preferred_entry_cols)
  )]
  
  # --------------------------------------------------------------------------
  # Write outputs
  # --------------------------------------------------------------------------
  fwrite(full_annotated, out_full, sep = "\t", quote = FALSE, na = "NA")
  fwrite(all_matches_out, out_entry, sep = "\t", quote = FALSE, na = "NA")
  fwrite(gene_level_annotated, out_gene, sep = "\t", quote = FALSE, na = "NA")
  fwrite(no_match, out_nomatch, sep = "\t", quote = FALSE, na = "NA")
  
  invisible(list(
    input = df2,
    matches = all_matches,
    gene_summary = gene_summary,
    full_annotated = full_annotated,
    gene_level_annotated = gene_level_annotated,
    no_match = no_match
  ))
}

# ------------------------------------------------------------------------------
# 4) Define columns to retain from the input coverage table
# ------------------------------------------------------------------------------
keep_cols_coverage <- c(
  "gene_stable_id", "gene_name", "gene_type", "chrom", "strand", "gene_start", "gene_end",
  "gene_body_bp", "exonic_bp", "n_overlapping_hervs", "herv_bp_in_gene_body", "herv_bp_in_exons",
  "gene_body_herv_fraction", "gene_body_herv_percent", "exonic_herv_fraction", "exonic_herv_percent",
  "overlapping_subfamilies"
)

# Add transcript column if present in the input table
cov_cols_present <- colnames(fread(IN_COVERAGE, nrows = 0, sep = "\t", header = TRUE))
if ("transcript_stable_id" %in% cov_cols_present) keep_cols_coverage <- c(keep_cols_coverage, "transcript_stable_id")
if ("transcript_id" %in% cov_cols_present) keep_cols_coverage <- c(keep_cols_coverage, "transcript_id")

# ------------------------------------------------------------------------------
# 5) Run the annotation workflow
# ------------------------------------------------------------------------------
res_cov <- annotate_lncrna_table(
  in_file = IN_COVERAGE,
  out_prefix = "lncRNA_herv_coverage",
  keep_input_cols = keep_cols_coverage,
  placeholder = "no lncRNAwiki annotation"
)
saveRDS(res_cov, file = paste0(OUT_DIR, "/res_cov.RDS"))

# ------------------------------------------------------------------------------
# 6) Build additional summary/QC tables
# ------------------------------------------------------------------------------
# Distribution of genes by fraction of gene body covered by HERV sequence
class_counts <- res_cov$full_annotated %>%
  filter(!is.na(gene_body_herv_fraction)) %>%
  distinct(gene_stable_id, .keep_all = TRUE) %>%
  mutate(
    herv_class = case_when(
      gene_body_herv_fraction < 0.10 ~ "<10%",
      gene_body_herv_fraction < 0.50 ~ "10–50%",
      gene_body_herv_fraction < 0.80 ~ "50–80%",
      TRUE                           ~ ">80%"
    ),
    herv_class = factor(herv_class, levels = c("<10%", "10–50%", "50–80%", ">80%"))
  ) %>%
  count(herv_class, name = "n") %>%
  mutate(
    fraction = n / sum(n),
    percent = round(100 * fraction, 1)
  )

# Number of annotated vs non-annotated genes
annotation_counts <- res_cov$full_annotated %>%
  distinct(gene_stable_id, .keep_all = TRUE) %>%
  summarise(
    annotated = sum(n_lncrnawiki_rows > 0, na.rm = TRUE),
    unannotated = sum(n_lncrnawiki_rows == 0, na.rm = TRUE),
    total = n()
  )

# Annotation status specifically among genes with very high HERV coverage
high_herv_annot <- res_cov$full_annotated %>%
  distinct(gene_stable_id, .keep_all = TRUE) %>%
  filter(gene_body_herv_fraction > 0.80) %>%
  summarise(
    total_high_herv = n(),
    annotated_high_herv = sum(n_lncrnawiki_rows > 0, na.rm = TRUE),
    unannotated_high_herv = sum(n_lncrnawiki_rows == 0, na.rm = TRUE)
  )

# ------------------------------------------------------------------------------
# 7) Write summary/QC outputs
# ------------------------------------------------------------------------------
fwrite(class_counts, file.path(OUT_DIR, "lncRNA_herv_coverage.class_counts.tsv"), sep = "\t")
fwrite(annotation_counts, file.path(OUT_DIR, "lncRNA_herv_coverage.annotation_counts.tsv"), sep = "\t")
fwrite(high_herv_annot, file.path(OUT_DIR, "lncRNA_herv_coverage.high_herv_annotation_counts.tsv"), sep = "\t")

# ------------------------------------------------------------------------------
# 8) Print summaries for quick inspection
# ------------------------------------------------------------------------------
print(class_counts)
print(annotation_counts)
print(high_herv_annot)
