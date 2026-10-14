#!/usr/bin/env Rscript
# ==============================================================================
# 05_lncrna step 4 (frozen global order 37, effective order 38), LOCAL:
# annotate the HERV-associated lncRNA genes with LncRNAWiki and GO.
#
# REFACTOR_SCRIPT of B_transcript_context/04_annotate_lncrna_overlap_tables.R (sha256 afd6fb2f...).
# The matching and collapsing logic (section 3, annotate_lncrna_table) and the helper
# functions are copied VERBATIM. Release changes, each marked "RELEASE":
#   * input  = the accepted order-35 coverage table (D027) joined 1:1 with the order-37 category
#              table (D029); the category columns are carried, and no cut point exists here
#   * LncRNAWiki = frozen SRC14a lncrnawiki.json (sha256-checked), the historical file
#   * GO     = the frozen, recovered SRC17 snapshot (D030) joined on the stripped gene id, in
#              place of the live Ensembl BioMart query (biomaRt is not loaded)
#   * removed: the Figure 3 class_counts block (obsolete; categories belong to effective order 37)
#              and the res_cov.RDS serialisation (consumers read the TSVs)
#   * high-HERV annotation counts select herv_coverage_category == "predominant" instead of
#     re-applying 'gene_body_herv_fraction > 0.80' (the same 144 genes)
#   * explicit arguments, checksum guards, staged output published atomically
# Matching (historical, unchanged): gene id, stripped gene id, transcript id and stripped
# transcript id (only if the input has a transcript column; the coverage table has none),
# normalised symbol, normalised synonym; union, distinct on (gene, entry, disease,
# biological_process, functional_mechanism, match_type); collapsed per gene.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(jsonlite)
  library(digest)
})

# ------------------------------------------------------------------------------
# RELEASE: arguments and input guards
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("missing argument ", flag)
  args[[i + 1]]
}
IN_COVERAGE   <- get_arg("--coverage")
IN_CATEGORIES <- get_arg("--categories")
IN_WIKI_JSON  <- get_arg("--lncrnawiki-json")
IN_GO         <- get_arg("--go-snapshot")
SNAPSHOT_SUMS <- get_arg("--snapshot-checksums")
FINAL_OUT     <- get_arg("--outdir")

sha <- function(p) digest(file = p, algo = "sha256")
verify_manifest <- function(path) {
  m <- read.table(file.path(dirname(path), "checksums.sha256"), col.names = c("sha", "file"), colClasses = "character")
  want <- m$sha[m$file %in% c(basename(path), paste0("./", basename(path)))]
  if (length(want) != 1 || sha(path) != want) stop(path, " does not match its accepted manifest")
  message("[ok] ", basename(path), " sha256 ", want, " (manifest)")
}
verify_snapshot <- function(path, source_id) {
  m <- read.delim(SNAPSHOT_SUMS, colClasses = "character")
  want <- m$sha256[m$source_id == source_id & m$filename == basename(path)]
  if (length(want) != 1 || sha(path) != want) stop(path, " is not the registered ", source_id, " snapshot")
  message("[ok] ", source_id, " ", basename(path), " sha256 ", want)
}
verify_manifest(IN_COVERAGE)
verify_manifest(IN_CATEGORIES)
verify_snapshot(IN_WIKI_JSON, "SRC14a")
verify_snapshot(IN_GO, "SRC17")


# ------------------------------------------------------------------------------
# Helper functions, VERBATIM from utils_transcript_context.R (sha256 a5e81361...)
# ------------------------------------------------------------------------------
collapse_chr <- function(x, sep = ";") {
  x <- as.character(x)
  x <- x[!is.na(x) & x != "" & x != "NA"]
  x <- unique(trimws(x))
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = sep)
}

collapse_unique_keep_na <- function(x, sep = ";") {
  x <- as.character(x)
  x <- x[!is.na(x) & x != "" & x != "NA"]
  x <- unique(trimws(x))
  x <- x[x != ""]
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = sep)
}

strip_version <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  sub("\\.[0-9]+$", "", x)
}

na_to_empty <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

norm_symbol <- function(x) {
  x <- na_to_empty(x)
  x <- str_trim(x)
  x <- str_replace_all(x, fixed('"'), "")
  x <- str_replace_all(x, "\\s+", "")
  toupper(x)
}

split_synonyms <- function(x) {
  x <- na_to_empty(x)
  str_split(x, pattern = ";", simplify = FALSE)
}

# ------------------------------------------------------------------------------
# RELEASE: canonical input = order-35 coverage (17 columns) + order-37 category (1:1)
# ------------------------------------------------------------------------------
coverage <- fread(IN_COVERAGE, sep = "\t", header = TRUE, quote = "", data.table = FALSE)
categories <- fread(IN_CATEGORIES, sep = "\t", header = TRUE, quote = "", data.table = FALSE, encoding = "UTF-8")
if (anyDuplicated(coverage$gene_stable_id) || anyDuplicated(categories$gene_stable_id) ||
    !setequal(coverage$gene_stable_id, categories$gene_stable_id)) stop("coverage and category tables must share one row per gene")
cat_cols <- c("herv_coverage_category", "herv_coverage_category_label", "herv_coverage_category_rank")
chk <- categories[match(coverage$gene_stable_id, categories$gene_stable_id), ]
if (!identical(as.numeric(chk$herv_bp_in_gene_body), as.numeric(coverage$herv_bp_in_gene_body)) ||
    !identical(as.numeric(chk$gene_body_bp), as.numeric(coverage$gene_body_bp))) stop("category table does not match the coverage table")
annot_input <- cbind(coverage, chk[, cat_cols])

GO_SNAPSHOT <- fread(IN_GO, sep = "\t", quote = "", na.strings = "NA", data.table = FALSE) %>% dplyr::select(-gene_stable_id)
if (anyDuplicated(GO_SNAPSHOT$gene_stable_id_stripped)) stop("duplicated gene ids in the GO snapshot")
missing_go <- setdiff(strip_version(coverage$gene_stable_id), GO_SNAPSHOT$gene_stable_id_stripped)
if (length(missing_go) > 0) stop(length(missing_go), " universe genes were not in the historical GO query")

# RELEASE: staged output, created only after every input check has passed
parent <- dirname(normalizePath(FINAL_OUT, mustWork = FALSE))
dir.create(parent, recursive = TRUE, showWarnings = FALSE)
OUT_DIR <- tempfile(pattern = ".stage_lncrna_annotation.", tmpdir = parent)
dir.create(OUT_DIR)

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
annotate_lncrna_table <- function(df,
                                  out_prefix,
                                  keep_input_cols = NULL,
                                  placeholder = "no lncRNAwiki annotation",
                                  add_go_annotations = TRUE) {
  cat("\n========================================\n")
  cat("Processing: order-35 coverage + order-37 categories\n")
  cat("========================================\n")
  
  # --------------------------------------------------------------------------
  # Read input table
  # --------------------------------------------------------------------------
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
    # RELEASE (D030): the recovered SRC17 snapshot replaces the live Ensembl BioMart query
    cat("Joining GO annotations from the frozen SRC17 snapshot...\n")
    go_summary <- GO_SNAPSHOT
    
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
  "overlapping_subfamilies",
  # RELEASE (D029): the canonical category from effective order 37, carried, never re-derived
  "herv_coverage_category", "herv_coverage_category_label", "herv_coverage_category_rank"
)

# ------------------------------------------------------------------------------
# 5) Run the annotation workflow
# ------------------------------------------------------------------------------
res_cov <- annotate_lncrna_table(
  df = annot_input,
  out_prefix = "lncRNA_herv_coverage",
  keep_input_cols = keep_cols_coverage,
  placeholder = "no lncRNAwiki annotation"
)

# ------------------------------------------------------------------------------
# 6) Summary tables (RELEASE: the class_counts block is removed; categories are order 37's)
# ------------------------------------------------------------------------------
# Number of annotated vs non-annotated genes
annotation_counts <- res_cov$full_annotated %>%
  distinct(gene_stable_id, .keep_all = TRUE) %>%
  summarise(
    annotated = sum(n_lncrnawiki_rows > 0, na.rm = TRUE),
    unannotated = sum(n_lncrnawiki_rows == 0, na.rm = TRUE),
    total = n()
  )

# Annotation status among the predominant (>80%) genes: RELEASE, from the canonical category
high_herv_annot <- res_cov$full_annotated %>%
  distinct(gene_stable_id, .keep_all = TRUE) %>%
  filter(herv_coverage_category == "predominant") %>%
  summarise(
    total_high_herv = n(),
    annotated_high_herv = sum(n_lncrnawiki_rows > 0, na.rm = TRUE),
    unannotated_high_herv = sum(n_lncrnawiki_rows == 0, na.rm = TRUE)
  )

fwrite(annotation_counts, file.path(OUT_DIR, "lncRNA_herv_coverage.annotation_counts.tsv"), sep = "\t")
fwrite(high_herv_annot, file.path(OUT_DIR, "lncRNA_herv_coverage.high_herv_annotation_counts.tsv"), sep = "\t")
print(annotation_counts)
print(high_herv_annot)

# RELEASE: atomic publication
Sys.chmod(list.files(OUT_DIR, full.names = TRUE), "0644")
if (dir.exists(FINAL_OUT)) unlink(FINAL_OUT, recursive = TRUE)
if (!file.rename(OUT_DIR, FINAL_OUT)) stop("could not publish ", FINAL_OUT)
Sys.chmod(FINAL_OUT, "0755")
message("[OK] wrote ", FINAL_OUT)
