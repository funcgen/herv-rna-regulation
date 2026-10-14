#!/usr/bin/env Rscript
# ==============================================================================
# 06_lncrna_function step 3 (frozen global order 40, effective order 41), LOCAL, offline:
# classify the v5 literature candidates as evidence for HERV-overlapping lncRNAs (D032).
#
# D032: the frozen historical source F_literature_mining/03_classify_lncrna.py does not read v5 (it
# builds an unused LncRNADisease atlas: DO_NOT_MIGRATE). The historical v5 classification is the
# literature-matching block of G_literature_integration/exploration_lit_mining.R (sha256 see
# PROVENANCE), which is moved here: sections 1, 3, 4, 5, 6 and 10 are copied VERBATIM.
#   1  clean the v5 candidates: str_trim, remove '"', str_squish, toupper -> Candidate_LncRNA_std
#   2  RELEASE: the HERV lncRNA table is the canonical order-38 full annotated table (TSV) instead of
#      res_cov.RDS$full_annotated
#   3-4 alias table: gene_name, matched_symbols, matched_synonyms, matched_gene_ids,
#      matched_transcript_ids, split on ',' / ';', squished, upper-cased, distinct
#   5  left join v5 x alias table on the upper-cased name (many-to-many); keep matched rows. A
#      candidate is evidence only if it equals an alias of one of the 4,074 HERV lncRNAs: this is the
#      whole false-positive rule (no tier, sentence-type or count threshold)
#   6  per-gene literature summary (n_literature_rows = matched rows, n_pmids = distinct PMIDs, ...)
#   10 sentence-level export
# Outputs (staged, atomic): lit_matches_herv_sentence_level.tsv (historical file name) and
# lit_summary_per_gene.tsv (the section-6 table, historically merged into herv_lncrna_with_literature.tsv).
# Needs R + readr/dplyr/stringr/tidyr/purrr/tibble/data.table/digest only (no spaCy, no network).
# ==============================================================================
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(data.table)
  library(digest)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("missing argument ", flag)
  args[[i + 1]]
}
V5 <- get_arg("--v5")
ANNOT <- get_arg("--annotation")
FINAL_OUT <- get_arg("--outdir")

# RELEASE: inputs must match their accepted manifests
verify_manifest <- function(path) {
  m <- read.table(file.path(dirname(path), "checksums.sha256"), col.names = c("sha", "file"), colClasses = "character")
  want <- m$sha[m$file %in% c(basename(path), paste0("./", basename(path)))]
  if (length(want) != 1 || digest(file = path, algo = "sha256") != want) stop(path, " does not match its accepted manifest")
  message("[ok] ", basename(path), " sha256 ", want, " (manifest)")
}
verify_manifest(V5)
verify_manifest(ANNOT)
if (basename(V5) != "lncrna_functions_extracted_v5.csv") stop("only the canonical v5 function table is accepted")

CFG <- list(lit_csv = V5)

# ============================================================
# 1. Load literature-derived table
# ============================================================

lit_raw <- read_csv(CFG$lit_csv, show_col_types = FALSE)

lit <- lit_raw %>%
  mutate(
    Candidate_LncRNA = str_trim(Candidate_LncRNA),
    Candidate_LncRNA = str_replace_all(Candidate_LncRNA, "\"", ""),
    Candidate_LncRNA = str_squish(Candidate_LncRNA),
    Candidate_LncRNA_std = toupper(Candidate_LncRNA)
  ) %>%
  filter(!is.na(Candidate_LncRNA_std), Candidate_LncRNA_std != "")

# ============================================================
# 2. HERV-overlapping lncRNA table (RELEASE: canonical order-38 TSV, NA restored)
# ============================================================

herv_lncrna <- fread(ANNOT, sep = "\t", quote = "", na.strings = "NA", data.table = FALSE, encoding = "UTF-8")
if (anyDuplicated(herv_lncrna$gene_stable_id)) stop("order-38 table must have one row per gene")

# ============================================================
# 3. Helper to split alias columns
# ============================================================

split_aliases <- function(x) {
  x <- unlist(x, use.names = FALSE)
  x <- x[!is.na(x)]
  
  if (length(x) == 0) return(character(0))
  
  x <- paste(x, collapse = ",")
  
  if (x %in% c("", "no lncRNAwiki annotation")) return(character(0))
  
  x %>%
    str_split(",|;") %>%
    unlist() %>%
    str_trim() %>%
    discard(~ .x %in% c("", "NA", "no lncRNAwiki annotation")) %>%
    unique()
}

# ============================================================
# 4. Build alias table from HERV lncRNA table
# ============================================================

# ============================================================
# 4. Build alias table from HERV lncRNA table (Vectorized)
# ============================================================

# ============================================================
# 4. Build alias table from HERV lncRNA table (Vectorized)
# ============================================================

alias_tbl <- herv_lncrna %>%
  dplyr::select(
    gene_stable_id, 
    gene_name, 
    matched_symbols, 
    matched_synonyms,
    matched_gene_ids, 
    matched_transcript_ids
  ) %>%
  # 1. Create a safe copy of gene_name to keep in the final table
  mutate(base_gene_name = gene_name) %>%
  # 2. Pivot the alias columns into a long format
  pivot_longer(
    cols = c(gene_name, matched_symbols, matched_synonyms, matched_gene_ids, matched_transcript_ids),
    names_to = "alias_type",
    values_to = "alias_raw"
  ) %>%
  # 3. Filter out missing data BEFORE splitting to save processing time
  filter(
    !is.na(alias_raw), 
    alias_raw != "no lncRNAwiki annotation", 
    alias_raw != ""
  ) %>%
  # 4. Split strings containing commas or semicolons into new rows
  separate_rows(alias_raw, sep = "[,;]") %>%
  # 5. Clean up whitespace and standardize formatting
  mutate(
    alias = str_squish(alias_raw),
    alias_std = toupper(alias)
  ) %>%
  # 6. Remove any artifacts created by trailing commas/semicolons
  filter(alias_std != "", alias_std != "NA") %>%
  # 7. Restore the base gene_name column and finalize structure
  dplyr::rename(gene_name = base_gene_name) %>%
  dplyr::select(gene_stable_id, gene_name, alias_type, alias, alias_std) %>%
  distinct()

# Optional QC
print(alias_tbl)
alias_tbl %>% count(alias_type)


alias_tbl %>%
  filter(gene_name %in% c("UCA1", "LINC01671", "CERNA2")) %>%
  arrange(gene_name, alias_type, alias)

# ============================================================
# 5. Match literature candidates to alias table
# ============================================================

lit_matches <- lit %>%
  left_join(
    alias_tbl,
    by = c("Candidate_LncRNA_std" = "alias_std"),
    relationship = "many-to-many"
  )

lit_matches_herv <- lit_matches %>%
  filter(!is.na(gene_stable_id))

# ============================================================
# 6. Summarise literature evidence per HERV-overlapping gene
# ============================================================

collapse_terms <- function(x, sep = "; ") {
  # Flatten + character
  x <- as.character(unlist(x, use.names = FALSE))
  
  # Split FIRST
  x <- str_split(x, ";")
  x <- unlist(x)
  
  # Clean AFTER splitting
  x <- str_trim(x)
  
  # Remove garbage
  x <- x[!is.na(x) & x != "" & x != "NA"]
  
  if (length(x) == 0) return(NA_character_)
  
  x %>%
    unique() %>%
    sort() %>%
    paste(collapse = sep)
}

collapse_unique <- function(x, sep = "; ") {
  x <- x[!is.na(x) & x != ""]
  x <- unique(x)
  if (length(x) == 0) return(NA_character_)
  paste(sort(x), collapse = sep)
}

lit_summary_per_gene <- lit_matches_herv %>%
  group_by(gene_stable_id, gene_name) %>%
  summarise(
    n_literature_rows = n(),
    n_pmids = n_distinct(PMID),
    pmids_literature = collapse_unique(as.character(PMID)),
    
    matched_aliases = collapse_unique(alias),
    alias_types = collapse_unique(alias_type),
    
    functions_detected = collapse_terms(Functions_Detected),
    mechanistic_terms = collapse_terms(Mechanistic_Terms),
    phenotypic_terms = collapse_terms(Phenotypic_Terms),
    associative_terms = collapse_terms(Associative_Terms),
    pathway_terms = collapse_terms(Pathway_Terms),
    regulatory_verbs = collapse_terms(Regulatory_Verbs),
    
    sentence_types = collapse_unique(Sentence_Type),
    evidence_tiers = collapse_unique(Evidence_Tier),
    
    min_candidates_in_sentence = suppressWarnings(min(n_candidates_in_sentence, na.rm = TRUE)),
    max_candidates_in_sentence = suppressWarnings(max(n_candidates_in_sentence, na.rm = TRUE)),
    
    evidence_sentences = paste(head(unique(Evidence_Sentence), 5), collapse = " || "),
    .groups = "drop"
  ) %>%
  mutate(
    min_candidates_in_sentence = ifelse(is.infinite(min_candidates_in_sentence), NA, min_candidates_in_sentence),
    max_candidates_in_sentence = ifelse(is.infinite(max_candidates_in_sentence), NA, max_candidates_in_sentence)
  )

head(lit_summary_per_gene)
nrow(lit_summary_per_gene)

# RELEASE: staged output, created after all inputs are verified
parent <- dirname(normalizePath(FINAL_OUT, mustWork = FALSE))
dir.create(parent, recursive = TRUE, showWarnings = FALSE)
STAGE <- tempfile(pattern = ".stage_literature_classification.", tmpdir = parent)
dir.create(STAGE)

# ============================================================
# 10. Optional: raw sentence-level table for matched HERV lncRNAs
# ============================================================

lit_matches_herv_export <- lit_matches_herv %>%
  dplyr::select(
    gene_stable_id,
    gene_name,
    Candidate_LncRNA,
    alias,
    alias_type,
    PMID,
    Year,
    Evidence_Tier,
    Sentence_Type,
    Mechanistic_Terms,
    Phenotypic_Terms,
    Associative_Terms,
    Pathway_Terms,
    Regulatory_Verbs,
    Functions_Detected,
    n_candidates_in_sentence,
    Evidence_Sentence
  )

write.table(
  lit_matches_herv_export,
  file = file.path(STAGE, "lit_matches_herv_sentence_level.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)


# RELEASE: the section-6 per-gene summary, written as its own table
write.table(
  lit_summary_per_gene,
  file = file.path(STAGE, "lit_summary_per_gene.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
message("[ok] ", nrow(lit_matches_herv_export), " matched sentence-level rows; ", nrow(lit_summary_per_gene), " genes with literature")
Sys.chmod(list.files(STAGE, full.names = TRUE), "0644")
if (dir.exists(FINAL_OUT)) unlink(FINAL_OUT, recursive = TRUE)
if (!file.rename(STAGE, FINAL_OUT)) stop("could not publish ", FINAL_OUT)
Sys.chmod(FINAL_OUT, "0755")
message("[OK] wrote ", FINAL_OUT)
