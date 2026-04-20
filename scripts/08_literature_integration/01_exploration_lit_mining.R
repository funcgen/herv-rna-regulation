# ==============================================================================
# 01_exploration_lit_mining.R
# ==============================================================================

# Purpose
# Integrate literature-derived annotations, LncRNAWiki, and GO information
# to characterize the functional landscape of HERV-associated lncRNAs.
#
# Overview
#   1. Map literature-extracted lncRNA candidates to HERV-overlapping genes
#      using a comprehensive alias matching strategy
#   2. Summarize functional evidence per gene (PMIDs, mechanisms, terms)
#   3. Integrate annotations into a unified lncRNA–HERV table
#   4. Generate simplified annotation layers for downstream analysis
#   5. Visualize:
#        - annotation coverage (Panel A)
#        - functional category distribution (Panel B)
#        - relationship between HERV contribution and annotation (Panel C)
#
# Outputs
#   - Integrated lncRNA–HERV annotation tables
#   - Literature-matched gene summaries
#   - Multi-panel figure summarizing functional annotations
#
# Notes
#   - Alias matching is critical (symbols, synonyms, IDs)
#   - Literature annotations are collapsed at gene level
#   - Functions are aggregated across multiple annotation sources
# ==============================================================================

# ============================================================
# Libraries
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(forcats)
  library(scales)
  library(patchwork)
  library(cowplot)
  library(openxlsx)
})




# ============================================================
# Config / paths
# ============================================================

CFG <- list(
  res_cov_rds = "../../../results/lncRNA/res_cov.RDS",
  lit_csv     = "../../../lncRNAs/lncrna_functions_extracted_v5.csv",
  
  out_dir     = "../../../results/lit_mining/",
  fig_dir     = "../../../results/lit_mining/figure_herv_lncrna_annotations/"
)

dir.create(CFG$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(CFG$fig_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load data
# ============================================================

res_cov <- readRDS(CFG$res_cov_rds)



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
# 2. HERV-overlapping lncRNA table
# ============================================================

herv_lncrna <- res_cov$full_annotated

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

alias_tbl <- pmap_dfr(
  herv_lncrna %>%
    dplyr::select(
      gene_stable_id, gene_name,
      matched_symbols, matched_synonyms,
      matched_gene_ids, matched_transcript_ids
    ),
  function(gene_stable_id, gene_name, matched_symbols, matched_synonyms, matched_gene_ids, matched_transcript_ids) {
    
    tibble(
      gene_stable_id = gene_stable_id,
      gene_name = gene_name,
      alias_type = c(
        rep("gene_name", length(split_aliases(gene_name))),
        rep("matched_symbols", length(split_aliases(matched_symbols))),
        rep("matched_synonyms", length(split_aliases(matched_synonyms))),
        rep("matched_gene_ids", length(split_aliases(matched_gene_ids))),
        rep("matched_transcript_ids", length(split_aliases(matched_transcript_ids)))
      ),
      alias = c(
        split_aliases(gene_name),
        split_aliases(matched_symbols),
        split_aliases(matched_synonyms),
        split_aliases(matched_gene_ids),
        split_aliases(matched_transcript_ids)
      )
    )
  }
) %>%
  mutate(
    alias = str_squish(alias),
    alias_std = toupper(alias)
  ) %>%
  filter(!is.na(alias_std), alias_std != "") %>%
  distinct()

# Optional QC
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
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  
  x %>%
    paste(collapse = ";") %>%
    str_split(";") %>%
    unlist() %>%
    str_trim() %>%
    discard(~ .x == "" || .x == "NA") %>%
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
    
    evidence_sentences = paste(unique(Evidence_Sentence)[1:min(5, dplyr::n())], collapse = " || "),
    .groups = "drop"
  ) %>%
  mutate(
    min_candidates_in_sentence = ifelse(is.infinite(min_candidates_in_sentence), NA, min_candidates_in_sentence),
    max_candidates_in_sentence = ifelse(is.infinite(max_candidates_in_sentence), NA, max_candidates_in_sentence)
  )

head(lit_summary_per_gene)
nrow(lit_summary_per_gene)

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

check_tbl %>%
  count(has_lncrnawiki, has_literature)

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
  file = file.path(CFG$out_dir, "lit_only_inspection.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

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
  file = file.path(CFG$out_dir, "lit_matches_herv_sentence_level.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

# ============================================================
# 11. Integrated final table
# ============================================================

write.table(
  herv_with_literature,
  file = file.path(CFG$out_dir, "herv_lncrna_with_literature.tsv"),
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
    has_lncRNAwiki, has_lit, has_go, has_any_ann
  ) %>% 
  mutate(
    herv_contribution = case_when(
      gene_body_herv_percent <= 10 ~ "minor", 
      gene_body_herv_percent > 10 & gene_body_herv_percent <= 50 ~ "substantial",
      gene_body_herv_percent > 50 & gene_body_herv_percent <= 80 ~ "dominant",
      gene_body_herv_percent > 80 ~ "predominant"
      ),
    herv_contribution = factor(
      herv_contribution, 
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

herv_simplified %>% count(herv_contribution)  
herv_simplified %>% count(has_any_ann)

table(herv_simplified$has_any_ann, herv_simplified$herv_contribution)

herv_simpl_pred <- herv_simplified %>% filter(herv_contribution == "predominant")

herv_simplified_list <- list("lncRNA_hervs" = herv_simplified)
openxlsx::write.xlsx(herv_simplified_list, file = "../results/lncRNA/herv_lncrna_simplified.xlsx")

nrow(herv_simplified)


herv_simplified$functions_detected[!is.na(herv_simplified$functions_detected)]
herv_simplified$functions_detected[!is.na(herv_simplified$functions_detected)]

library(dplyr)
library(purrr)

herv_simplified <- herv_simplified %>%
  mutate(
    functions = pmap_chr(
      list(functions_detected,
           functional_mechanisms,
           biological_processes,
           molecular_functions),
      ~ collapse_terms(paste(na.omit(c(...)), collapse = ";"))
    )
  )

head(herv_simplified$functions[!is.na(herv_simplified$functions)])


library(stringr)

# Convert to lowercase
herv_simplified <- herv_simplified %>%
  mutate(functions = str_to_lower(functions))


# create a table of most common terms
functions_popularity <- herv_simplified %>%
  filter(!is.na(functions)) %>%
  separate_rows(functions, sep = ";") %>%
  mutate(functions = str_trim(functions)) %>%
  count(functions, sort = TRUE)

head(functions_popularity, 20)




# ============================================================
# Panel A — Functional annotation coverage of HERV-associated lncRNAs
# ============================================================

library(dplyr)
library(ggplot2)
library(scales)

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

pA <- ggplot(panelA_df, aes(x = evidence_class, y = n, fill = evidence_class)) +
  geom_col(width = 0.72, color = NA) +
  geom_text(
    aes(label = label),
    vjust = -0.35,
    size = 3.8,
    lineheight = 0.95
  ) +
  scale_fill_manual(values = c(
    "No annotation"    = "#D9D9D9",
    "GO only"          = "#BFD3E6",
    "LncRNAWiki only"  = "#9ECAE1",
    "Literature only"  = "#6BAED6",
    "Multiple sources" = "#3182BD"
  )) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.10)),
    labels = comma
  ) +
  labs(
    x = NULL,
    y = "Number of HERV-associated lncRNAs"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13),
    axis.text.x = element_text(size = 11, angle = 20, hjust = 1),
    axis.text.y = element_text(size = 11),
    plot.margin = margin(8, 8, 8, 8)
  )

pA





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
      category_differentiation = "Differentiation / stemness",
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

pB <- ggplot(panelB_df, aes(x = n_genes, y = fct_reorder(category, n_genes))) +
  geom_col(width = 0.72, fill = "#4C78A8") +
  geom_text(
    aes(
      label = label,
      hjust = ifelse(n_genes > 50, 1.05, -0.1)  # inside if large, outside if small
    ),
    color = ifelse(panelB_df$n_genes > 50, "white", "black"),
    size = 3.8
  ) +
  scale_x_continuous(
    limits = c(0, max(panelB_df$n_genes) * 1.2),
    expand = c(0, 0),
    labels = comma
  ) +
  labs(
    x = "Number of annotated HERV-associated lncRNAs",
    y = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11)
  )




library(dplyr)
library(ggplot2)
library(scales)

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

pC <- ggplot(panelC_df, aes(x = herv_contribution, y = frac_annotated, group = 1)) +
  geom_line(linewidth = 0.8, color = "black") +
  geom_point(size = 3.5, color = "#4C78A8") +
  geom_point(
    data = panelC_df %>% filter(herv_contribution == "predominant"),
    size = 4.5,
    color = "#C21F1F"
  ) +
  geom_text(aes(label = label), vjust = -0.8, size = 3.5, lineheight = 0.95) +
  scale_x_discrete(
    labels = c(
      "minor" = "<10%",
      "substantial" = "10–50%",
      "dominant" = "50–80%",
      "predominant" = ">80%"
    )
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, max(panelC_df$frac_annotated) * 1.25),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "HERV contribution to lncRNA gene body",
    y = "Annotated lncRNAs"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    plot.margin = margin(8, 8, 8, 8)
  )

pC


pA
pB
pC

library(readr)
library(dplyr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(scales)
library(forcats)
library(patchwork)
library(cowplot)
plot_margin <- theme(plot.margin = margin(10, 10, 10, 10)) # top, right, bottom, left (in pts)
p_multi <- plot_grid(
  pA + plot_margin,
  pB + plot_margin, 
  pC + plot_margin,
  labels = c("A", "B", "C"),
  ncol = 3,
  label_size = 11,
  label_fontface = "bold"
)
p_multi


ggsave(
  filename = file.path(CFG$fig_dir, "figure_herv_lncrna_annotations.png"),
  plot = p_multi,
  width = 20,
  height = 9
)

ggsave(
  filename = file.path(CFG$fig_dir, "figure_herv_lncrna_annotations.svg"),
  plot = p_multi,
  width = 20,
  height = 9
)
