suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
IN_CANDIDATES <- "../../../data/ltr_candidates_antisense_last_exon_or_utr3.protein_coding_or_lncRNA.tsv"
IN_U3R        <- "../../../data/U3R_U5_catalogue.tsv"

OUTDIR <- "../../../results/sparcs_like/"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------
# Helper
# ------------------------------------------------------------------
infer_subfamily <- function(x) stringr::str_replace(x, "_merged.*$", "")

# ------------------------------------------------------------------
# Read data
# ------------------------------------------------------------------
cand <- read_tsv(IN_CANDIDATES, show_col_types = FALSE) %>% as.data.table()
u3r  <- read_tsv(IN_U3R,        show_col_types = FALSE) %>% as.data.table()

# ------------------------------------------------------------------
# Standardize
# ------------------------------------------------------------------
cand[, subfamily := infer_subfamily(a_name)]
u3r[, subfamily := infer_subfamily(name)]

# Make a unified LTR name column
cand[, name := a_name]

# ------------------------------------------------------------------
# Define SPARCS-like candidates
# ------------------------------------------------------------------
# Conservative definition:
#   - protein_coding only if overlap is in 3'UTR
#   - lncRNA only if overlap is in last exon
cand_pc  <- cand[gene_type == "protein_coding" & feature_class == "three_prime_UTR"]
cand_lnc <- cand[gene_type == "lncRNA" & feature_class == "last_exon"]

sparcs_like <- rbindlist(list(cand_pc, cand_lnc), use.names = TRUE, fill = TRUE)

# Optional: if you want the older broader lncRNA definition, uncomment:
# cand_lnc <- cand[gene_type == "lncRNA"]

# Collapse to one row per LTR, because one LTR can map to multiple transcripts
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

# ------------------------------------------------------------------
# Join SPARCS-like set to U3R/U5 catalogue
# ------------------------------------------------------------------
u3r_dt <- copy(u3r)

# derived metrics
u3r_dt[, ltr_len := end - start]
u3r_dt[, has_u5_boundary := !is.na(u5_start_rel)]
u3r_dt[, has_u3_boundary := !is.na(u3_end_rel)]
u3r_dt[, reconstructed_both := has_u3_boundary & has_u5_boundary]
u3r_dt[, R_len := ifelse(reconstructed_both, u5_start_rel - u3_end_rel - 1, NA_real_)]

# simplify status if useful
u3r_dt[, status_simple := fifelse(status == "OK", "OK", "LOW_CONF")]

# annotate group membership
u3r_dt[, group := "background_all"]
u3r_dt[name %in% sparcs_like_uniq$name, group := "SPARCS_like"]

# for solo-only comparisons
u3r_dt[, group_solo := fifelse(role == "solo", "solo_background", NA_character_)]
u3r_dt[role == "solo" & name %in% sparcs_like_uniq$name, group_solo := "SPARCS_like"]

# merged candidate table
sparcs_u3r <- merge(
  sparcs_like_uniq,
  u3r_dt,
  by = "name",
  all.x = TRUE
)

# Save merged table
fwrite(sparcs_u3r, file.path(OUTDIR, "SPARCS_like_candidates_with_U3R_U5_catalogue.tsv"), sep = "\t")

# ------------------------------------------------------------------
# Basic checks
# ------------------------------------------------------------------
cat("Total candidate transcript-level rows:", nrow(sparcs_like), "\n")
cat("Unique SPARCS-like LTRs:", nrow(sparcs_like_uniq), "\n")
cat("SPARCS-like LTRs found in U3R/U5 catalogue:", sum(!is.na(sparcs_u3r$chrom)), "\n")

# ------------------------------------------------------------------
# Summary tables
# ------------------------------------------------------------------
summary_by_group <- u3r_dt[
  ,
  .(
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
  ),
  by = group
]

summary_solo <- u3r_dt[role == "solo",
                       .(
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
                       ),
                       by = group_solo
]

fwrite(summary_by_group, file.path(OUTDIR, "summary_vs_all_LTRs.tsv"), sep = "\t")
fwrite(summary_solo,     file.path(OUTDIR, "summary_vs_soloLTRs.tsv"), sep = "\t")

print(summary_by_group)
print(summary_solo)

# ------------------------------------------------------------------
# Statistical tests
# ------------------------------------------------------------------

# 1) SPARCS-like vs all LTRs: proportion of OK
tab_ok_all <- table(
  u3r_dt$group,
  u3r_dt$status == "OK",
  useNA = "no"
)
print(tab_ok_all)
fisher_ok_all <- fisher.test(tab_ok_all)
print(fisher_ok_all)

# 2) SPARCS-like vs solo LTRs: proportion of OK
tab_ok_solo <- table(
  u3r_dt[role == "solo"]$group_solo,
  u3r_dt[role == "solo"]$status == "OK",
  useNA = "no"
)
print(tab_ok_solo)
fisher_ok_solo <- fisher.test(tab_ok_solo)
print(fisher_ok_solo)

# 3) SPARCS-like vs all LTRs: confidence distribution
wilcox_conf_all <- wilcox.test(
  confidence ~ group,
  data = as.data.frame(u3r_dt[group %in% c("SPARCS_like", "background_all")]),
  exact = FALSE
)
print(wilcox_conf_all)

# 4) SPARCS-like vs solo LTRs: confidence distribution
wilcox_conf_solo <- wilcox.test(
  confidence ~ group_solo,
  data = as.data.frame(u3r_dt[role == "solo" & group_solo %in% c("SPARCS_like", "solo_background")]),
  exact = FALSE
)
print(wilcox_conf_solo)

# 5) reconstructed_both proportions
tab_both_all <- table(
  u3r_dt$group,
  u3r_dt$reconstructed_both,
  useNA = "no"
)
fisher_both_all <- fisher.test(tab_both_all)
print(tab_both_all)
print(fisher_both_all)

tab_both_solo <- table(
  u3r_dt[role == "solo"]$group_solo,
  u3r_dt[role == "solo"]$reconstructed_both,
  useNA = "no"
)
fisher_both_solo <- fisher.test(tab_both_solo)
print(tab_both_solo)
print(fisher_both_solo)

# Save a compact stats report
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

write_tsv(stats_report, file.path(OUTDIR, "stats_report.tsv"))

# ------------------------------------------------------------------
# Per-subfamily enrichment of SPARCS-like candidates among OK LTRs
# ------------------------------------------------------------------
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

fwrite(subfam_tab, file.path(OUTDIR, "subfamily_summary.tsv"), sep = "\t")

# ------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------

# 1) Confidence vs groups (all)
p1 <- u3r_dt[group %in% c("SPARCS_like", "background_all")] %>%
  as_tibble() %>%
  ggplot(aes(x = group, y = confidence)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 1) +
  theme_bw(base_size = 12) +
  labs(
    title = "Reconstruction confidence: SPARCS-like vs all LTRs",
    x = "",
    y = "Confidence score"
  )

ggsave(file.path(OUTDIR, "confidence_SPARCS_vs_all.pdf"), p1, width = 5, height = 4)

# 2) Confidence vs groups (solo only)
p2 <- u3r_dt[role == "solo" & group_solo %in% c("SPARCS_like", "solo_background")] %>%
  as_tibble() %>%
  ggplot(aes(x = group_solo, y = confidence)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 1) +
  theme_bw(base_size = 12) +
  labs(
    title = "Reconstruction confidence: SPARCS-like vs solo LTRs",
    x = "",
    y = "Confidence score"
  )

ggsave(file.path(OUTDIR, "confidence_SPARCS_vs_solo.pdf"), p2, width = 5, height = 4)

# 3) OK fractions
plot_ok_df <- bind_rows(
  u3r_dt[group %in% c("SPARCS_like", "background_all")] %>%
    as_tibble() %>%
    count(group, status) %>%
    group_by(group) %>%
    mutate(frac = n / sum(n), comparison = "All LTRs"),
  u3r_dt[role == "solo" & group_solo %in% c("SPARCS_like", "solo_background")] %>%
    as_tibble() %>%
    count(group = group_solo, status) %>%
    group_by(group) %>%
    mutate(frac = n / sum(n), comparison = "Solo LTRs")
)

p3 <- ggplot(plot_ok_df, aes(x = group, y = frac, fill = status)) +
  geom_col(position = "stack") +
  facet_wrap(~comparison, scales = "free_x") +
  theme_bw(base_size = 12) +
  labs(
    title = "Status distribution",
    x = "",
    y = "Fraction"
  )

ggsave(file.path(OUTDIR, "status_fraction_plots.pdf"), p3, width = 7, height = 4)

# 4) Reconstructed both boundaries
plot_both_df <- bind_rows(
  u3r_dt[group %in% c("SPARCS_like", "background_all")] %>%
    as_tibble() %>%
    count(group, reconstructed_both) %>%
    group_by(group) %>%
    mutate(frac = n / sum(n), comparison = "All LTRs"),
  u3r_dt[role == "solo" & group_solo %in% c("SPARCS_like", "solo_background")] %>%
    as_tibble() %>%
    count(group = group_solo, reconstructed_both) %>%
    group_by(group) %>%
    mutate(frac = n / sum(n), comparison = "Solo LTRs")
)

p4 <- ggplot(plot_both_df, aes(x = group, y = frac, fill = reconstructed_both)) +
  geom_col(position = "stack") +
  facet_wrap(~comparison, scales = "free_x") +
  theme_bw(base_size = 12) +
  labs(
    title = "Fraction with both U3 and U5 boundaries reconstructed",
    x = "",
    y = "Fraction"
  )

ggsave(file.path(OUTDIR, "reconstructed_both_fraction_plots.pdf"), p4, width = 7, height = 4)

# ------------------------------------------------------------------
# Candidate-focused output ranked by confidence
# ------------------------------------------------------------------
sparcs_ranked <- sparcs_u3r %>%
  as_tibble() %>%
  arrange(desc(status == "OK"), desc(confidence), desc(reconstructed_both), desc(n_tx))

write_tsv(sparcs_ranked, file.path(OUTDIR, "SPARCS_like_candidates_ranked.tsv"))

cat("Done.\n")