#!/usr/bin/env Rscript
# ==============================================================================
# build_gene_overlap_tables.R
# 03_transcript_context step 11 (frozen global order 18, effective order 19), LOCAL
# ==============================================================================
# HERV x gene overlap summaries and the Figure 2A source data, computed ONLY from
# the accepted integrated table results/integration/transcript_context.with_herv_id.tsv
# (53,176 x 42, effective order 18). No annotation is recomputed: host genes, gene
# types and exon/intron evidence are parsed from integrated fields.
#
# Historical lineage (reproduced byte for byte; see PROVENANCE.md, step 11):
#   B_transcript_context/01_build_herv_gene_overlap_tables.R  (+ utils_transcript_context.R)
#     ran first (2026-04-28 12:34) and wrote 8 tables to results/lncRNA/; only its two
#     *_feature_overlap_summary.tsv files survive as its output.
#   B_transcript_context/old/03_hervs_and_lncRNA.R
#     ran next (12:35-12:37), OVERWROTE the six locus/gene tables with its own schemas
#     and wrote the Figure 2A source (results/figure_herv_transcript_overlap/).
# Each canonical output follows the producer of the corresponding historical file.
#
# Rules (historical, unchanged)
#   host gene      : closest_gene_same_strand_distance_bp == 0 and gene_type not NA/"."/"";
#                    every tied host gene in closest_gene_same_strand_name is kept
#   lncRNA host    : name token matching  [^,|]+\|[^,|]+\|[^,|]*lncRNA[^,|]*
#   protein-coding : name token matching  [^,|]+\|[^,|]+\|protein_coding
#   locus overlap  : str_detect(host gene_type, "lncRNA" | "protein_coding")
#   exon / intron  : some host gene_id g has an entry matching ^g.*:.*exon (or intron)
#                    in the per-gene feature field (';'-separated). Not exclusive: a locus
#                    can be exonic and intronic, lncRNA and protein-coding.
#   evidence column: old/03 (tables, Figure 2A) reads feature_overlap_per_gene (any
#                    strand); 01_build (feature summaries) reads
#                    feature_overlap_per_gene_same_strand. Both are computed; under the
#                    D022 per-gene assets they must agree on every locus (checked).
#   panel A        : exclusive category, priority PC exon > lncRNA exon > PC intron >
#                    lncRNA intron > other; the plotted counts are NOT exclusive:
#                    sum(is_X & overlap_Y) over all loci, percent = n / all loci.
#
# Canonical change: panelA_global_overlap_context_df.tsv keeps only the columns the
# Figure 2A classification uses. The historical file dumped all 69 atlas/derived
# columns, re-parsing RBP, LTR and domain fields and passing ltr3_status through; the
# canonical step reads none of them (D017 / D018 impact NONE).
#
# Usage:
#   Rscript build_gene_overlap_tables.R --input <transcript_context.with_herv_id.tsv> --outdir <dir>
#     [--expected-sha256 <hex>]            default: the accepted order-18 table
#     [--allow-evidence-disagreement]      validation re-runs on the April table only
# ==============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(purrr)
})

ACCEPTED_SHA256 <- "cbaead1273c916fbcfe4337c7602a9d12c1d5b8bdacfd2a0dd05c0b51f662995"
N_LOCI <- 53176L

# ------------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------------
parse_args <- function(argv) {
  a <- list(input = NULL, outdir = NULL, expected_sha256 = ACCEPTED_SHA256,
            allow_evidence_disagreement = FALSE)
  i <- 1
  while (i <= length(argv)) {
    k <- argv[[i]]
    if (k == "--input") { a$input <- argv[[i + 1]]; i <- i + 2 }
    else if (k == "--outdir") { a$outdir <- argv[[i + 1]]; i <- i + 2 }
    else if (k == "--expected-sha256") { a$expected_sha256 <- argv[[i + 1]]; i <- i + 2 }
    else if (k == "--allow-evidence-disagreement") { a$allow_evidence_disagreement <- TRUE; i <- i + 1 }
    else stop("unknown argument: ", k)
  }
  if (is.null(a$input) || is.null(a$outdir)) stop("usage: --input <table> --outdir <dir>")
  a
}

sha256_file <- function(path) {
  out <- system2("sha256sum", shQuote(path), stdout = TRUE)
  substr(out[[1]], 1, 64)
}

# ------------------------------------------------------------------------------
# Parsing helpers (old/03_hervs_and_lncRNA.R, verbatim semantics)
# ------------------------------------------------------------------------------
as_num <- function(x) suppressWarnings(as.numeric(x))
infer_subfamily <- function(x) str_remove(x, "-int_.*$")
infer_class <- function(subfamily) case_when(str_detect(subfamily, "^HERV") ~ "ERV", TRUE ~ "Other")

# Host gene type is trusted only for a true same-strand overlap (distance 0).
flag_overlap_gene_type_same_strand <- function(distance_bp, gene_type) {
  case_when(
    is.na(distance_bp) ~ NA_character_,
    distance_bp != 0 ~ NA_character_,
    is.na(gene_type) ~ NA_character_,
    gene_type %in% c(".", "") ~ NA_character_,
    TRUE ~ gene_type
  )
}

HOST_PATTERN <- c(
  lncRNA         = "[^,|]+\\|[^,|]+\\|[^,|]*lncRNA[^,|]*",
  protein_coding = "[^,|]+\\|[^,|]+\\|protein_coding"
)

# All host genes of one target class named in closest_gene_same_strand_name
# (ties are kept: the field lists every same-strand gene at the minimum distance).
extract_overlapping_genes <- function(x, target) {
  empty <- tibble(gene_stable_id = character(0), gene_name = character(0), gene_type = character(0))
  if (is.na(x) || x == "") return(empty)
  hits <- regmatches(x, gregexpr(HOST_PATTERN[[target]], x, perl = TRUE))[[1]]
  if (length(hits) == 0 || (length(hits) == 1 && hits[1] == "-1")) return(empty)
  parts <- str_split_fixed(hits, "\\|", 3)
  distinct(tibble(gene_stable_id = parts[, 1], gene_name = parts[, 2], gene_type = parts[, 3]))
}

collapse_gene_field <- function(df, col) {
  vals <- unique(df[[col]])
  vals <- vals[!is.na(vals) & vals != ""]
  paste(vals, collapse = ";")
}

# TRUE if some entry of the ';'-separated per-gene field matches ^<gene_id>.*:.*<feature>
is_feature_for_gene <- function(feature_per_gene, gene_id, feature) {
  if (is.na(feature_per_gene) || is.na(gene_id)) return(FALSE)
  entries <- str_split(feature_per_gene, ";")[[1]]
  any(str_detect(entries, paste0("^", gene_id, ".*:.*", feature)))
}

any_host_feature <- function(host_lists, per_gene, feature) {
  map2_lgl(host_lists, per_gene, function(h, fpg) {
    any(vapply(h$gene_stable_id, function(g) is_feature_for_gene(fpg, g, feature), logical(1)))
  })
}

# ------------------------------------------------------------------------------
# Load the accepted integrated table and derive the classification inputs
# ------------------------------------------------------------------------------
load_integrated_table <- function(path) {
  raw <- read_tsv(path, col_types = cols(.default = col_character()), progress = FALSE)
  needed <- c("HERV_id", "locid", "chrom", "start", "end", "strand", "length",
              "closest_gene_same_strand_name", "closest_gene_same_strand_gene_type",
              "closest_gene_same_strand_distance_bp", "feature_overlap",
              "feature_overlap_per_gene", "feature_overlap_per_gene_same_strand")
  miss <- setdiff(needed, names(raw))
  if (length(miss)) stop("input lacks columns: ", paste(miss, collapse = ", "))
  raw %>%
    select(all_of(needed)) %>%
    mutate(
      subfamily    = infer_subfamily(HERV_id),
      class        = infer_class(subfamily),
      coord_chr    = chrom,
      coord_start  = start,            # character, as historically
      coord_end    = end,
      coord_strand = strand,
      length_num   = as_num(length),
      overlap_gene_same_strand_gene_type = flag_overlap_gene_type_same_strand(
        as_num(closest_gene_same_strand_distance_bp), closest_gene_same_strand_gene_type)
    )
}

# Per-locus host genes, overlap flags and exon/intron flags under one evidence column.
classify_loci <- function(tab, evidence_col) {
  lnc <- map(tab$closest_gene_same_strand_name, extract_overlapping_genes, target = "lncRNA")
  pc  <- map(tab$closest_gene_same_strand_name, extract_overlapping_genes, target = "protein_coding")
  gt  <- tab$overlap_gene_same_strand_gene_type
  ev  <- tab[[evidence_col]]
  tibble(
    overlap_lncrna         = !is.na(gt) & str_detect(gt, "lncRNA"),
    overlap_protein_coding = !is.na(gt) & str_detect(gt, "protein_coding"),
    is_exonic_lncrna   = any_host_feature(lnc, ev, "exon"),
    is_intronic_lncrna = any_host_feature(lnc, ev, "intron"),
    is_exonic_pc       = any_host_feature(pc, ev, "exon"),
    is_intronic_pc     = any_host_feature(pc, ev, "intron"),
    hosts_lncrna = lnc,
    hosts_pc     = pc
  )
}

# ------------------------------------------------------------------------------
# Table builders (old/03 schemas: the historical files downstream modules read)
# ------------------------------------------------------------------------------
HERV_LEVEL_LEAD <- c("HERV_id", "subfamily", "class", "length_num", "coord_strand", "coord_chr",
                     "coord_start", "coord_end", "closest_gene_same_strand_distance_bp",
                     "closest_gene_same_strand_name", "closest_gene_same_strand_gene_type",
                     "feature_overlap")

build_herv_level_table <- function(sub, hosts, prefix, n_col, with_per_gene = FALSE) {
  sub %>%
    mutate(
      !!paste0("overlap_", prefix, "_gene_ids")   := map_chr(hosts, collapse_gene_field, "gene_stable_id"),
      !!paste0("overlap_", prefix, "_gene_names") := map_chr(hosts, collapse_gene_field, "gene_name"),
      !!paste0("overlap_", prefix, "_gene_types") := map_chr(hosts, collapse_gene_field, "gene_type"),
      !!n_col := map_int(hosts, nrow)
    ) %>%
    select(all_of(HERV_LEVEL_LEAD), any_of(if (with_per_gene) "feature_overlap_per_gene"),
           starts_with(paste0("overlap_", prefix, "_gene_")), all_of(n_col))
}

# One row per unique host gene; loci counted once per gene (n_distinct), ids sorted,
# coordinates in table order. `extra` = columns kept through the unnest step.
build_gene_level_table <- function(sub, hosts, extra = character(0), exonic = FALSE, with_features = FALSE) {
  n_col  <- if (exonic) "n_exonic_overlapping_hervs" else "n_overlapping_hervs"
  id_col <- if (exonic) "exonic_overlapping_herv_ids" else "overlapping_herv_ids"
  long <- sub %>%
    mutate(row_id_tmp = row_number(), extracted_genes = hosts) %>%
    select(row_id_tmp, HERV_id, subfamily, class, coord_chr, coord_start, coord_end, coord_strand,
           all_of(extra), extracted_genes) %>%
    unnest(extracted_genes) %>%
    distinct()
  out <- long %>%
    group_by(gene_stable_id, gene_name, gene_type) %>%
    summarise(
      !!n_col  := n_distinct(HERV_id),
      !!id_col := paste(sort(unique(HERV_id)), collapse = ";"),
      overlapping_subfamilies = paste(sort(unique(subfamily)), collapse = ";"),
      overlapping_features = if (with_features) paste(sort(unique(feature_overlap)), collapse = ";") else NA_character_,
      overlapping_coords = paste(paste0(coord_chr, ":", coord_start, "-", coord_end, "(", coord_strand, ")"),
                                 collapse = ";"),
      .groups = "drop"
    )
  if (!with_features) out <- select(out, -overlapping_features)
  arrange(out, desc(.data[[n_col]]), gene_name)
}

# 01_build_herv_gene_overlap_tables.R: exon/intron counts within a class subset
# (same-strand evidence), percent of that subset.
feature_summary <- function(exonic, intronic) {
  tibble(feature = c("exon", "intron"), n = c(sum(exonic, na.rm = TRUE), sum(intronic, na.rm = TRUE))) %>%
    mutate(percent = 100 * n / length(exonic))
}

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------
main <- function() {
  a <- parse_args(commandArgs(trailingOnly = TRUE))
  got <- sha256_file(a$input)
  if (got != a$expected_sha256) stop("input sha256 ", got, " != expected ", a$expected_sha256)
  message("[ok] input sha256 ", got)

  tab <- load_integrated_table(a$input)
  if (nrow(tab) != N_LOCI || n_distinct(tab$locid) != N_LOCI || n_distinct(tab$HERV_id) != N_LOCI)
    stop("expected ", N_LOCI, " loci, one row each")

  any_ev  <- classify_loci(tab, "feature_overlap_per_gene")              # old/03
  same_ev <- classify_loci(tab, "feature_overlap_per_gene_same_strand")  # 01_build
  flag_cols <- c("is_exonic_lncrna", "is_intronic_lncrna", "is_exonic_pc", "is_intronic_pc")
  disagree <- vapply(flag_cols, function(f) sum(any_ev[[f]] != same_ev[[f]]), integer(1))
  message("[qc] any- vs same-strand evidence, loci disagreeing: ",
          paste(names(disagree), disagree, sep = "=", collapse = ", "))
  if (any(disagree > 0) && !a$allow_evidence_disagreement)
    stop("any-strand and same-strand exon/intron evidence disagree; D022 assets expected")

  lnc_rows <- which(any_ev$overlap_lncrna)
  pc_rows  <- which(any_ev$overlap_protein_coding)
  lnc_ex   <- lnc_rows[any_ev$is_exonic_lncrna[lnc_rows]]
  lnc <- tab[lnc_rows, ]; pc <- tab[pc_rows, ]; lncx <- tab[lnc_ex, ]

  out <- list()
  out[["tables/lncRNA_overlapping_HERVs.herv_level.tsv"]] <-
    build_herv_level_table(lnc, any_ev$hosts_lncrna[lnc_rows], "lnc", "n_overlapping_lncrnas")
  out[["tables/protein_coding_overlapping_HERVs.herv_level.tsv"]] <-
    build_herv_level_table(pc, any_ev$hosts_pc[pc_rows], "pc", "n_overlapping_protein_coding")
  out[["tables/lncRNA_overlapping_HERVs.gene_level.tsv"]] <-
    build_gene_level_table(lnc, any_ev$hosts_lncrna[lnc_rows])
  out[["tables/protein_coding_overlapping_HERVs.gene_level.tsv"]] <-
    build_gene_level_table(pc, any_ev$hosts_pc[pc_rows], extra = "feature_overlap", with_features = TRUE)
  out[["tables/lncRNA_overlapping_HERVs.exonic.herv_level.tsv"]] <-
    build_herv_level_table(lncx, any_ev$hosts_lncrna[lnc_ex], "lnc", "n_overlapping_lncrnas", with_per_gene = TRUE)
  out[["tables/lncRNA_overlapping_HERVs.exonic.gene_level.tsv"]] <-
    build_gene_level_table(lncx, any_ev$hosts_lncrna[lnc_ex], extra = c("feature_overlap", "feature_overlap_per_gene"),
                           exonic = TRUE, with_features = TRUE)
  out[["tables/lncRNA_feature_overlap_summary.tsv"]] <-
    feature_summary(same_ev$is_exonic_lncrna[same_ev$overlap_lncrna], same_ev$is_intronic_lncrna[same_ev$overlap_lncrna])
  out[["tables/protein_coding_feature_overlap_summary.tsv"]] <-
    feature_summary(same_ev$is_exonic_pc[same_ev$overlap_protein_coding], same_ev$is_intronic_pc[same_ev$overlap_protein_coding])

  # Figure 2A source data (old/03 panel A)
  panel_df <- tab %>%
    select(HERV_id, locid, chrom, start, end, strand, subfamily, class,
           closest_gene_same_strand_distance_bp, closest_gene_same_strand_name,
           closest_gene_same_strand_gene_type, overlap_gene_same_strand_gene_type,
           feature_overlap_per_gene, feature_overlap_per_gene_same_strand) %>%
    bind_cols(select(any_ev, all_of(c("overlap_lncrna", "overlap_protein_coding", flag_cols)))) %>%
    mutate(panelA_category = case_when(
      overlap_protein_coding & is_exonic_pc       ~ "Protein-coding exon",
      overlap_lncrna         & is_exonic_lncrna   ~ "lncRNA exon",
      overlap_protein_coding & is_intronic_pc     ~ "Protein-coding intron",
      overlap_lncrna         & is_intronic_lncrna ~ "lncRNA intron",
      TRUE                                        ~ "Other / no same-strand overlap"))
  panel_counts <- tibble(
    category = c("lncRNA intron", "lncRNA exon", "PC intron", "PC exon"),
    n = c(sum(panel_df$is_intronic_lncrna & panel_df$overlap_lncrna, na.rm = TRUE),
          sum(panel_df$is_exonic_lncrna   & panel_df$overlap_lncrna, na.rm = TRUE),
          sum(panel_df$is_intronic_pc     & panel_df$overlap_protein_coding, na.rm = TRUE),
          sum(panel_df$is_exonic_pc       & panel_df$overlap_protein_coding, na.rm = TRUE))
  ) %>% mutate(percent = n / nrow(tab))
  out[["figure2a_source/panelA_global_overlap_context_df.tsv"]] <- panel_df
  out[["figure2a_source/panelA_global_overlap_context.tsv"]] <- panel_counts

  # Headline counts, derived from the outputs above (for manuscript synchronisation later)
  lnc_gene <- out[["tables/lncRNA_overlapping_HERVs.gene_level.tsv"]]
  pc_gene  <- out[["tables/protein_coding_overlapping_HERVs.gene_level.tsv"]]
  out[["summary_counts.tsv"]] <- tibble(
    metric = c("lncRNA_same_strand_loci", "lncRNA_unique_genes", "protein_coding_same_strand_loci",
               "protein_coding_unique_genes", "lncRNA_intron_loci", "lncRNA_exon_loci",
               "protein_coding_intron_loci", "protein_coding_exon_loci", "lncRNA_exonic_genes"),
    n = c(nrow(out[["tables/lncRNA_overlapping_HERVs.herv_level.tsv"]]), nrow(lnc_gene),
          nrow(out[["tables/protein_coding_overlapping_HERVs.herv_level.tsv"]]), nrow(pc_gene),
          panel_counts$n, nrow(out[["tables/lncRNA_overlapping_HERVs.exonic.gene_level.tsv"]]))
  ) %>% mutate(percent_of_all_loci = 100 * n / nrow(tab))

  # Staged, atomic publication
  stage <- file.path(dirname(normalizePath(a$outdir, mustWork = FALSE)),
                     paste0(".stage_gene_overlap.", Sys.getpid()))
  unlink(stage, recursive = TRUE)
  for (p in names(out)) {
    dir.create(dirname(file.path(stage, p)), recursive = TRUE, showWarnings = FALSE)
    write_tsv(out[[p]], file.path(stage, p))
    message("[ok] ", p, ": ", nrow(out[[p]]), " rows x ", ncol(out[[p]]), " cols")
  }
  writeLines(capture.output(sessionInfo()), file.path(stage, "sessionInfo.txt"))
  unlink(a$outdir, recursive = TRUE)
  if (!file.rename(stage, a$outdir)) stop("could not publish ", a$outdir)
  message("[OK] wrote ", a$outdir)
}

main()
