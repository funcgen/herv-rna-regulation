#!/usr/bin/env Rscript
# ==============================================================================
# 05_lncrna step 3 (frozen global order 36, effective order 37), LOCAL:
# assign ONE canonical HERV-coverage category to every lncRNA gene (D029).
#
# Input   the accepted D027-corrected order-35 table
#         results/herv_coverage/lncRNA_herv_coverage.tsv (checked against its manifest)
# Value   x = herv_bp_in_gene_body / gene_body_bp, the EXACT integer ratio. The rounded
#         gene_body_herv_fraction / _percent columns are carried for display only and
#         never classify.
# Rule    config/parameters.yaml lncrna.coverage_bins, boundary_rule "right_closed"
#         (the manuscript Methods definition):
#           minor        x <= minor_max                       (<= 10%)
#           substantial  minor_max < x <= substantial_max     (> 10-50%)
#           dominant     substantial_max < x <= dominant_max  (> 50-80%)
#           predominant  x > dominant_max                     (> 80%)
#         A threshold t = p/q (read as an exact decimal) is compared as
#         herv_bp * q <= p * gene_body_bp, in integers (all products < 2^53).
# Output  results/coverage_categories/lncRNA_herv_coverage_categories.tsv (one row per gene)
#         results/coverage_categories/lncRNA_herv_coverage_category_counts.tsv
# This is the ONLY place categories are assigned; consumers read herv_coverage_category.
# Historical equivalents: the binning blocks of 04_annotate_lncrna_overlap_tables.R and
# exploration_lit_mining.R (the frozen plan's 06_plot_herv_lncrna_overlap.R bins nothing).
# ==============================================================================
suppressPackageStartupMessages({
  library(yaml)
  library(digest)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("missing argument ", flag)
  args[[i + 1]]
}
coverage_path <- get_arg("--coverage")
manifest_path <- get_arg("--coverage-manifest")
params_path   <- get_arg("--params")
outdir        <- get_arg("--outdir")

CATEGORIES <- c("minor", "substantial", "dominant", "predominant")
LABELS     <- c(minor = "≤10%", substantial = ">10–50%", dominant = ">50–80%", predominant = ">80%")

# ---------------------------------------------------------------- guards
verify_against_manifest <- function(path, manifest) {
  rel <- basename(path)
  m <- read.table(manifest, header = FALSE, col.names = c("sha", "file"), colClasses = "character")
  want <- m$sha[m$file %in% c(rel, paste0("./", rel))]
  got <- digest(file = path, algo = "sha256")
  if (length(want) != 1 || got != want) stop(path, " does not match its accepted manifest ", manifest)
  got
}

# exact decimal "0.10" -> c(p = 10, q = 100)
as_decimal_ratio <- function(x) {
  s <- as.character(x)
  if (!grepl("^[0-9]+(\\.[0-9]+)?$", s)) stop("threshold is not a plain decimal: ", s)
  parts <- strsplit(s, ".", fixed = TRUE)[[1]]
  frac <- if (length(parts) == 2) parts[2] else ""
  c(p = as.numeric(paste0(parts[1], frac)), q = 10^nchar(frac))
}

load_rule <- function(path) {
  # thresholds are read as text so that 0.10 stays exactly 10/100
  txt <- readLines(path, warn = FALSE)
  block <- txt[seq(grep("^lncrna:", txt), length(txt))]
  block <- block[seq_len(if (length(e <- grep("^[a-z]", block[-1])) > 0) e[1] else length(block))]
  pick <- function(key) {
    l <- grep(paste0("^\\s+", key, ":"), block, value = TRUE)
    if (length(l) != 1) stop("parameter lncrna.", key, " not found exactly once")
    gsub("[\" ]", "", sub("#.*$", "", sub(paste0("^\\s+", key, ":"), "", l)))
  }
  p <- yaml::read_yaml(path)$lncrna
  rule <- list(boundary_rule = p$boundary_rule, classification_value = p$classification_value,
               minor_max = as_decimal_ratio(pick("minor_max")),
               substantial_max = as_decimal_ratio(pick("substantial_max")),
               dominant_max = as_decimal_ratio(pick("dominant_max")))
  if (!identical(rule$boundary_rule, "right_closed")) stop("lncrna.boundary_rule must be 'right_closed' (D029)")
  if (!identical(rule$classification_value, "exact_ratio")) stop("lncrna.classification_value must be 'exact_ratio' (D029)")
  v <- sapply(rule[c("minor_max", "substantial_max", "dominant_max")], function(r) r[["p"]] / r[["q"]])
  if (!all(diff(v) > 0) || v[1] <= 0 || v[3] >= 1) stop("coverage_bins thresholds must satisfy 0 < minor < substantial < dominant < 1")
  rule
}

# ---------------------------------------------------------------- classification
# x <= p/q  <=>  num * q <= p * den   (exact in double while every product is < 2^53)
at_most <- function(num, den, t) num * t[["q"]] <= t[["p"]] * den

assign_category <- function(num, den, rule) {
  if (any(den <= 0) || any(num < 0) || any(num > den)) stop("coverage numerator/denominator out of range")
  if (max(num * 10^6, den * 10^6) >= 2^53) stop("integer products would exceed 2^53")
  out <- ifelse(at_most(num, den, rule$minor_max), "minor",
         ifelse(at_most(num, den, rule$substantial_max), "substantial",
         ifelse(at_most(num, den, rule$dominant_max), "dominant", "predominant")))
  factor(out, levels = CATEGORIES)
}

interval_text <- function(rule) {
  f <- function(t) format(t[["p"]] / t[["q"]], nsmall = 2)
  c(minor = paste0("x <= ", f(rule$minor_max)),
    substantial = paste0(f(rule$minor_max), " < x <= ", f(rule$substantial_max)),
    dominant = paste0(f(rule$substantial_max), " < x <= ", f(rule$dominant_max)),
    predominant = paste0("x > ", f(rule$dominant_max)))
}

# ---------------------------------------------------------------- main
sha <- verify_against_manifest(coverage_path, manifest_path)
message("[ok] coverage table sha256 ", sha, " (manifest)")
rule <- load_rule(params_path)
cov <- read.delim(coverage_path, colClasses = "character", check.names = FALSE, quote = "", na.strings = character(0))
if (anyDuplicated(cov$gene_stable_id)) stop("duplicated gene_stable_id in the coverage table")
num <- as.numeric(cov$herv_bp_in_gene_body)
den <- as.numeric(cov$gene_body_bp)
if (anyNA(num) || anyNA(den) || any(num != round(num)) || any(den != round(den))) stop("non-integer coverage counts")
category <- assign_category(num, den, rule)

out <- data.frame(
  gene_stable_id = cov$gene_stable_id, gene_name = cov$gene_name, gene_type = cov$gene_type,
  chrom = cov$chrom, strand = cov$strand, gene_start = cov$gene_start, gene_end = cov$gene_end,
  herv_bp_in_gene_body = cov$herv_bp_in_gene_body, gene_body_bp = cov$gene_body_bp,
  exact_ratio = paste0(cov$herv_bp_in_gene_body, "/", cov$gene_body_bp),
  gene_body_herv_fraction_display = cov$gene_body_herv_fraction,
  gene_body_herv_percent_display = cov$gene_body_herv_percent,
  herv_coverage_category = as.character(category),
  herv_coverage_category_label = unname(LABELS[as.character(category)]),
  herv_coverage_category_rank = as.integer(category),
  stringsAsFactors = FALSE
)
n <- as.vector(table(category))
counts <- data.frame(
  herv_coverage_category = CATEGORIES, herv_coverage_category_label = unname(LABELS),
  interval_on_exact_ratio = unname(interval_text(rule)[CATEGORIES]), n = n,
  percent = sprintf("%.1f", 100 * n / sum(n)), stringsAsFactors = FALSE
)
counts <- rbind(counts, data.frame(herv_coverage_category = "total", herv_coverage_category_label = "",
                                   interval_on_exact_ratio = "", n = sum(n), percent = "100.0"))

parent <- dirname(normalizePath(outdir, mustWork = FALSE))
dir.create(parent, recursive = TRUE, showWarnings = FALSE)
stage <- tempfile(pattern = ".stage_coverage_categories.", tmpdir = parent)
dir.create(stage)
w <- function(df, name) {
  con <- file(file.path(stage, name), open = "wb")
  write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
  close(con)
}
w(out, "lncRNA_herv_coverage_categories.tsv")
w(counts, "lncRNA_herv_coverage_category_counts.tsv")
Sys.chmod(list.files(stage, full.names = TRUE), "0644")
if (dir.exists(outdir)) unlink(outdir, recursive = TRUE)
if (!file.rename(stage, outdir)) stop("could not publish ", outdir)
Sys.chmod(outdir, "0755")
message("[ok] ", nrow(out), " genes: ", paste(sprintf("%s %d", CATEGORIES, n), collapse = ", "))
message("[OK] wrote ", outdir)
