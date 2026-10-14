#!/usr/bin/env Rscript
# 08_figures_supplement, effective order 54.1 (D035): Supplemental Figure 2, GO BP of the TFs depleted in lncRNA SPARCS-like LTRs.
#
# PLOTTING ONLY. The accepted order-45 table (07_sparcs/results/enrichment/go/sparcs_lnc_depleted_GO_BP_over_depleted_TFsymbols.tsv,
# 7 terms, D033) is read as is and placed in a DOSE enrichResult container (no enrichGO, no test, no adjustment: the container
# only holds the accepted rows). The figure is drawn by the historical call of E_sparcs/01_analysis_sparcs.R run_tf_go(), verbatim:
#   enrichplot::dotplot(ego, showCategory = 20, font.size = 12) + ggtitle("GO: BP — Depleted TF motifs (SPARCS-like lncRNA LTRs)")
#   ggsave(png, width = 10, height = 7, dpi = 300, bg = "white")
# so x = GeneRatio, point size = Count, colour = p.adjust, terms ordered by enrichplot (GeneRatio), as historically.
#
# --mode canonical               the accepted order-45 table -> figures + panel sources
# --mode historical_reproduction the historical 5-term table (07's checksum-recorded historical reproduction copy) -> a PNG
#                                 for the historical comparison only (never a canonical output)
# Usage: plot_supplementary_figure2.R --mode <mode> --release <ncRNA root> --figures-out <dir> --sources-out <dir>
suppressPackageStartupMessages({
  library(DOSE); library(enrichplot); library(ggplot2); library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) { if (is.null(default)) stop("missing ", flag); default } else args[[i + 1]] }
MODE <- arg("--mode"); REL <- normalizePath(arg("--release"), mustWork = TRUE)
FIG <- arg("--figures-out"); SRCOUT <- arg("--sources-out", "")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------------------------------------- parameters
TITLE <- "GO: BP — Depleted TF motifs (SPARCS-like lncRNA LTRs)"   # historical title, verbatim
SHOW_CATEGORY <- 20; FONT_SIZE <- 12                                 # historical dotplot arguments
WIDTH_IN <- 10; HEIGHT_IN <- 7; DPI <- 300                           # historical ggsave arguments
CONTAINER_P_CUTOFF <- 0.10; CONTAINER_Q_CUTOFF <- 0.20               # order-45 tf_go cutoffs (config/parameters.yaml), recorded in the container
FIXED_PDF_DATE <- "20261001000000"
GO_TABLE <- "go/sparcs_lnc_depleted_GO_BP_over_depleted_TFsymbols.tsv"
SOURCE <- switch(MODE,
  canonical = list(path = file.path(REL, "07_sparcs/results/enrichment", GO_TABLE),
                   records = c(file.path(REL, "07_sparcs/results/enrichment/checksums.sha256"), file.path(REL, "07_sparcs/checksums.sha256")),
                   keys = c(GO_TABLE, paste0("results/enrichment/", GO_TABLE))),
  historical_reproduction = list(path = file.path(REL, "07_sparcs/validation/enrichment/historical_reproduction", GO_TABLE),
                   records = file.path(REL, "07_sparcs/validation/enrichment/historical_reproduction/checksums.sha256"), keys = GO_TABLE),
  stop("unknown --mode ", MODE))

# ------------------------------------------------------------------------------------------------------------- source guard
recorded <- unlist(lapply(seq_along(SOURCE$records), function(i) {
  l <- readLines(SOURCE$records[i]); sha <- sub("  .*$", "", l); p <- sub("^\\./", "", sub("^[0-9a-f]{64}  ", "", l))
  sha[p == SOURCE$keys[i]]
}))
actual <- digest::digest(file = SOURCE$path, algo = "sha256")
if (length(unique(recorded)) != 1 || recorded[1] != actual) stop("source table is not the accepted (checksum-recorded) table: ", SOURCE$path)

tab <- fread(SOURCE$path, sep = "\t", colClasses = "character", na.strings = NULL)
res <- as.data.frame(tab, stringsAsFactors = FALSE)
for (cl in c("pvalue", "p.adjust", "qvalue")) res[[cl]] <- as.numeric(res[[cl]])
res$Count <- as.integer(res$Count)
rownames(res) <- res$ID

# ------------------------------------------------------------------------------------------------------------- container + plot
ego <- new("enrichResult", result = res, pvalueCutoff = CONTAINER_P_CUTOFF, pAdjustMethod = "BH", qvalueCutoff = CONTAINER_Q_CUTOFF,
           organism = "Homo sapiens", ontology = "BP", gene = character(), keytype = "SYMBOL", readable = TRUE)
if (nrow(as.data.frame(ego)) != nrow(res)) stop("the container would drop accepted rows")
p <- enrichplot::dotplot(ego, showCategory = SHOW_CATEGORY, font.size = FONT_SIZE) + ggtitle(TITLE)

stem <- if (MODE == "canonical") "Supplemental_Figure_2" else "Supplemental_Figure_2_historical_reproduction"
ggsave(file.path(FIG, paste0(stem, ".png")), p, width = WIDTH_IN, height = HEIGHT_IN, dpi = DPI, bg = "white")
if (MODE == "canonical") {
  ggsave(file.path(FIG, paste0(stem, ".svg")), p, width = WIDTH_IN, height = HEIGHT_IN, bg = "white")
  pdf_path <- file.path(FIG, paste0(stem, ".pdf"))
  ggsave(pdf_path, p, width = WIDTH_IN, height = HEIGHT_IN, bg = "white")
  # grDevices::pdf stamps the wall-clock time: overwrite the 14 digits in place (binary-safe, same length, offsets unchanged)
  b <- readBin(pdf_path, "raw", file.size(pdf_path))
  for (key in c("/CreationDate (D:", "/ModDate (D:")) {
    at <- grepRaw(key, b, fixed = TRUE, all = TRUE)
    if (length(at) != 1) stop("expected one ", key, " in the PDF")
    i <- at + nchar(key) + 0:13
    if (!all(b[i] >= charToRaw("0") & b[i] <= charToRaw("9"))) stop("unexpected date format after ", key)
    b[i] <- charToRaw(FIXED_PDF_DATE)
  }
  writeBin(b, pdf_path)
}

# ------------------------------------------------------------------------------------------------------------- panel sources
if (MODE == "canonical" && nzchar(SRCOUT)) {
  dir.create(SRCOUT, recursive = TRUE, showWarnings = FALSE)
  fwrite(tab, file.path(SRCOUT, "supplementary_figure2_source_data.tsv"), sep = "\t", quote = FALSE)  # the accepted fields, unparsed
  bd <- ggplot_build(p)$data[[1]]; lv <- levels(p$data$Description)
  map <- data.table(term_position_from_bottom = bd$y, Description = lv[bd$y], x_GeneRatio_plotted = sprintf("%.10f", bd$x),
                    point_size_aesthetic = sprintf("%.6f", bd$size), point_fill = bd$fill)
  map <- merge(map, tab[, .(Description, ID, GeneRatio, Count, p.adjust)], by = "Description", sort = FALSE)
  fwrite(map[order(-term_position_from_bottom)], file.path(SRCOUT, "supplementary_figure2_plot_mapping.tsv"), sep = "\t", quote = FALSE)
}
cat("[OK] Supplemental Figure 2 (", MODE, "): ", nrow(res), " GO terms from ", sub(paste0(REL, "/"), "", SOURCE$path), "; sans -> ",
    systemfonts::match_fonts("sans")$path, "\n", sep = "")
