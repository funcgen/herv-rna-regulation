#!/usr/bin/env Rscript
# ==============================================================================
# Freeze SRC18 (D034): the MSigDB C7 collection exactly as msigdbr 26.1.0 serves it offline in the release
# environment (MSigDB 2026.1.Hs, from the local msigdbr cache). NOT the historical release: the historical
# order-45 run (2026-04-27) used an msigdbr/MSigDB version that is UNKNOWN (this msigdbr was installed on
# 2026-08-09). Written exactly as E_sparcs/01_analysis_sparcs.R builds its enricher() inputs:
#   TERM2GENE = distinct(gs_name, gene_symbol)     TERM2NAME = distinct(gs_name, gs_description)
# in msigdbr's native row order. Run offline (unshare -rn); refuses any other msigdbr / MSigDB version.
# Usage: freeze_msigdb_c7.R <outdir>
# ==============================================================================
suppressPackageStartupMessages({ library(msigdbr); library(dplyr); library(data.table) })
out <- commandArgs(trailingOnly = TRUE)[1]
if (as.character(packageVersion("msigdbr")) != "26.1.0") stop("msigdbr 26.1.0 required, found ", packageVersion("msigdbr"))
msig_c7 <- suppressWarnings(msigdbr(species = "Homo sapiens", category = "C7"))   # the historical call
dbv <- unique(msig_c7$db_version)
if (!identical(dbv, "2026.1.Hs")) stop("MSigDB 2026.1.Hs required, found ", paste(dbv, collapse = ","))
term2gene <- msig_c7 %>% dplyr::select(gs_name, gene_symbol) %>% distinct()
term2name <- msig_c7 %>% dplyr::select(gs_name, gs_description) %>% distinct()
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fwrite(term2gene, file.path(out, "msigdb_2026.1.Hs_C7_term2gene.tsv"), sep = "\t")
fwrite(term2name, file.path(out, "msigdb_2026.1.Hs_C7_term2name.tsv"), sep = "\t")
cache <- file.path(tools::R_user_dir("msigdbr", "cache"), "msigdb.2026.1.zip")
meta <- data.table(field = c("source", "msigdbr_version", "db_version", "species", "collection", "subcollections", "msigdbr_rows",
                             "gene_sets", "term2gene_pairs", "term2gene_symbols", "msigdbr_cache_zip_sha256", "historical_release", "erratum"),
                   value = c("msigdbr offline (local cache) in the release R environment", as.character(packageVersion("msigdbr")), dbv,
                             "Homo sapiens", "C7", paste(sort(unique(msig_c7$gs_subcollection)), collapse = ","), nrow(msig_c7),
                             uniqueN(term2gene$gs_name), nrow(term2gene), uniqueN(term2gene$gene_symbol),
                             if (file.exists(cache)) digest::digest(file = cache, algo = "sha256") else "NA",
                             "UNKNOWN (not this snapshot)", "D034"))
fwrite(meta, file.path(out, "msigdb_2026.1.Hs_C7.metadata.tsv"), sep = "\t")
print(meta)
