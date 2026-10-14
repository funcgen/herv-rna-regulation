# Workflow

The analysis is a fixed sequence of steps in eight modules (`01`–`08`). This page gives the order, how to run each step
and the recorded HPC resources. Installation, configuration and input data are described in the
[README](../README.md).

## Module dependencies

```text
01_sources ──┬─> 02_rbp ─────────────────────────────────────────────┐
             ├─> 03_transcript_context ─┬─> 04_domains_terminal_exons ├─> 08_figures_supplement
             │                          ├─> 05_lncrna ─> 06_lncrna_function
             │                          └─> 07_sparcs ────────────────┘
             └─> 07_sparcs (LTR projection, step 16, is read by 03 step 17)
```

## Running the steps

- **HPC steps** have a SLURM wrapper (`*.slurm`) next to the producer. Replace the placeholders first (see the README),
  then submit from the wrapper's directory:

  ```bash
  sbatch --chdir=03_transcript_context/scripts/overlaps run_overlap_internal_vs_features.slurm
  ```

- **Local steps** are run by the `run_local.sh` of their module, one step or all steps in order:

  ```bash
  05_lncrna/run_local.sh                      # all local steps of the module, in order
  05_lncrna/run_local.sh coverage_categories  # one step
  ```

Most steps check that their inputs match the checksums recorded by the step that produced them (or the published ones)
and stop otherwise, so steps must run in the order below.

## Execution order

`Where`: `HPC` = SLURM cluster, `LOCAL` = workstation.

| Order | Module | Step | Where | How to run |
|---|---|---|---|---|
| 1 | `01_sources` | Fetch and verify published upstream releases | HPC | `01_sources/scripts/fetch_upstream_releases.sh`, then `extract_ltr_fimo.sh` (or its `.slurm`) and `verify_sources.sh` |
| 2 | `01_sources` | Register manual snapshots | LOCAL | not included: registers the dated snapshots from the original project |
| 3 | `01_sources` | Build HERV locus identifiers | HPC | not needed (the HERV GTF is published in Zenodo 19661036) |
| 4 | `02_rbp` | Build RBP motif library | HPC | `02_rbp/scripts/motif_library/build_attract_meme.slurm` |
| 5 | `02_rbp` | Estimate strand-specific Markov background | HPC | `02_rbp/scripts/scanning/make_markov_background.slurm` |
| 6 | `02_rbp` | Primary FIMO scan | HPC | `02_rbp/scripts/scanning/run_fimo_primary.slurm` |
| 7 | `02_rbp` | Parse motif hits to genomic coordinates | HPC | `02_rbp/scripts/scanning/run_parse_fimo_primary.slurm` |
| 8 | `03_transcript_context` | Extract GENCODE features | HPC | `03_transcript_context/scripts/gencode_features/run_build_feature_beds.slurm` |
| 9 | `03_transcript_context` | Build per-gene merged exons and per-gene introns | HPC | `03_transcript_context/scripts/gencode_features/run_build_introns.slurm` |
| 10 | `03_transcript_context` | Sort assets in genome order | HPC | `03_transcript_context/scripts/gencode_features/run_sort_assets.slurm` |
| 11 | `03_transcript_context` | Define transcript-aware UTR split and terminal exons | HPC | `03_transcript_context/scripts/gencode_features/run_build_terminal_exons.slurm` |
| 12 | `03_transcript_context` | Overlap HERV internal regions with gene features | HPC | `03_transcript_context/scripts/overlaps/run_overlap_internal_vs_features.slurm` |
| 13 | `03_transcript_context` | Nearest gene and TSS | HPC | `03_transcript_context/scripts/overlaps/run_closest_gene_tss.slurm` |
| 14 | `03_transcript_context` | Overlap domains with gene and transcript features | HPC | `03_transcript_context/scripts/overlaps/run_overlap_domains_vs_features.slurm` |
| 15 | `03_transcript_context` | Overlap LTRs with transcript features | HPC | `03_transcript_context/scripts/overlaps/run_overlap_ltrs_vs_features.slurm` |
| 16 | `07_sparcs` | Build the minimal LTR feature projection | HPC | `07_sparcs/scripts/ltr_structure/run_build_ltr_projection.slurm` |
| 17 | `03_transcript_context` | Build integrated locus-level table | HPC | `03_transcript_context/scripts/integration/run_build_integrated_locus_table.slurm` |
| 18 | `03_transcript_context` | Attach HERV locus identifiers | HPC | `03_transcript_context/scripts/integration/run_add_herv_locus_id.slurm` |
| 19 | `03_transcript_context` | Build HERV x gene and transcript overlap summaries | LOCAL | `03_transcript_context/run_local.sh gene_overlap_tables` |
| 20 | `02_rbp` | Generate dinucleotide-preserving shuffles | HPC | `02_rbp/scripts/null_model/generate_dinuc_shuffles.slurm` |
| 21 | `02_rbp` | QC the shuffles | HPC | `02_rbp/scripts/null_model/run_qc_dinuc_shuffles.sh` (run on the HPC node, after step 20) |
| 22 | `02_rbp` | FIMO scan of observed plus five shuffles | HPC | `02_rbp/scripts/null_model/run_fimo_null_array.slurm` |
| 23 | `02_rbp` | Global observed versus null summary | HPC | `02_rbp/scripts/enrichment/run_summarise_global_null.slurm` |
| 24 | `02_rbp` | Subfamily composition-controlled enrichment | HPC | `02_rbp/scripts/enrichment/run_subfamily_enrichment.slurm` |
| 25 | `02_rbp` | HERVH/HERVK profiles, GO ORA and pre-ranked GSEA | HPC | `02_rbp/scripts/enrichment/run_hervh_hervk_go_gsea.slurm` |
| 26 | `02_rbp` | Build strict-threshold datasets | HPC | `02_rbp/scripts/sensitivity_strict/run_filter_fimo_strict.slurm` |
| 27 | `02_rbp` | Strict global, subfamily and GO/GSEA re-analysis | HPC | `02_rbp/scripts/sensitivity_strict/run_strict_analysis.slurm` (steps 23-25 at the strict threshold) |
| 28 | `02_rbp` | Threshold comparison and sensitivity summary | HPC | `02_rbp/scripts/sensitivity_strict/run_compare_thresholds.slurm` |
| 29 | `02_rbp` | FMR1 locus screen | LOCAL | `02_rbp/run_local.sh fmr1_locus_screen` |
| 30 | `04_domains_terminal_exons` | Build domain-centred transcript table | HPC | `04_domains_terminal_exons/scripts/domain_transcript_table/run_build_domain_table.slurm` |
| 31 | `04_domains_terminal_exons` | Transcript-aware positional null | HPC | `04_domains_terminal_exons/scripts/positional_null/run_terminal_exon_null.slurm` |
| 32 | `04_domains_terminal_exons` | Conserved versus lower-coverage comparison | HPC | `04_domains_terminal_exons/scripts/conservation_comparison/run_compare_conservation.slurm` |
| 33 | `04_domains_terminal_exons` | Locus-representative sensitivity | HPC | `04_domains_terminal_exons/scripts/conservation_comparison/run_compare_locus_representative.slurm` |
| 34 | `04_domains_terminal_exons` | Terminal-exon architecture | HPC | `04_domains_terminal_exons/scripts/architecture/run_terminal_exon_architecture.slurm` |
| 35 | `05_lncrna` | Per-gene HERV coverage for lncRNAs | LOCAL | `05_lncrna/run_local.sh lncrna_coverage` |
| 36 | `05_lncrna` | Per-gene HERV coverage for protein-coding genes | LOCAL | `05_lncrna/run_local.sh pc_coverage` |
| 37 | `05_lncrna` | Assign HERV coverage categories | LOCAL | `05_lncrna/run_local.sh coverage_categories` |
| 38 | `05_lncrna` | Annotate lncRNA overlap tables | LOCAL | `05_lncrna/run_local.sh lncrna_annotation` |
| 39 | `06_lncrna_function` | Retrieve PubMed abstracts | LOCAL | `06_lncrna_function/run_local.sh pubmed_corpus` |
| 40 | `06_lncrna_function` | Extract lncRNA function statements | LOCAL | `06_lncrna_function/run_local.sh function_statements` |
| 41 | `06_lncrna_function` | Classify lncRNAs by evidence | LOCAL | `06_lncrna_function/run_local.sh literature_classification` |
| 42 | `06_lncrna_function` | Integrate curated, GO and literature annotation | LOCAL | `06_lncrna_function/run_local.sh integrated_annotation` |
| 43 | `07_sparcs` | Build antisense terminal LTR candidates | HPC | `07_sparcs/scripts/candidates/run_build_sparcs_candidates.slurm` |
| 44 | `07_sparcs` | Process interferon STAT and IRF motifs | LOCAL | `07_sparcs/run_local.sh interferon_motifs` |
| 45 | `07_sparcs` | SPARCS-like analysis, GO, MSigDB and TF motifs | LOCAL | `07_sparcs/run_local.sh enrichment` |
| 46 | `07_sparcs` | U3-R-U5 structural comparison | LOCAL | `07_sparcs/run_local.sh u3ru5_integrity` |
| 47 | `08_figures_supplement` | Plot Figure 1 panels | LOCAL | `08_figures_supplement/run_local.sh figure1` |
| 48 | `08_figures_supplement` | Plot Figure 2 panels | LOCAL | `08_figures_supplement/run_local.sh figure2` |
| 49 | `08_figures_supplement` | Plot Figure 3 panels | LOCAL | `08_figures_supplement/run_local.sh figure3` |
| 50 | `08_figures_supplement` | Plot Figure 4 | LOCAL | `08_figures_supplement/run_local.sh figure4` |
| 51 | `08_figures_supplement` | Plot Figure 5 panels | LOCAL | `08_figures_supplement/run_local.sh figure5` |
| 54 | `08_figures_supplement` | Build supplement workbooks and tables (Supplemental Tables 1-8) | LOCAL | `08_figures_supplement/run_local.sh supplement` |
| 54.1 | `08_figures_supplement` | Plot Supplemental Figure 2 | LOCAL | `08_figures_supplement/run_local.sh supplementary_figure2` |

Not included: the HERVarium browser panels (Figures 1C, 2E, 3C), the Figure 5A schematic, Supplemental Figure 1 and the
final multi-panel layouts, which are manual assets; and the release validation of the original project.

## Recorded HPC resources

Requests in the `#SBATCH` headers of the wrappers. Most steps finish in minutes; the largest requests are the primary
FIMO scan (50 GB) and the subfamily enrichment (32 GB).

| SLURM wrapper | CPUs | Memory | Time limit | Array |
|---|---|---|---|---|
| `01_sources/scripts/extract_ltr_fimo.slurm` | 2 | 8G | 02:00:00 | — |
| `02_rbp/scripts/enrichment/run_hervh_hervk_go_gsea.slurm` | 1 | 8G | 01:00:00 | — |
| `02_rbp/scripts/enrichment/run_subfamily_enrichment.slurm` | 1 | 32G | 02:00:00 | — |
| `02_rbp/scripts/enrichment/run_summarise_global_null.slurm` | 1 | 16G | 00:30:00 | — |
| `02_rbp/scripts/motif_library/build_attract_meme.slurm` | 1 | 4G | 00:20:00 | — |
| `02_rbp/scripts/null_model/generate_dinuc_shuffles.slurm` | 1 | 8G | 01:00:00 | 1-5 |
| `02_rbp/scripts/null_model/run_fimo_null_array.slurm` | 1 | 8G | 02:00:00 | 1-5 |
| `02_rbp/scripts/scanning/make_markov_background.slurm` | 1 | 8G | 00:30:00 | — |
| `02_rbp/scripts/scanning/run_fimo_primary.slurm` | 1 | 50G | 02:00:00 | — |
| `02_rbp/scripts/scanning/run_parse_fimo_primary.slurm` | 1 | 8G | 00:20:00 | — |
| `02_rbp/scripts/sensitivity_strict/run_compare_thresholds.slurm` | 1 | 2G | 00:30:00 | — |
| `02_rbp/scripts/sensitivity_strict/run_filter_fimo_strict.slurm` | 1 | 2G | 00:30:00 | — |
| `02_rbp/scripts/sensitivity_strict/run_strict_analysis.slurm` | 1 | 8G | 01:00:00 | — |
| `03_transcript_context/scripts/gencode_features/run_build_feature_beds.slurm` | 1 | 4G | 01:00:00 | — |
| `03_transcript_context/scripts/gencode_features/run_build_introns.slurm` | 1 | 8G | 01:00:00 | — |
| `03_transcript_context/scripts/gencode_features/run_build_terminal_exons.slurm` | 1 | 8G | 01:00:00 | — |
| `03_transcript_context/scripts/gencode_features/run_build_tss.slurm` | 1 | 4G | 00:30:00 | — |
| `03_transcript_context/scripts/gencode_features/run_sort_assets.slurm` | 1 | 16G | 01:00:00 | — |
| `03_transcript_context/scripts/integration/run_add_herv_locus_id.slurm` | 1 | 2G | 00:15:00 | — |
| `03_transcript_context/scripts/integration/run_build_integrated_locus_table.slurm` | 1 | 16G | 01:00:00 | — |
| `03_transcript_context/scripts/overlaps/run_closest_gene_tss.slurm` | 1 | 4G | 00:30:00 | — |
| `03_transcript_context/scripts/overlaps/run_overlap_domains_vs_features.slurm` | 1 | 4G | 00:30:00 | — |
| `03_transcript_context/scripts/overlaps/run_overlap_internal_vs_features.slurm` | 1 | 8G | 01:00:00 | — |
| `03_transcript_context/scripts/overlaps/run_overlap_ltrs_vs_features.slurm` | 1 | 4G | 00:30:00 | — |
| `04_domains_terminal_exons/scripts/architecture/run_terminal_exon_architecture.slurm` | 1 | 16G | 01:00:00 | — |
| `04_domains_terminal_exons/scripts/conservation_comparison/run_compare_conservation.slurm` | 1 | 2G | 00:10:00 | — |
| `04_domains_terminal_exons/scripts/conservation_comparison/run_compare_locus_representative.slurm` | 1 | 2G | 00:10:00 | — |
| `04_domains_terminal_exons/scripts/domain_transcript_table/run_build_domain_table.slurm` | 1 | 4G | 00:20:00 | — |
| `04_domains_terminal_exons/scripts/positional_null/run_terminal_exon_null.slurm` | 1 | 4G | 00:30:00 | — |
| `07_sparcs/scripts/candidates/run_build_sparcs_candidates.slurm` | 1 | 8G | 00:30:00 | — |
| `07_sparcs/scripts/ltr_structure/run_build_ltr_projection.slurm` | 1 | 8G | 01:00:00 | — |

## Outputs

Each module writes to its own `results/` directory (`08_figures_supplement` to `figures/`, `panel_sources/` and
`supplement/tables/`). `RESULTS_MANIFEST.tsv` in a module (for `02_rbp`, `checksums.sha256`) lists the SHA-256 of every
published output, so a re-run can be compared with the published results byte for byte.
