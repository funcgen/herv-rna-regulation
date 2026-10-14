# herv-rna-regulation

Scripts to reproduce the analyses of **"Genome-wide characterization of HERV-derived RNA regulatory potential in the
human transcriptome"** (Montserrat-Ayuso T, Pujol A, Esteve-Codina A; manuscript in preparation).

Version 1.0.0 · License: [MIT](LICENSE)

The workflow characterises the RNA regulatory potential of 53,176 human endogenous retrovirus (HERV) internal regions
and 560,594 HERV long terminal repeats (LTRs) in GRCh38 (GENCODE release 48). All results are sequence- and
annotation-based predictions.

## Contents

| Directory | Analysis |
|---|---|
| `01_sources/` | Download and checksum verification of the inputs |
| `02_rbp/` | RNA-binding protein motif scans (FIMO; 151 ATtRACT/CISBP-RNA motifs, 93 RBPs), dinucleotide-shuffle null model, composition-controlled subfamily enrichment, HERVH/HERVK GO and GSEA, strict-threshold sensitivity, FMR1 locus screen |
| `03_transcript_context/` | GENCODE v48 features; overlap of HERV internal regions, domains and LTRs with transcript features; integrated locus table |
| `04_domains_terminal_exons/` | Conserved retroviral domains in transcripts, transcript-aware positional null model, terminal-exon architecture |
| `05_lncrna/` | HERV coverage of lncRNA and protein-coding genes, coverage categories, LncRNAWiki/GO annotation |
| `06_lncrna_function/` | PubMed literature mining and functional annotation of HERV-associated lncRNAs |
| `07_sparcs/` | SPARCS-like antisense LTR configurations; GO, MSigDB C7, interferon (STAT/IRF) and TF motif analyses; U3–R–U5 comparison |
| `08_figures_supplement/` | Figure panels, Supplemental Figure 2 and Supplemental Tables 1–8 from the module results |
| `config/` | `parameters.yaml` (every scientific threshold and seed), `paths.yaml`, `lib/common.sh`, `environments/` (software versions) |
| `docs/WORKFLOW.md` | Step order, how to run each step, HPC resources |

In each module, `scripts/` holds the producers and the SLURM wrappers (`*.slurm`) of the HPC steps; `run_local.sh` runs
the module's workstation steps.

## Requirements

- Linux; a SLURM cluster for the large steps (up to 50 GB of memory for one job); about 12 GB of disk for inputs.
- `bash`, `curl`, `xz`, `sha256sum`, `/usr/bin/time` and `unshare` (util-linux), used by the run scripts.
- Software, at the versions the results were produced with (specifications in [`config/environments/`](config/environments/)):

| Software | Used by | Versions |
|---|---|---|
| MEME Suite | `02` | 5.5.9 (`fasta-get-markov`, `fimo`) |
| bedtools | `03`, `07` | 2.31.1 |
| Python | `03`, `04`, `05`, `07` | 3.11 with pandas 2.0.3, numpy 1.26.4, scipy 1.17.1, PyYAML |
| R | `02`, `03`, `05`–`08` | 4.3.3 with clusterProfiler 4.10.1, enrichplot 1.22.0, fgsea 1.28.0, org.Hs.eg.db 3.18.0, GO.db 3.18.0, AnnotationDbi 1.64.1, DOSE, data.table, tidyverse (readr, dplyr, tidyr, stringr, purrr, tibble, forcats), scales, digest, jsonlite, yaml, zip; ggplot2 4.0.2, patchwork 1.3.2, cowplot, svglite 2.1.3, ragg 1.2.6, systemfonts, pheatmap 1.0.12 and openxlsx 4.2.5 for `08`; msigdbr 26.1.0 for the MSigDB snapshot |
| Python NLP | `06` | Python 3.11.16, spaCy 3.8.11, en_core_web_sm 3.8.0 (full pin list: `config/environments/python_nlp.lock.txt`) |

The org.Hs.eg.db and GO.db versions define the GO gene sets; other versions can change the enrichment results.

## Installation and configuration

```bash
git clone https://github.com/funcgen/herv-rna-regulation.git
cd herv-rna-regulation
export NCRNA_RELEASE_ROOT="$PWD"
export NCRNA_SYSTEM=cnag                          # the scripts' profile name for a full run (all inputs; see config/paths.yaml)
export NCRNA_GENCODE_V48_GTF=/path/to/gencode.v48.primary_assembly.annotation.gtf
export NCRNA_GRCH38_FA=/path/to/GRCh38.primary_assembly.genome.fa
export NCRNA_RSCRIPT=/path/to/Rscript            # R environment above (default: Rscript on PATH)
export NCRNA_MEME_BIN=/path/to/meme-5.5.9/bin
export NCRNA_BEDTOOLS=/path/to/bedtools           # 2.31.1
export NCRNA_NLP_PYTHON=/path/to/nlp-env/bin/python
```

The SLURM wrappers contain placeholders for site-specific values; replace them before submitting:
`your_slurm_account`, `your_partition`, `/path/to/herv-rna-regulation` (log paths), `/path/to/conda/envs/<env>`
(interpreters and tools) and `/path/to/references/` (GRCh38 and GENCODE files).

## Input data

No data are stored in this repository.

- **GRCh38 primary assembly and GENCODE release 48 annotation**:
  <https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/>.
- Every input, with its version, source and checksum, is listed in `01_sources/SOURCE_STATUS.tsv`.
- **HERV annotations** (internal regions, domains, LTRs, U3–R–U5 catalogue, LTR motif scan, ATtRACT files, HERV GTF), from
  the authors' Zenodo records [10.5281/zenodo.17662456](https://doi.org/10.5281/zenodo.17662456),
  [10.5281/zenodo.19593457](https://doi.org/10.5281/zenodo.19593457) and
  [10.5281/zenodo.19661036](https://doi.org/10.5281/zenodo.19661036). Download and verify them with:

  ```bash
  01_sources/scripts/fetch_upstream_releases.sh
  01_sources/scripts/extract_ltr_fimo.sh
  01_sources/scripts/verify_sources.sh
  ```

- **Dated snapshots** of LncRNAWiki (2026-03-10), the PubMed abstract corpus (2026-03-11), an Ensembl BioMart GO table
  (2026-04-16) and MSigDB C7 (2026.1.Hs; regenerated by `01_sources/scripts/freeze_msigdb_c7.R`), read by `05_lncrna`
  (annotation step), `06_lncrna_function` and `07_sparcs` (C7 analysis). They are not redistributed; their file names and
  checksums are in `01_sources/manifests/manual_snapshot_checksums.tsv`, and they are expected under `01_sources/frozen/`.

## Running the workflow

Run the steps in the order given in [docs/WORKFLOW.md](docs/WORKFLOW.md): HPC steps with their SLURM wrapper, local
steps with the module's `run_local.sh`. For example:

```bash
sbatch --chdir=02_rbp/scripts/scanning run_fimo_primary.slurm
05_lncrna/run_local.sh
08_figures_supplement/run_local.sh supplement
```

Most steps verify their inputs against the checksums recorded by the step that produced them and stop if they differ.

## Outputs

Each module writes to its `results/` directory; `08_figures_supplement` writes the figure panels and their source data
(`figures/`, `panel_sources/`) and Supplemental Tables 1–8 (`supplement/tables/`). `RESULTS_MANIFEST.tsv` in each module
(for `02_rbp`, `checksums.sha256`) lists the SHA-256 of every published output, so a re-run can be compared with the
published results.

Not produced by scripts: the HERVarium browser panels (Figures 1C, 2E, 3C), the Figure 5A schematic, Supplemental
Figure 1 and the final multi-panel layouts of the figures.

## Validation

Before release, every value reported in the manuscript was checked against the outputs of this workflow, and repeated
runs of each step produced byte-identical outputs.

## Citation, license and contact

If you use this code, please cite:

Montserrat-Ayuso et al., Genome-wide characterization of HERV-derived RNA regulatory potential in the
human transcriptome (manuscript in preparation)

and the corresponding Zenodo dataset: [https://doi.org/10.5281/zenodo.23361077](https://doi.org/10.5281/zenodo.23361077)

Contact: Anna Esteve-Codina, Centre Nacional d'Anàlisi Genòmica (CNAG), Barcelona (anna.esteve@cnag.eu).

This repository is part of HERVarium (https://hervarium.cnag.eu/), an integrative framework for exploring HERV structure and function across multiple regulatory layers.

### Funding

This work was supported by the AGAUR-FI predoctoral grant programme (2025 FI-1 00642) Joan Oró of the Secretariat for
Universities and Research, Department of Research and Universities, Government of Catalonia, and by the European Social
Fund Plus. AP received funding from the Secretariat for Universities and Research of the Government of Catalonia
(2021SGR00899) and the Instituto de Salud Carlos III (FIS PI23/00835; FEDER, "Una manera de hacer Europa").
Institutional support to CNAG came from the Spanish Ministry of Science, Innovation and Universities, the Fondo de
Investigaciones Sanitarias co-funded with ERDF funds (PI19/01772), the 2014–2020 Smart Growth Operating Program, and the
Generalitat de Catalunya through the Departament de Recerca i Universitats and the Departament de Salut.
