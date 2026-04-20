# 🧬 HERV RNA Regulatory Annotation Pipeline
This repository contains the code used to generate a genome-wide, multi-layer annotation of RNA regulatory features in human endogenous retroviruses (HERVs) in the GRCh38 human genome.

The pipeline integrates sequence-level, structural, and transcriptomic information to characterize the RNA-centered regulatory potential of HERV internal regions and long terminal repeats (LTRs).

## 🔍 Overview

This framework combines:

- RNA-binding protein binding site prediction (RBPBS)
- Transcript context integration (lncRNAs and protein-coding genes)
- Mapping of conserved retroviral domains to transcript features
- Identification of antisense LTR configurations (SPARCS-like, dsRNA potential)
- Functional annotation of HERV-associated lncRNAs via literature mining
- Interferon-related regulatory motif analysis (STAT / IRF)

Together, these analyses define a multi-layer regulatory landscape encoded by HERV sequences.

This repository corresponds to the computational framework supporting:

> Montserrat-Ayuso et al., Human endogenous retroviruses shape the RNA regulatory landscape of the human genome (manuscript in preparation)

## 📁 Repository structure
```
scripts/
├── 01_genome_annotation/
├── 02_rbp_motif_analysis/
├── 03_integrated_annotation/
├── 04_transcript_context/
├── 05_interferon_analysis/
├── 06_sparcs_analysis/
├── 07_literature_mining/
├── 08_literature_integration/
├── 09_overview/
```

## ⚙️ Repository description

### 01_genome_annotation

Preparation of genomic features and overlap calculations using BEDTools:

- GENCODE feature extraction (genes, exons, introns, UTRs)
- Overlaps between HERV regions, domains, LTRs, and transcript features
- Closest gene and TSS annotations

### 02_rbp_motif_analysis

Prediction of RNA-binding protein binding sites (RBPBS):

- Processing of ATtRACT motifs
- Background model generation
- FIMO motif scanning
- Conversion of motif hits to genomic coordinates

### 03_integrated_annotation

Construction of a unified annotation table:

- Integration of genomic, structural, and regulatory features
- Assignment of unique HERV locus identifiers
- Mapping of domains and LTR features to loci
- Generation of HERVarium-ready annotations

### 04_transcript_context

Integration of HERVs into transcript architectures:

- Overlap with lncRNAs and protein-coding genes
- Mapping of domains to exons, CDS, and UTRs
- Identification of terminal exon enrichment
- Annotation of HERV-associated lncRNAs

### 05_interferon_analysis

Analysis of interferon-responsive regulatory features:

- STAT and IRF motif enrichment
- Motif clustering within LTRs
- Classification of interferon regulatory potential

### 06_sparcs_analysis

Identification of SPARCS-like architectures:

- Detection of antisense LTRs in transcript termini
- Characterization of dsRNA-forming potential
- Integration with structural (U3–R–U5) LTR annotations

### 07_literature_mining

Automated extraction of lncRNA functional information:

- Retrieval of PubMed abstracts
- NLP-based extraction of functional annotations
- Classification of lncRNA-associated functions

### 08_literature_integration

Integration of curated and mined annotations:

- Merging literature-derived and database annotations
- Gene-level functional summaries

### 09_overview

High-level analysis and summary:

- ```01_analysis.R``` generates global summaries and figures used for interpretation

## 📦 Inputs

The pipeline relies on:

- GRCh38 genome assembly
- GENCODE gene annotation
- HERV internal regions and LTR annotations (from previous Zenodo releases)
- ATtRACT RBP motif database

## 📊 Outputs

Key outputs include:

- Genome-wide RBP binding site annotations
- Integrated HERV locus annotation tables
- Transcript-context mappings (lncRNAs, protein-coding genes)
- SPARCS-like candidate loci
- Interferon motif summaries
- Functionally annotated HERV-associated lncRNAs

These outputs are released separately via Zenodo: 

## 🧠 Conceptual framework

This pipeline implements a multi-layer model in which HERV sequences act as:

- RNA regulatory modules (RBP binding, RNA motifs)
- Structural components of transcripts (lncRNAs, UTRs, exons)
- Sources of residual coding potential (micropeptides)
- Drivers of dsRNA-mediated immune signaling (SPARCS-like architectures)

## 📖 Citation

If you use this code, please cite:

> Montserrat-Ayuso et al., Human endogenous retroviruses shape the RNA regulatory landscape of the human genome (manuscript in preparation)

and the corresponding Zenodo dataset: 

## 🚀 Future directions

This repository is part of HERVarium ([https://hervarium.cnag.eu/](https://hervarium.cnag.eu/)), an integrative framework for exploring HERV structure and function across multiple regulatory layers.
