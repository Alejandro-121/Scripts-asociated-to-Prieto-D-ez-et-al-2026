
# Scripts Associated to Prieto-Díez et al.

**Title:** Transcriptional reprogramming of the eIF5A silent paralogue gene bypasses mutations on the main eIF5A isoform
**Citation:** Prieto-Díez et al.
**Prepublished DOI:** In porcess

This repository contains all the bioinformatic scripts used in the above publication. The code is provided to facilitate reproducibility of the analyses described in the paper.

---

## Repository structure

```
.
├── GATK_haploid/          # Variant calling pipeline (haploid genomes)
├── parse-mutations-eif5a/ # Mutation parsing and genotype filtering
├── snpeff/                # VCF annotation with SnpEff
├── sppIDer/               # Genomic composition analysis in hybrid strains
└── grouth_rate/           # Curve analisys script
```

---

## GATK_haploid

GATK-based variant calling pipeline designed for haploid genomes, prepared for execution on HPC systems (tested on HPC-Drago). The workflow is split into two stages with a manual QC checkpoint in between.

Key scripts:

- `variant_calling_1.py` — alignment, duplicate marking, and pre-BQSR variant discovery.
- `variant_calling_2.py` — BQSR, final variant calling, and hard filtering.
- `variant_calling_fun.py` — helper functions used by the two scripts above.
- `slurm_variant_calling.sh` — SLURM wrapper for cluster submission.
- `plot_metrix.rmd` — R Markdown report to visualize QC metrics and guide filter thresholds.

Dependencies: GATK 4.x, BWA, SAMtools, Picard, R (`ggplot2`, `cowplot`, `gridExtra`).

---

## parse-mutations-eif5a

Pipeline to analyze variants from a consolidated VCF (`merged.vcf`) and generate per-sample mutation and allele tables, as well as genotype-based filters. Hardcoded to the samples used in the paper.

Key scripts:

- `01_variant_extraction.py` — reads `merged.vcf`, extracts `ANN` annotations, and outputs `mutation_wide_table_with_alleles.csv`.
- `02_genotype_based_filtering.R` — filters by genotype patterns, adds SGD gene annotations, and generates tables and plots.
- `run_pipeline.sh` — orchestrates both steps and validates dependencies.
- `environment.yml` — Conda environment definition (Python + R + Bioconductor packages).

Dependencies: Python 3.10+ (`pandas`), R 4.3+ (`readr`, `dplyr`, `ggplot2`, `pheatmap`, `org.Sc.sgd.db`, ...). See `environment.yml` for the full list.

---

## snpeff

Scripts to annotate VCF files using SnpEff against the *S. cerevisiae* R64-1-1 SGD database.

Key scripts:

- `launch_snpeff.sh` — loops over all VCF files in the directory, compresses and indexes them with `bgzip`/`tabix`, and annotates each one with SnpEff.
- `change_head.py` — renames FASTA headers of the S288C reference genome to standard chromosome names (`chrI`, `chrII`, …, `chrMito`) required by SnpEff.

Dependencies: SnpEff, bgzip, tabix (HTSlib), Python 3.

---

## sppIDer

Genomic composition analysis pipeline for detecting hybrid strains and species contributions from short-read sequencing data. Reads are mapped against a combined multi-species reference genome; coverage depth per species, chromosome, and sliding window is then used to identify hybrids, introgressions, and contamination.

This implementation is adapted for HPC-Drago and extends the original pipeline with a SLURM array launcher, automatic SE/PE detection, multi-run merging, and an aggregated HTML report across all samples. For full usage details see [`sppIDer/README.md`](sppIDer/README.md).

The pipeline is based on the original sppIDer tool developed by GLBRC. For methodology, citation, and upstream documentation see [https://github.com/GLBRC/sppIDer](https://github.com/GLBRC/sppIDer).

Key scripts:

- `run_sppIDer_array_se.sh` — SLURM array launcher; handles PE and SE samples, merges multi-run data, and distributes work across jobs.
- `sppIDer.py` — core pipeline per sample (BWA mapping → coverage → depth statistics → plots)
- `combineRefGenomes.py` — builds the combined multi-species reference FASTA and indexes.
- `aggregate_sppIDer_report.py` — run after all samples are processed; produces a summary TSV and an interactive HTML report with per-sample species calls and quality flags.

Dependencies: BWA, SAMtools, BEDTools, Python 3, R (`ggplot2`, `data.table`, `modes`).

This pipeline covers diferents uses in the asociated paper of this repo it was only used to verify the relative ploidy of the samples. 

---
## grouth rate
R script for calculating the maximum specific growth rate (μ) and doubling time of bacterial cultures from OD measurements in microplate format. For each biological replicate, the script identifies the optimal exponential growth window via exhaustive linear regression on log-transformed OD values and exports the results to Excel.

---

## Contact

For questions about the code, please open an issue in this repository.
