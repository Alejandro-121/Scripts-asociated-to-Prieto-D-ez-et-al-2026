# Scripts Associated to Prieto-Díez et al. XXX
> **Title** XXX
> **Citation:** Prieto-Díez et al. (XXX) Prepublished
> **DOI:** Prepublished

This repository contains all the bioinformatic scripts used in the above publication. The code is provided to facilitate reproducibility of the analyses described in the paper.

---

## Repository structure

```
.
├── GATK_haploid/          # Variant calling pipeline (haploid genomes)
├── parse-mutations-eif5a/ # Mutation parsing and genotype filtering
└── snpeff/                # VCF annotation with SnpEff
```

---

## GATK_haploid

GATK-based variant calling pipeline designed for **haploid genomes**, prepared for execution on HPC systems (tested on HPC Drago). The workflow is split into two stages with a manual QC checkpoint in between.

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

Scripts to annotate VCF files using **SnpEff** against the *S. cerevisiae* R64-1-1 SGD database.

Key scripts:
- `launch_snpeff.sh` — loops over all VCF files in the directory, compresses and indexes them with `bgzip`/`tabix`, and annotates each one with SnpEff.
- `change_head.py` — renames FASTA headers of the S288C reference genome to standard chromosome names (`chrI`, `chrII`, …, `chrMito`) required by SnpEff.

Dependencies: SnpEff, bgzip, tabix (HTSlib), Python 3.

---

## Contact

For questions about the code, please open an issue in this repository.
