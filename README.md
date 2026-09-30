
# Scripts Associated to Prieto-Díez et al.

**Title:** Transcriptional reprogramming of the eIF5A silent paralogue gene bypasses mutations on the main eIF5A isoform

**Citation:** Prieto-Díez et al.

**Prepublished DOI:** In porcess

This repository contains all the bioinformatic scripts used in the above publication. The code is provided to facilitate reproducibility of the analyses described in the paper.

---

## Repository structure

```
.
├── wgs_variant_calling/   # WGS pipeline: download, trimming, mapping, BQSR, joint variant calling (SLURM)
├── snpeff/                # Custom SnpEff database (R64-1-1) and VCF annotation
├── suppressor_variants/   # Suppressor-specific variants, Supplementary Tables S6-S7 and Figure 3a
├── sppider/               # Genomic composition / ploidy analysis
└── growth_rate/           # Growth curve analysis
```

---

## WGS workflow

`wgs_variant_calling` → `snpeff` → `suppressor_variants`

1. **`wgs_variant_calling/`** — SLURM pipeline for the 11 sequenced strains (wild-type BY4741, parental *tif51A-1* and *tif51A-3*, and eight suppressors): SRA download, Trim Galore, BWA-MEM2 + MarkDuplicates, bootstrapped BQSR, and GATK HaplotypeCaller in per-base GVCF mode (haploid) followed by joint genotyping with GenotypeGVCFs and hard filtering. Joint genotyping gives every strain a genotype and its read depth at every variant site. See [`wgs_variant_calling/README.md`](wgs_variant_calling/README.md).
2. **`snpeff/`** — builds a SnpEff database from the S288C R64-1-1 reference and SGD annotation (validated against SGD coding and protein sequences) and annotates the joint VCF. See [`snpeff/README.md`](snpeff/README.md).
3. **`suppressor_variants/`** — defines suppressor-specific variants by comparison with the three sequenced control strains (a control must be genotyped as reference with sufficient coverage; no-calls are never treated as reference), and produces Supplementary Table S6 (variants), Supplementary Table S7 (NCBI accessions) and Figure 3a. See [`suppressor_variants/README.md`](suppressor_variants/README.md).

Data: BioProject PRJNA1418127 (SRA runs SRR37083317–SRR37083327).

> This workflow replaces the per-sample variant calling and genotype-parsing scripts of the first submission (`GATK_haploid/`, `parse-mutations-eif5a/`, `snpeff/launch_snpeff.sh`, `snpeff/change_head.py`), which remain available in the git history.

---

## sppIDer

Genomic composition analysis pipeline for detecting hybrid strains and species contributions from short-read sequencing data. Reads are mapped against a combined multi-species reference genome; coverage depth per species, chromosome, and sliding window is then used to identify hybrids, introgressions, and contamination.

This implementation is adapted for HPC-Drago and extends the original pipeline with a SLURM array launcher, automatic SE/PE detection, multi-run merging, and an aggregated HTML report across all samples. For full usage details see [`sppIDer/README.md`](sppIDer/README.md).

The pipeline is based on the original sppIDer tool developed by GLBRC. For methodology, citation, and upstream documentation see [https://github.com/GLBRC/sppIDer](https://github.com/GLBRC/sppIDer).

Key scripts:

- `run_sppIDer_array+se.sh` — SLURM array launcher; handles PE and SE samples, merges multi-run data, and distributes work across jobs.
- `sppIDer.py` — core pipeline per sample (BWA mapping → coverage → depth statistics → plots)
- `combineRefGenomes.py` — builds the combined multi-species reference FASTA and indexes.
- `aggregate_sppIDer_report.py` — run after all samples are processed; produces a summary TSV and an interactive HTML report with per-sample species calls and quality flags.

Dependencies: BWA, SAMtools, BEDTools, Python 3, R (`ggplot2`, `data.table`, `modes`).

This pipeline covers diferents uses in the asociated paper of this repo it was only used to verify the relative ploidy of the samples. 

---
## growth rate
R script for calculating the maximum specific growth rate (μ) and doubling time of cultures from OD measurements in microplate format. For each biological replicate, the script identifies the optimal exponential growth window via exhaustive linear regression on log-transformed OD values and exports the results to Excel.

---

## Contact

For questions about the code, please open an issue in this repository.
