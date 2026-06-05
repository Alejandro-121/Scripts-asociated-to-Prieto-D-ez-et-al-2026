# Parse Mutations Pipeline

This repository contains a pipeline to analyze variants from a consolidated VCF (`merged.vcf`) and generate per-sample mutation/allele tables, as well as filters based on genotype patterns. It is hardcoded to the data used in the paper. This repo is created to facilitate reproducibility

## Current pipeline status

- The pipeline is defined in `run_pipeline.sh` and runs two steps: variant extraction (Python) and genotype-based filtering (R).
- The expected input file is `merged.vcf` located at the project root. If it does not exist, the script fails at startup.

## Project structure

- `01_variant_extraction.py`: reads `merged.vcf`, extracts `ANN` annotations, and generates `mutation_wide_table_with_alleles.csv` with per-sample columns (mutation and alleles).
- `02_genotipe_based_filtering.R`: filters the table by genotype patterns, adds SGD gene annotations, and generates tables/plots of common mutations and reversions.
- `run_pipeline.sh`: orchestrates the steps and validates dependencies (Python/R + file existence).
- `environment.yml`: Conda environment definition with Python and R dependencies (includes Bioconductor packages).

## Requirements

- Python 3.10+ with `pandas`.
- R 4.3+ with packages: `readr`, `dplyr`, `ggplot2`, `reshape2`, `pheatmap`, `tidyr`, `org.Sc.sgd.db`, `AnnotationDbi`.
- The `merged.vcf` file at the repository root.
