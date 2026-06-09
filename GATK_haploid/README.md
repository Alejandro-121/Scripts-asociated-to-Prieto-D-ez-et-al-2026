# GATK Haploid Variant Calling Pipeline

This repository contains a two-stage GATK-based variant calling workflow designed
for haploid genomes. The pipeline is tailored for HPC execution and emphasizes a
manual QC checkpoint before base quality score recalibration (BQSR). All the pipeline is prepared to run in HPC-DRAGO

## Overview

The workflow is split into two Python entry points:

1. **`variant_calling_1.py`**: Alignment and pre-BQSR variant discovery.
2. **`variant_calling_2.py`**: BQSR, final variant calling, and hard filtering.

A SLURM wrapper (`slurm_variant_calling.sh`) is provided for cluster execution, and
an R Markdown report (`plot_metrix.rmd`) summarizes QC metrics to guide filtering
thresholds.

## Pipeline Stages

### Stage 1: Alignment and Pre-BQSR Discovery
Runs the following steps in order:

- Build reference indexes (`bwa`, `samtools`, `gatk CreateSequenceDictionary`).
- Align reads with `bwa mem` and sort/index alignments.
- Mark duplicates.
- Call preliminary variants with `HaplotypeCaller`.
- Split preliminary calls into SNPs and indels.

**Entry point:** `variant_calling_1.py`

### QC Checkpoint (Manual)
Review preliminary variant metrics to set appropriate hard-filter thresholds.
Use `plot_metrix.rmd` to plot distributions of FS, QD, MQ, SOR, MQRankSum,
ReadPosRankSum, and DP.

### Stage 2: BQSR and Final Variant Calls
Runs the following steps in order:

- Apply hard filters to pre-BQSR SNPs and indels.
- Select unfiltered variants as known-sites inputs for BQSR.
- Perform BQSR and generate covariate plots.
- Call final variants from recalibrated BAM.
- Split final calls into SNPs and indels.
- Apply hard filters to final SNPs and indels.

**Entry point:** `variant_calling_2.py`

## Requirements

- **GATK 4.x**
- **BWA**
- **SAMtools**
- **Picard** (only for `CollectAlignmentSummaryMetrics` in the SLURM wrapper)
- **R** with `ggplot2`, `cowplot`, and `gridExtra` (for QC plots)
- **bcftools** (used in `plot_metrix.rmd`)

On HPC systems, modules similar to those in `slurm_variant_calling.sh` are
expected.

## Inputs

- Reference genome FASTA (indexed during the pipeline).
- Illumina paired-end FASTQ files (or single-end with `read2` omitted).

## Outputs (Key Files)

- `*.sam`, `*sorted.bam`, `*sorted_dedup.bam`: alignment files.
- `*pre_bqsr_variants.vcf`: preliminary calls.
- `*bqsr_report.txt`, `*plots_post_bqsr_report.pdf`: BQSR outputs.
- `*bqsred_variants.vcf`: final calls from recalibrated BAM.
- `*bqsred_filtered_snps.vcf`, `*bqsred_filtered_indels.vcf`: final filtered calls.

## Usage

### Local Execution

```bash
python variant_calling_1.py \
  -p sample01 \
  -numth 8 \
  -r /path/to/reference.fasta \
  -r1 /path/to/sample01_R1.fastq.gz \
  -r2 /path/to/sample01_R2.fastq.gz
```

After reviewing QC metrics and selecting filter thresholds, run:
In the paper the default thresholds encoded in this pipeline were used, this reflets the recomeneded by GATK.

```bash
python variant_calling_2.py \
  -p sample01 \
  -r /path/to/reference.fasta
```

### HPC Execution (SLURM)

Edit `slurm_variant_calling.sh` with the correct sample inputs and submit:

```bash
sbatch slurm_variant_calling.sh
```
