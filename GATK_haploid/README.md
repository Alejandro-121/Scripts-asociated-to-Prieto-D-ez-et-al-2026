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

### SLURM modules

All bioinformatics tools are loaded via `module load` in `slurm_variant_calling.sh`. The modules currently configured for HPC-Drago are:

| Module | Version | Role |
|--------|---------|------|
| GCCcore | 11.2.0 | Dependency for BWA, GATK, and R |
| BWA | 0.7.17 | Read alignment |
| GATK | 4.2.5.0-Java-11 | Variant calling, BQSR, filtering |
| Picard | 2.25.1-Java-11 | `CollectAlignmentSummaryMetrics` (called directly in the SLURM script, not from Python — see note below) |
| GCC | 11.2.0 | Dependency for SAMtools |
| SAMtools | 1.14 | BAM sorting, indexing, filtering |
| OpenMPI | 4.1.1 | Dependency for R |
| R | 4.1.2 | QC plots (`plot_metrix.rmd`) |

> ⚠️ **Picard note:** `CollectAlignmentSummaryMetrics` cannot be called from a Python subprocess on Drago due to the Java environment path. It is called directly in `slurm_variant_calling.sh` between Stage 1 and Stage 2 using `java -jar $EBROOTPICARD/picard.jar`.

> ⚠️ **Updating modules:** Replace the module names above with those available on your cluster:
> ```bash
> sed -i 's|module load GCCcore/11.2.0|module load <your_GCCcore_module>|g'  slurm_variant_calling.sh
> sed -i 's|module load BWA/0.7.17|module load <your_BWA_module>|g'           slurm_variant_calling.sh
> sed -i 's|module load GATK/4.2.5.0-Java-11|module load <your_GATK_module>|g' slurm_variant_calling.sh
> sed -i 's|module load picard/2.25.1-Java-11|module load <your_picard_module>|g' slurm_variant_calling.sh
> sed -i 's|module load GCC/11.2.0|module load <your_GCC_module>|g'           slurm_variant_calling.sh
> sed -i 's|module load SAMtools/1.14|module load <your_SAMtools_module>|g'    slurm_variant_calling.sh
> sed -i 's|module load OpenMPI/4.1.1|module load <your_OpenMPI_module>|g'     slurm_variant_calling.sh
> sed -i 's|module load R/4.1.2|module load <your_R_module>|g'                 slurm_variant_calling.sh
> ```

### R packages

Required by `plot_metrix.rmd` and must be available in the R module loaded above:

| Package | Role |
|---------|------|
| `ggplot2` | QC metric plots |
| `cowplot` | Plot composition |
| `gridExtra` | Multi-panel layout |
| `bcftools` | Called as a system command from within the Rmd |

### Python dependencies

The pipeline scripts (`variant_calling_1.py`, `variant_calling_2.py`, `variant_calling_fun.py`) use **only Python standard library modules** — no additional packages need to be installed:

| Module | Role |
|--------|------|
| `subprocess` | Calls all external tools (bwa, gatk, samtools) |
| `os` | Working directory handling |
| `argparse` | Command-line argument parsing |

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

Note: In the paper, the default thresholds encoded in this pipeline were used; this reflects those recommended by GATK.

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
