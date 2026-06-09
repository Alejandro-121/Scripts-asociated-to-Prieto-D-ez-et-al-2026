# sppIDer — Usage Guide

> **Credits:** sppIDer was originally developed by the Gasch Lab at the Great Lakes Bioenergy Research Center (GLBRC). For full details on the pipeline design, methodology, and citation, please refer to the original repository and documentation at [https://github.com/GLBRC/sppIDer](https://github.com/GLBRC/sppIDer).

> **HPC-Drago:** This version of the pipeline is prepared to run on the HPC-Drago cluster. All paths, module loads, and SLURM configuration are set up for that environment.

---

## What this pipeline produces

sppIDer takes short-read Illumina data from a sample of unknown or mixed genomic origin and maps it against a **combined reference genome** built from multiple species. The goal is to determine how much of each reference species is represented in the sample — detecting hybrids, introgressions, contamination, or aneuploidies.

For each sample the pipeline produces:

- **Mapping quality report** — percentage of reads mapped per species, mapping quality (MQ) score distributions, and a chi-squared test comparing species contributions. Output: `_plotMQ.pdf`, `_MQsummary.txt`.
- **Species-level coverage** — mean sequencing depth per species across the combined genome, expressed both in absolute depth and as log2 relative to the whole-genome mean. Output: `_speciesDepth.pdf`, `_speciesAvgDepth-*.txt`.
- **Window-level coverage** — the genome is tiled into sliding windows; each window gets a mean depth and a log2 coverage value. This is the core data used to visualise which genomic regions belong to which species. Output: `_winAvgDepth-*.txt`.
- **Coverage bin analysis** — windows are grouped into discrete coverage bins by detecting anti-modes in the log2 depth distribution. Each species is flagged as "contributing" if enough of its windows fall above the baseline bin, which filters out noise from non-contributing references. Output: `_covDistPlots.pdf`, `_sppIDerDepthPlot-covBins.pdf`, `_covBinsSummary.txt`.

A mitochondrial variant (`mitoSppIDer`) runs the same logic on a mitochondrial combined reference, with optional CDS region shading on the coverage plots when a combined GFF file is supplied.

---

## Running the pipeline on a SLURM cluster: `run_sppIDer_array_se.sh`

This is the **recommended launcher** for processing multiple samples on the Drago cluster. It handles paired-end (PE) and single-end (SE) samples in the same run, merges samples that were sequenced across multiple runs, and distributes work across a SLURM array job.

### Basic usage

```bash
bash run_sppIDer_array_se.sh \
    -i /path/to/reads/ \
    -r /path/to/combined/SaccharomycesCombo.fasta \
    -o /path/to/output/
```

All arguments are optional — defaults are hardcoded at the top of the script and should be edited to match your paths before first use.

### Options

| Flag | Description | Default |
|------|-------------|---------|
| `-i` | Directory containing trimmed read files | hardcoded path |
| `-r` | Path to the combined reference FASTA | hardcoded path |
| `-o` | Output directory (created if absent) | hardcoded path |
| `-j` | Maximum number of SLURM array jobs | 60 |
| `-s` | R1 file suffix (e.g. `_R1.trimmed.fastq.gz`) | auto-detected |
| `-S` | R2 file suffix (e.g. `_R2.trimmed.fastq.gz`) | inferred from `-s` |
| `-b` | Pass `--byGroup` to sppIDer (faster for large genomes) | off |
| `-k` | Pass `--keep-intermediates` to sppIDer (retain SAM/BAM) | off |

### What the script does automatically

**Suffix auto-detection.** If `-s` is not provided, the script scans the reads directory and tries a prioritised list of known PE suffixes (`_R1.trimmed.fastq.gz`, `_1.trimmed.fastq.gz`, `_R1.fastq.gz`, `_1.fastq.gz`) until it finds matching files.

**SE/PE detection.** For every R1 file found, the script looks for the corresponding R2. If R2 is missing the sample is marked as SE and `sppIDer.py` is called without `--r2`. SE samples that do not match any PE suffix are also collected directly.

**Multi-run merging.** Samples whose filenames end with `_run1`, `_run2`, etc. are grouped under the same base name. Before analysis their read files are concatenated into a single temporary FASTQ pair, which is deleted after the run (unless `-k` is passed).

**SLURM array creation.** The script writes a `sppIDer_samples.tsv` file listing all detected samples with their mode and file paths, then generates and submits a `.slurm` worker script as an array job. Each array task processes a sequential batch of samples. If there are fewer samples than the requested number of jobs, the array size is reduced automatically.

**Logs.** Per-job stdout and stderr go to `{outdir}/logs/sppIDer_{jobID}_{taskID}.out/.err`. The script prints `squeue`, `scancel`, and `tail` commands at the end for convenience.

### Example with manual suffixes and run merging

```bash
# Reads directory contains:
#   strain42_run1_R1.trimmed.fastq.gz  strain42_run1_R2.trimmed.fastq.gz
#   strain42_run2_R1.trimmed.fastq.gz  strain42_run2_R2.trimmed.fastq.gz
#   strain99_R1.trimmed.fastq.gz       strain99_R2.trimmed.fastq.gz
#   strainOld.trimmed.fastq.gz         (no R2 → SE)

bash run_sppIDer_array_se.sh \
    -i /lustre/reads/trimmed/ \
    -r /lustre/ref/SaccharomycesCombo.fasta \
    -o /lustre/results/batch3 \
    -s _R1.trimmed.fastq.gz \
    -S _R2.trimmed.fastq.gz \
    -j 30
```

The script will report:
```
  [PE][2 runs] strain42
  [PE][1 run ] strain99
  [SE][1 run ] strainOld

  Total: 3 samples (2 PE, 1 SE)
```

---

## Changes relative to the original sppIDer

The original pipeline (GLBRC/sppIDer) was designed to run inside a Docker container with all inputs and outputs in a single working directory. This version adapts it to a shared HPC environment and extends it in several ways.

### Update
Python scripts were refactored from python2 to python3.

### Structural changes

**Decoupled input paths.** The original required all input files (reference FASTA, fastq reads) to be in the same directory as the working output. Here every input — reference, R1, R2 — takes its own independent path, so reads and references can live anywhere on the filesystem.

**`--sample` / `--outdir` model.** The original used a single `--out` prefix that doubled as both the output filename prefix and an implicit working directory. The rewrite splits this into `--sample` (filename prefix) and `--outdir` (output directory), making it easier to organise results per sample.

**`utils.py` shared module.** Logging, subprocess management, path resolution, and elapsed-time formatting were factored out into a shared module used by all pipeline scripts, replacing ad-hoc inline code in the original.

### Dependency changes

**`doBy` replaced by `splitBy_local.R`.** The original R scripts depended on `doBy::splitBy()`. That package pulls in a long chain of dependencies (`Deriv`, `broom`, `modelr`, `microbenchmark`) which caused version conflicts with `ggplot2` on the cluster R installation. A minimal local reimplementation of the one function actually used was written instead.

**`data.table` added for `-d` mode.** The per-basepair bedgraph files produced by `bedtools -d` can be tens of gigabytes for nuclear genomes. `meanDepth_sppIDer-d.R` now uses `data.table::fread()` for efficient columnar reading and processes one species at a time via an `awk | grep` pipe to keep memory usage bounded.

### SLURM launcher (entirely new)

The original repository provides no cluster submission scripts. `run_sppIDer_array_se.sh` was written for Drago: it distributes N samples across a SLURM array job with automatic suffix detection, SE/PE classification, and multi-run merging.

### Aggregated report (entirely new)

`aggregate_sppIDer_report.py` is run once after the pipeline has finished, against the top-level output directory:

```bash
python3 aggregate_sppIDer_report.py /path/to/sppider_out
```

It walks every per-sample subdirectory, reads `_MQsummary.txt` and `_speciesAvgDepth-[d|g].txt`, and produces two files in the output root:

- `species_assignment_summary.tsv` — one row per sample with read counts, per-species mapping percentages, per-species mean coverage, a species call, and any quality flags.
- `sppIDer_report.html` — an interactive report with summary cards, a call distribution table, and a filterable/sortable per-sample table with tab-based filtering, a coverage threshold slider, and TSV export.

Each sample receives a **call** based on configurable thresholds:

| Call | Condition |
|------|-----------|
| `pure_{species}` | Best species ≥ 70% of mapped reads |
| `hybrid_{sp1}_x_{sp2}` | Best ≥ 30% and second ≥ 20% |
| `ambiguous` | Neither condition met |

Quality **flags** are raised independently of the call:

| Flag | Condition |
|------|-----------|
| `LOW_READS` | Total reads < 500,000 |
| `HIGH_UNMAPPED` | Unmapped reads > 20% |
| `LOW_COV` | Mean depth of best species < 5× |
| `CONTAMINATION` | Any unexpected species ≥ 10% of mapped reads |

All thresholds are defined as constants at the top of the script (`MIN_READS`, `MAX_UNMAPPED_PCT`, `PURE_THRESHOLD`, `HYBRID_BEST_MIN`, `HYBRID_SEC_MIN`, `CONTAM_THRESHOLD`, `LOW_COV`) and can be adjusted before running.
