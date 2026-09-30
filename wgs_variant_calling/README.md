# WGS variant calling

Whole-genome sequencing pipeline used to identify suppressor mutations. It runs on an HPC cluster with SLURM (tested on HPC-Drago, CSIC): each launcher script checks its inputs, writes the SLURM job scripts and submits them with the right dependencies.

The 11 strains (wild-type BY4741, the parental *tif51A-1* and *tif51A-3* strains, and eight suppressors) are processed together: variants are called per strain in GVCF mode and then **genotyped jointly**, so that every strain has a genotype and its own read depth at every variant site. This is what allows suppressor-specific mutations to be defined by comparison with the sequenced control strains (see `../suppressor_variants`).

## Workflow

| Step | Script | What it does |
|---|---|---|
| 1 | `run_download_sra.sh` | Downloads the 11 runs of BioProject PRJNA1418127 from SRA (`fasterq-dump`) and checks R1/R2 pairing |
| 2 | `run_trim_galore.sh` | Trim Galore (Nextera adapters, first 15 bp clipped from R1 and R2, NextSeq quality trimming) |
| 3 | `run_bwa_gatk.sh` | Reference indexes; BWA-MEM2 alignment, sorting, GATK MarkDuplicates, flagstat |
| 4 | `run_variant_calling.sh` | **Initial call set**: HaplotypeCaller (`-ploidy 1 -ERC BP_RESOLUTION`) → CombineGVCFs → GenotypeGVCFs (`-ploidy 1`) → hard filters |
| 5 | `run_bqsr.sh` | Base quality score recalibration (BaseRecalibrator / ApplyBQSR) using the PASS variants of step 4 as known sites (bootstrapping); AnalyzeCovariates before/after |
| 6 | `run_variant_calling.sh` | **Final call set**: same as step 4 on the recalibrated BAMs |

The final VCF (`<output>/vcf/cohort.filtered.vcf.gz`) is annotated with `../snpeff/annotate_snpeff.sh`.

### Key choices
- **Haploid ploidy** (`-ploidy 1`) in HaplotypeCaller and GenotypeGVCFs.
- **Per-base reference confidence** (`-ERC BP_RESOLUTION`): with the default GVCF reference blocks, strains without a variant get the minimum depth of the whole block instead of the depth at the site.
- **Joint genotyping** of the 11 strains: a strain that cannot be confidently genotyped at a site is reported as no-call (`.`), never silently as reference.
- **Hard filters** (GATK VariantFiltration; failing sites are flagged, not removed):
  - SNPs: QD < 2.0, QUAL < 30, SOR > 3.0, FS > 60.0, MQ < 40.0, MQRankSum < −12.5, ReadPosRankSum < −8.0
  - Indels (and mixed sites): QD < 2.0, QUAL < 30, FS > 200.0, ReadPosRankSum < −20.0

## Reference genome
*S. cerevisiae* S288C, assembly R64-1-1 (SGD release 20110203), nuclear chromosomes I–XVI with headers `chrI`…`chrXVI`. The reference is only a coordinate system: the strains are BY4741-derived and differ from S288C at many positions, so suppressor mutations are defined relative to the sequenced control strains.

## Usage

Run every launcher from a login node with the `gatk` conda environment active; `-h` prints the full help.

```bash
conda activate gatk
REF=/path/to/ref/nuclear.fasta

bash run_download_sra.sh  -o data
bash run_trim_galore.sh   -i data/fastq -o trimmed
bash run_bwa_gatk.sh      -r $REF -i trimmed/trimmed -o mapped
bash run_variant_calling.sh -r $REF -i mapped/bam -o vc_round0
bash run_bqsr.sh          -r $REF -k vc_round0/vcf/cohort.filtered.PASS.vcf.gz -i mapped/bam -o bqsr
bash run_variant_calling.sh -r $REF -i bqsr/bam -o vc_final
```

Wait for each step to finish before launching the next one (steps inside a script are chained with SLURM dependencies; steps across scripts are not). All scripts skip samples whose outputs already exist, so they can be re-launched after a failure; `-s N` restarts a script from step N.

Common options: `-j` maximum number of queued jobs (default 60), `-t` threads per job.

## Outputs of `run_variant_calling.sh`
- `gvcf/SAMPLE.g.vcf.gz` — per-strain GVCF (per-base)
- `vcf/cohort.raw.vcf.gz` — joint-genotyped variants of the 11 strains
- `vcf/cohort.filtered.vcf.gz` — same, with hard-filter results in FILTER
- `vcf/cohort.filtered.PASS.vcf.gz` — PASS sites only (used as known sites for BQSR)
- `vcf/cohort.filtered.table.tsv` — GATK VariantsToTable export

## Dependencies
- SLURM
- Environment modules: `GATK/4.6.0.0-Java-17`, `SRA-Toolkit/3.2.0`, `Trim_Galore/0.6.10`
- Conda environment `gatk`: `bwa-mem2`, `samtools`
- Conda environment `download`: `pigz`
