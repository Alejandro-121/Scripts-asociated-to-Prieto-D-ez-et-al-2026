# SnpEff annotation

`annotate_snpeff.sh` builds a custom SnpEff database for *S. cerevisiae* S288C **R64-1-1** and annotates the joint VCF produced by `../wgs_variant_calling`.

## Why a custom database
The database is built from exactly the files used for the alignment, so that coordinates, sequence and annotation all come from the same release:

- **Genome:** the reference FASTA used for mapping (nuclear chromosomes `chrI`…`chrXVI`, identical to SGD `S288C_reference_sequence_R64-1-1_20110203.fsa`).
- **Annotation:** SGD GFF3 `saccharomyces_cerevisiae_R64-1-1_20110208.gff`, restricted to the nuclear chromosomes and without its embedded `##FASTA` section. The GFF uses the systematic name in `Name` and the standard name in `gene=`; the script sets `Name` to the standard name so that SnpEff reports e.g. `ROX1` instead of `YPR065W`.
- **Check:** coding and protein sequences are validated against SGD `orf_coding_all` and `orf_trans_all` (R64-1-1: 6,575 transcripts, 0 errors). The script stops if the check reports errors.

Because the reference FASTA already uses the GFF chromosome names, no header renaming is needed.

## Usage
Download and unpack the SGD R64-1-1 bundle (`S288C_reference_genome_R64-1-1_20110203.tgz`) and run:

```bash
conda create -n snpeff -c conda-forge -c bioconda snpeff bcftools htslib
conda activate snpeff

bash annotate_snpeff.sh \
    -r ref/nuclear.fasta \
    -g S288C_reference_genome_R64-1-1_20110203 \
    -v vc_final/vcf/cohort.filtered.vcf.gz \
    -d snpeff_db
```

- `-r` reference FASTA used for mapping
- `-g` SGD R64-1-1 directory (GFF, `orf_coding_all`, `orf_trans_all`)
- `-v` VCF to annotate (bgzip)
- `-d` directory for `snpEff.config` and the database
- `-s 2` annotate only, reusing an existing database

## Outputs (next to the input VCF)
- `*.ann.vcf.gz` (+ `.tbi`) — annotated VCF (`ANN`, `LOF`, `NMD` fields)
- `snpeff_summary.html`, `snpeff_stats.csv`, `snpeff.log`

Tested with SnpEff 5.4c.
