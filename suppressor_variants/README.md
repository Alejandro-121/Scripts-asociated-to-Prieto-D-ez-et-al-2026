# Suppressor-specific variants

Scripts that turn the annotated joint VCF into the supplementary tables and Figure 3a. They use only the VCF (no BAM files).

## `find_suppressor_variants.py`
Identifies the variants that appeared in the suppressors and writes Supplementary Table S6 (`.xlsx` with sheets *Variants*, *Summary* and *Criteria*, plus a `.tsv`).

**Candidates:** alternative alleles genotyped in at least one suppressor and in none of the three control strains (wild-type BY4741, *tif51A-1*, *tif51A-3*). Each candidate is classified by what the controls show at that site:

| Class | Criterion | Meaning |
|---|---|---|
| background | a control has > 2 reads supporting the allele, or carries another alternative allele at the site | present in the controls (strain background, repetitive or duplicated regions) |
| not testable | a control is a no-call (`.`, GQ 0) or has < 10 reads | the controls cannot be assessed |
| suppressor-specific | all three controls genotyped as reference (`0`) with ≥ 10 reads and ≤ 2 reads supporting the allele | absent from the controls |

A no-call in a control is never treated as reference.

Suppressor-specific variants are listed in the table only when the **evidence in the suppressor** is sufficient: the site passes the GATK hard filters, the allele is present in ≥ 80 % of the reads of every carrier (haploid, clonal strains) and at least one carrier has ≥ 10 reads. The number of candidates excluded at each step is given in the *Summary* sheet.

Annotation (gene, effect, impact, HGVS c. and p.) comes from SnpEff; gene names are updated to the current SGD standard names using the current SGD GFF. The *Experimental validation* column records the alleles reconstructed by CRISPR-Cas9 in the paper (ROX1 W53L, L73P, W64C; MOT3 INS1).

```bash
python find_suppressor_variants.py \
    -v vc_final/vcf/cohort.filtered.ann.vcf.gz \
    --gff saccharomyces_cerevisiae.gff.gz \
    -o Supplementary_Table_WGS_variants
```

Thresholds can be changed with `--min-ctrl-dp` (10), `--max-ctrl-alt` (2), `--min-af` (0.8) and `--min-carrier-dp` (10). Sample names are the ones in the VCF read groups (`WT`, `2-1` = *tif51A-1*, `2-3` = *tif51A-3*, `sup.N`).

## `make_fig3a.py`
Draws Figure 3a (`Figure3a_corrected.pdf/.svg/.png`) from the table: protein-altering suppressor-specific mutations, SNPs and INDELs in separate panels, suppressors grouped by parental strain. Run it in the directory that contains `Supplementary_Table_WGS_variants.xlsx`.

## `make_accession_table.py`
Writes Supplementary Table S7 with the BioProject, BioSample and SRA run accessions of the 11 sequenced strains.

## Dependencies
Python 3.10+, `openpyxl`, `matplotlib`.
