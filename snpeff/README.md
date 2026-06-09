 # Helper Scripts

A small collection of utility scripts for preparing and annotating yeast (*S. cerevisiae*) genomic data.

---

## Scripts

### `rename_headers.py`

Renames the FASTA headers of the S288C reference genome to short chromosome names (`chrI`, `chrII`, … `chrMito`).

**Usage:**
```bash
python rename_headers.py
```

Edit the `os.chdir(...)` path and `in_fasta` variable at the top of the script to point to your local copy of the reference genome before running.

**Input:** `S288C_reference_sequence_R64-1-1_20110203.fsa`  
**Output:** `new_head.fasta`

---

### `run_snpeff.sh`

Batch-annotates VCF files in the current directory using [SnpEff](https://pcingola.github.io/SnpEff/) against the `R64-1-1_sgd` database.

For each `.vcf` file it:
1. Compresses it with `bgzip` (in parallel)
2. Indexes it with `tabix`
3. Runs SnpEff annotation
4. Organises outputs (annotated VCF, stats, CSV) into a per-sample directory
5. Copies the annotated VCF to `out_vcf/`

**Usage:**
```bash
# Place run_snpeff.sh and snpEff.jar in the same directory as your .vcf files, then:
bash run_snpeff.sh
```

**Dependencies:** `bgzip`, `tabix`, `parallel`, `java`

---

## Installing SnpEff

Download the latest version from the official site:

> https://pcingola.github.io/SnpEff/#download

Direct download (latest stable):
```bash
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
```

Then download the yeast database:
```bash
java -jar snpEff.jar download R64-1-1_sgd
```

Place `snpEff.jar` in the same directory as `run_snpeff.sh`, or adjust the path inside the script.

---
