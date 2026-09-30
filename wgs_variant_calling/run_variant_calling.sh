#!/bin/bash
# =============================================================================
# run_variant_calling.sh — Joint variant calling with GATK (HPC Drago, SLURM)
#
# Step 1: HaplotypeCaller -ploidy 1 -ERC BP_RESOLUTION per sample (array)
# Step 2: CombineGVCFs + GenotypeGVCFs -ploidy 1 (single job)
# Step 3: hard filters (SNP / non-SNP), MergeVcfs and VariantsToTable (single job)
#
# Output: gvcf/SAMPLE.g.vcf.gz, vcf/cohort.raw.vcf.gz, vcf/cohort.filtered.vcf.gz,
#         vcf/cohort.filtered.PASS.vcf.gz, vcf/cohort.filtered.table.tsv
#
# Joint genotyping gives every sample a genotype at every variant site:
# 0 = reference with coverage, ./. = no confident call (never assumed to be reference).
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
REF=""
INPUT_DIRS=()
OUTPUT_DIR=""
MAX_JOBS=60
THREADS=4
START_STEP=1

usage() {
    cat << USAGE
Usage: $(basename "$0") -r REF.fa -i DIR_BAM [-i DIR2 ...] -o OUTPUT_DIR [-j MAX_JOBS] [-t THREADS] [-s STEP] [-h]

Calls variants in deduplicated (run_bwa_gatk.sh) or recalibrated (run_bqsr.sh) BAMs with HaplotypeCaller
in per-base GVCF mode (BP_RESOLUTION, haploid), genotypes all samples jointly and applies the GATK-recommended
hard filters. Failing sites are flagged in FILTER, not removed.
All work is submitted to SLURM; this script only writes and submits the jobs.

Flags:
  -r  Reference FASTA (required). Must have .fai and .dict (created by run_bwa_gatk.sh)
  -i  Directory with *_recal.bam or *_sorted_dedup.bam + .bai (required, repeatable).
      If a directory has both, the *_recal.bam files are used
  -o  Output directory (required). GVCF in OUTPUT_DIR/gvcf, VCF in OUTPUT_DIR/vcf,
      logs in OUTPUT_DIR/logs
  -j  Maximum number of queued jobs (default: ${MAX_JOBS} = MaxJobsPU)
  -t  Threads per HaplotypeCaller job (default: ${THREADS})
  -s  Step to start from (default: 1)
        1 = HaplotypeCaller per sample
        2 = CombineGVCFs + GenotypeGVCFs
        3 = filters + table
  -h  Show this help

Examples:
  conda activate gatk
  $(basename "$0") -r /lustre/home/iata/aaguilar/paper/ref/nuclear.fasta \\
      -i /lustre/home/iata/aaguilar/paper/mapping/bam -o /lustre/home/iata/aaguilar/paper/calling
  $(basename "$0") -r ref.fa -i bam -o calling -s 3      # redo the filters only

Samples whose GVCF and index already exist are skipped when re-launching.
USAGE
}

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
while getopts "r:i:o:j:t:s:h" opt; do
    case ${opt} in
        r) REF=${OPTARG} ;;
        i) INPUT_DIRS+=("${OPTARG}") ;;
        o) OUTPUT_DIR=${OPTARG} ;;
        j) MAX_JOBS=${OPTARG} ;;
        t) THREADS=${OPTARG} ;;
        s) START_STEP=${OPTARG} ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

if [[ -z "${REF}" || -z "${OUTPUT_DIR}" || ${#INPUT_DIRS[@]} -eq 0 ]]; then
    echo "ERROR: missing required arguments (-r, -i, -o)."
    usage
    exit 1
fi
if ! [[ "${MAX_JOBS}" =~ ^[0-9]+$ && "${THREADS}" =~ ^[0-9]+$ && "${START_STEP}" =~ ^[123]$ ]]; then
    echo "ERROR: -j and -t must be integers and -s must be 1, 2 or 3."
    exit 1
fi
[[ -f "${REF}" ]] || { echo "ERROR: reference ${REF} not found"; exit 1; }
for DIR in "${INPUT_DIRS[@]}"; do
    [[ -d "${DIR}" ]] || { echo "ERROR: directory ${DIR} not found"; exit 1; }
done

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
if [[ "${CONDA_DEFAULT_ENV:-}" != "gatk" ]]; then
    echo "WARNING: activate the conda environment 'gatk' before launching."
    echo "  conda activate gatk"
    exit 1
fi

# ---------------------------------------------------------------------------
# Directories
# ---------------------------------------------------------------------------
REF=$(realpath "${REF}")
OUTPUT_DIR=$(realpath -m "${OUTPUT_DIR}")
GVCF_DIR="${OUTPUT_DIR}/gvcf"
VCF_DIR="${OUTPUT_DIR}/vcf"
mkdir -p "${GVCF_DIR}" "${VCF_DIR}" "${OUTPUT_DIR}/logs" "${OUTPUT_DIR}/tmp"

for f in "${REF}.fai" "${REF%.*}.dict"; do
    [[ -f "${f}" ]] || { echo "ERROR: ${f} not found. Create it with step 1 of run_bwa_gatk.sh."; exit 1; }
done

# ---------------------------------------------------------------------------
# Sample list: SAMPLE<TAB>BAM
# ---------------------------------------------------------------------------
LIST_FILE="${OUTPUT_DIR}/samples.tsv"
: > "${LIST_FILE}"
for DIR in "${INPUT_DIRS[@]}"; do
    DIR=$(realpath "${DIR}")
    BAMS=$(find "${DIR}" -maxdepth 1 -type f -name "*_recal.bam" | sort)
    [[ -z "${BAMS}" ]] && BAMS=$(find "${DIR}" -maxdepth 1 -type f -name "*_sorted_dedup.bam" | sort)
    while read -r BAM; do
        [[ -z "${BAM}" ]] && continue
        SAMPLE=$(basename "${BAM}" .bam)
        SAMPLE=${SAMPLE%_recal}
        SAMPLE=${SAMPLE%_sorted_dedup}
        if [[ ! -f "${BAM}.bai" ]]; then
            echo "WARNING: no index for ${BAM}, skipped."
            continue
        fi
        if cut -f1 "${LIST_FILE}" | grep -qxF "${SAMPLE}"; then
            echo "ERROR: sample duplicated across directories: ${SAMPLE}"
            exit 1
        fi
        printf '%s\t%s\n' "${SAMPLE}" "${BAM}" >> "${LIST_FILE}"
    done <<< "${BAMS}"
done

TOTAL=$(wc -l < "${LIST_FILE}")
[[ ${TOTAL} -gt 0 ]] || { echo "ERROR: no *_recal.bam or *_sorted_dedup.bam found in ${INPUT_DIRS[*]}"; exit 1; }

if [[ "${START_STEP}" -eq 2 ]]; then
    while IFS=$'\t' read -r SAMPLE BAM; do
        [[ -s "${GVCF_DIR}/${SAMPLE}.g.vcf.gz.tbi" ]] \
            || { echo "ERROR: ${GVCF_DIR}/${SAMPLE}.g.vcf.gz(.tbi) not found. Start from step 1."; exit 1; }
    done < "${LIST_FILE}"
fi
if [[ "${START_STEP}" -eq 3 && ! -s "${VCF_DIR}/cohort.raw.vcf.gz.tbi" ]]; then
    echo "ERROR: ${VCF_DIR}/cohort.raw.vcf.gz(.tbi) not found. Start from step 2."
    exit 1
fi

# Steps 2 and 3 are single jobs: they are subtracted from the array pool
NUM_SEQ_JOBS=$(( 4 - START_STEP ))
[[ ${NUM_SEQ_JOBS} -gt 2 ]] && NUM_SEQ_JOBS=2
NUM_JOBS=$(( MAX_JOBS - NUM_SEQ_JOBS ))
[[ ${NUM_JOBS} -gt ${TOTAL} ]] && NUM_JOBS=${TOTAL}
[[ ${NUM_JOBS} -lt 1 ]] && { echo "ERROR: -j demasiado bajo."; exit 1; }

# ---------------------------------------------------------------------------
# SLURM script — Step 1: HaplotypeCaller per sample (array)
# ---------------------------------------------------------------------------
SCRIPT1="${OUTPUT_DIR}/haplotypecaller.slurm"
cat > "${SCRIPT1}" << 'EOF'
#!/bin/bash
#SBATCH --job-name=hc_gvcf
#SBATCH --partition=generic
#SBATCH --array=0-__LAST_TASK__
#SBATCH --cpus-per-task=__THREADS__
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=__OUTPUT_DIR__/logs/hc_gvcf_%A_%a.out
#SBATCH --error=__OUTPUT_DIR__/logs/hc_gvcf_%A_%a.err

set -eo pipefail
module purge
module load rama0.4 GCCcore/13.2.0 GATK/4.6.0.0-Java-17
eval "$(conda shell.bash hook)"
conda activate gatk
set -u

REF="__REF__"
LIST_FILE="__LIST_FILE__"
GVCF_DIR="__GVCF_DIR__"
TMP_DIR="__OUTPUT_DIR__/tmp"
THREADS=__THREADS__
NUM_JOBS=__NUM_JOBS__
TOTAL=$(wc -l < "${LIST_FILE}")

echo "=== HaplotypeCaller (per-base GVCF, ploidy 1) ==="
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID}, Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Reference: ${REF}"
echo "Start: $(date)"
echo "========================"

ITEMS_PER_JOB=$(( (TOTAL + NUM_JOBS - 1) / NUM_JOBS ))
START_LINE=$(( SLURM_ARRAY_TASK_ID * ITEMS_PER_JOB + 1 ))
END_LINE=$(( START_LINE + ITEMS_PER_JOB - 1 ))
[[ ${END_LINE} -gt ${TOTAL} ]] && END_LINE=${TOTAL}
[[ ${START_LINE} -gt ${TOTAL} ]] && { echo "No items assigned"; exit 0; }

OK=0
FAILED=0
N=0
N_BATCH=$(( END_LINE - START_LINE + 1 ))
while IFS=$'\t' read -r SAMPLE BAM; do
    N=$(( N + 1 ))
    GVCF="${GVCF_DIR}/${SAMPLE}.g.vcf.gz"
    SAMPLE_TMP="${TMP_DIR}/${SAMPLE}"
    echo
    echo "[${N}/${N_BATCH}] Processing: ${SAMPLE}"
    echo "  BAM: ${BAM}"
    echo "  Time: $(date +%H:%M:%S)"

    if [[ -s "${GVCF}" && -s "${GVCF}.tbi" ]]; then
        echo "  Already called, skipped."
        OK=$(( OK + 1 ))
        continue
    fi
    mkdir -p "${SAMPLE_TMP}"

    # BP_RESOLUTION: one line per base, no reference blocks, so the DP of samples
    # without the variant is the depth at that position, not the minimum of a block
    # set -e does not act inside an 'if': steps are chained with && to stop at the first failure
    if echo "  [1/1] HaplotypeCaller..." \
        && gatk --java-options "-Xmx12g" HaplotypeCaller \
            -R "${REF}" \
            -I "${BAM}" \
            -O "${GVCF}" \
            -ploidy 1 \
            -ERC BP_RESOLUTION \
            --native-pair-hmm-threads "${THREADS}" \
            --tmp-dir "${SAMPLE_TMP}" \
        && [[ -s "${GVCF}.tbi" ]] \
        && rm -rf "${SAMPLE_TMP:?}"; then
        echo "  OK: ${SAMPLE}"
        OK=$(( OK + 1 ))
    else
        echo "  ERROR: ${SAMPLE}" >&2
        rm -f "${GVCF}" "${GVCF}.tbi"
        rm -rf "${SAMPLE_TMP:?}"
        FAILED=$(( FAILED + 1 ))
    fi
done < <(sed -n "${START_LINE},${END_LINE}p" "${LIST_FILE}")

echo
echo "========================"
echo "Task ${SLURM_ARRAY_TASK_ID} summary:"
echo "  Processed: ${OK}"
echo "  Failed:    ${FAILED}"
echo "  End: $(date)"
echo "========================"
[[ ${FAILED} -eq 0 ]]
EOF

# ---------------------------------------------------------------------------
# SLURM script — Step 2: CombineGVCFs + GenotypeGVCFs (single job)
# ---------------------------------------------------------------------------
SCRIPT2="${OUTPUT_DIR}/genotype_gvcfs.slurm"
cat > "${SCRIPT2}" << 'EOF'
#!/bin/bash
#SBATCH --job-name=genotype
#SBATCH --partition=generic
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=__OUTPUT_DIR__/logs/genotype_%j.out
#SBATCH --error=__OUTPUT_DIR__/logs/genotype_%j.err

set -eo pipefail
module purge
module load rama0.4 GCCcore/13.2.0 GATK/4.6.0.0-Java-17
eval "$(conda shell.bash hook)"
conda activate gatk
set -u

REF="__REF__"
LIST_FILE="__LIST_FILE__"
GVCF_DIR="__GVCF_DIR__"
VCF_DIR="__VCF_DIR__"
TMP_DIR="__OUTPUT_DIR__/tmp"
COMBINED="${VCF_DIR}/cohort.g.vcf.gz"
RAW="${VCF_DIR}/cohort.raw.vcf.gz"

echo "=== CombineGVCFs + GenotypeGVCFs ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Reference: ${REF}"
echo "Samples: $(wc -l < "${LIST_FILE}")"
echo "Start: $(date)"
echo "========================"

VARIANT_ARGS=()
while IFS=$'\t' read -r SAMPLE BAM; do
    GVCF="${GVCF_DIR}/${SAMPLE}.g.vcf.gz"
    [[ -s "${GVCF}.tbi" ]] || { echo "ERROR: ${GVCF}(.tbi) not found" >&2; exit 1; }
    VARIANT_ARGS+=(-V "${GVCF}")
done < "${LIST_FILE}"

echo "  [1/2] CombineGVCFs..."
gatk --java-options "-Xmx12g" CombineGVCFs \
    -R "${REF}" \
    "${VARIANT_ARGS[@]}" \
    -O "${COMBINED}" \
    --tmp-dir "${TMP_DIR}"

echo "  [2/2] GenotypeGVCFs (ploidy 1)..."
gatk --java-options "-Xmx12g" GenotypeGVCFs \
    -R "${REF}" \
    -V "${COMBINED}" \
    -O "${RAW}" \
    -ploidy 1 \
    --tmp-dir "${TMP_DIR}"

echo "  Raw variants: $(zgrep -vc '^#' "${RAW}")"

echo
echo "========================"
echo "  End: $(date)"
echo "========================"
EOF

# ---------------------------------------------------------------------------
# SLURM script — Step 3: hard filters + table (single job)
# ---------------------------------------------------------------------------
SCRIPT3="${OUTPUT_DIR}/filter_variants.slurm"
cat > "${SCRIPT3}" << 'EOF'
#!/bin/bash
#SBATCH --job-name=filter_vcf
#SBATCH --partition=generic
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=__OUTPUT_DIR__/logs/filter_vcf_%j.out
#SBATCH --error=__OUTPUT_DIR__/logs/filter_vcf_%j.err

set -eo pipefail
module purge
module load rama0.4 GCCcore/13.2.0 GATK/4.6.0.0-Java-17
eval "$(conda shell.bash hook)"
conda activate gatk
set -u

REF="__REF__"
VCF_DIR="__VCF_DIR__"
RAW="${VCF_DIR}/cohort.raw.vcf.gz"
P="${VCF_DIR}/cohort"

echo "=== GATK hard filters ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "VCF: ${RAW}"
echo "Start: $(date)"
echo "========================"

# Split into SNP / rest (INDEL + MIXED) so that mixed sites are not lost
echo "  [1/5] SelectVariants SNP / no-SNP..."
gatk SelectVariants -R "${REF}" -V "${RAW}" --select-type-to-include SNP -O "${P}.snps.raw.vcf.gz"
gatk SelectVariants -R "${REF}" -V "${RAW}" --select-type-to-exclude SNP -O "${P}.indels.raw.vcf.gz"

# GATK-recommended hard-filter thresholds. LowDP flags the genotype (FT), not the site
echo "  [2/5] VariantFiltration SNP..."
gatk VariantFiltration -R "${REF}" -V "${P}.snps.raw.vcf.gz" -O "${P}.snps.filtered.vcf.gz" \
    -filter "QD < 2.0"              --filter-name "QD2" \
    -filter "QUAL < 30.0"           --filter-name "QUAL30" \
    -filter "SOR > 3.0"             --filter-name "SOR3" \
    -filter "FS > 60.0"             --filter-name "FS60" \
    -filter "MQ < 40.0"             --filter-name "MQ40" \
    -filter "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -G-filter "DP < 5"              -G-filter-name "LowDP"

echo "  [3/5] VariantFiltration INDEL/MIXED..."
gatk VariantFiltration -R "${REF}" -V "${P}.indels.raw.vcf.gz" -O "${P}.indels.filtered.vcf.gz" \
    -filter "QD < 2.0"               --filter-name "QD2" \
    -filter "QUAL < 30.0"            --filter-name "QUAL30" \
    -filter "FS > 200.0"             --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -G-filter "DP < 5"               -G-filter-name "LowDP"

echo "  [4/5] MergeVcfs + PASS..."
gatk MergeVcfs -I "${P}.snps.filtered.vcf.gz" -I "${P}.indels.filtered.vcf.gz" -O "${P}.filtered.vcf.gz"
gatk SelectVariants -R "${REF}" -V "${P}.filtered.vcf.gz" --exclude-filtered -O "${P}.filtered.PASS.vcf.gz"

echo "  [5/5] VariantsToTable..."
gatk VariantsToTable -V "${P}.filtered.vcf.gz" -O "${P}.filtered.table.tsv" \
    -F CHROM -F POS -F REF -F ALT -F TYPE -F QUAL -F FILTER \
    -GF GT -GF AD -GF DP -GF GQ -GF FT \
    --show-filtered

rm -f "${P}".snps.raw.vcf.gz* "${P}".indels.raw.vcf.gz* "${P}".snps.filtered.vcf.gz* "${P}".indels.filtered.vcf.gz*

echo
echo "  Total sites: $(zgrep -vc '^#' "${P}.filtered.vcf.gz")"
echo "  PASS sites:  $(zgrep -vc '^#' "${P}.filtered.PASS.vcf.gz")"
echo "========================"
echo "  End: $(date)"
echo "========================"
EOF

for S in "${SCRIPT1}" "${SCRIPT2}" "${SCRIPT3}"; do
    sed -i \
        -e "s|__LAST_TASK__|$(( NUM_JOBS - 1 ))|g" \
        -e "s|__THREADS__|${THREADS}|g" \
        -e "s|__OUTPUT_DIR__|${OUTPUT_DIR}|g" \
        -e "s|__REF__|${REF}|g" \
        -e "s|__LIST_FILE__|${LIST_FILE}|g" \
        -e "s|__GVCF_DIR__|${GVCF_DIR}|g" \
        -e "s|__VCF_DIR__|${VCF_DIR}|g" \
        -e "s|__NUM_JOBS__|${NUM_JOBS}|g" \
        "${S}"
done

# ---------------------------------------------------------------------------
# Submission
# ---------------------------------------------------------------------------
PREV_DEP=""
JOB_LINES=()
N_SUBMITTED=0
if [[ "${START_STEP}" -le 1 ]]; then
    OUTPUT=$(sbatch "${SCRIPT1}")
    JOB1_ID=$(echo "${OUTPUT}" | awk '{print $NF}')
    PREV_DEP="--dependency=afterok:${JOB1_ID}"
    JOB_LINES+=(" Job 1 (${JOB1_ID}): HaplotypeCaller GVCF, array 0-$(( NUM_JOBS - 1 ))")
    N_SUBMITTED=$(( N_SUBMITTED + NUM_JOBS ))
fi
if [[ "${START_STEP}" -le 2 ]]; then
    OUTPUT=$(sbatch ${PREV_DEP} "${SCRIPT2}")
    JOB2_ID=$(echo "${OUTPUT}" | awk '{print $NF}')
    JOB_LINES+=(" Job 2 (${JOB2_ID}): CombineGVCFs + GenotypeGVCFs${PREV_DEP:+ (${PREV_DEP#--dependency=})}")
    PREV_DEP="--dependency=afterok:${JOB2_ID}"
    N_SUBMITTED=$(( N_SUBMITTED + 1 ))
fi
OUTPUT=$(sbatch ${PREV_DEP} "${SCRIPT3}")
JOB3_ID=$(echo "${OUTPUT}" | awk '{print $NF}')
JOB_LINES+=(" Job 3 (${JOB3_ID}): filters + table${PREV_DEP:+ (${PREV_DEP#--dependency=})}")
N_SUBMITTED=$(( N_SUBMITTED + 1 ))

echo "============================================="
echo " Summary"
echo "============================================="
echo " Reference:    ${REF}"
for DIR in "${INPUT_DIRS[@]}"; do
echo " Input:        $(realpath "${DIR}")"
done
echo " Samples:      ${TOTAL} (${LIST_FILE})"
echo " Output GVCF:  ${GVCF_DIR}"
echo " Output VCF:   ${VCF_DIR}"
echo " Threads/job:  ${THREADS}"
echo " Start step:   ${START_STEP}"
echo "---------------------------------------------"
printf '%s\n' "${JOB_LINES[@]}"
echo "---------------------------------------------"
echo " Total jobs: ${N_SUBMITTED} / ${MAX_JOBS}"
echo "============================================="
echo
echo "Useful commands:"
echo "  squeue -u \$(whoami)"
echo "  tail -F ${OUTPUT_DIR}/logs/hc_gvcf_*.out"
echo "  grep -h 'OK:\|ERROR' ${OUTPUT_DIR}/logs/hc_gvcf_*.out"
echo "  grep -l 'ERROR' ${OUTPUT_DIR}/logs/*.err"
echo "  tail ${OUTPUT_DIR}/logs/filter_vcf_*.out"
