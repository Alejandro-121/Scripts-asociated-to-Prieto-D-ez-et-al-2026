#!/bin/bash
# =============================================================================
# run_bqsr.sh — Base quality score recalibration (BQSR) with GATK (HPC Drago, SLURM)
#
# Single step (array), per sample:
#   BaseRecalibrator (before) → ApplyBQSR → BaseRecalibrator (after) → AnalyzeCovariates
#
# There is no catalogue of known variants for the strain: the PASS variants of an initial
# call set are used as known sites (bootstrapping; cohort.filtered.PASS.vcf.gz from
# run_variant_calling.sh on the non-recalibrated BAMs).
#
# Output per sample: SAMPLE_recal.bam(.bai), SAMPLE_bqsr_before.table,
#                    SAMPLE_bqsr_after.table, SAMPLE_bqsr_covariates.csv (and .pdf if R is available)
# Then: run_variant_calling.sh -i OUTPUT_DIR/bam for the final call set.
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
REF=""
KNOWN=""
INPUT_DIRS=()
OUTPUT_DIR=""
MAX_JOBS=60
THREADS=2

usage() {
    cat << USAGE
Usage: $(basename "$0") -r REF.fa -k KNOWN.vcf.gz -i DIR_BAM [-i DIR2 ...] -o OUTPUT_DIR [-j MAX_JOBS] [-t THREADS] [-h]

Recalibrates base qualities of deduplicated BAMs (output of run_bwa_gatk.sh) with
GATK BaseRecalibrator + ApplyBQSR, using the high-confidence variants of an initial
call set as known sites. The recalibration is evaluated with a second BaseRecalibrator
pass and AnalyzeCovariates.
All work is submitted to SLURM; this script only writes and submits the jobs.

Flags:
  -r  Reference FASTA (required). Must have .fai and .dict (created by run_bwa_gatk.sh)
  -k  VCF of known sites, bgzip + .tbi (required). Usually
        <output of run_variant_calling.sh>/vcf/cohort.filtered.PASS.vcf.gz
  -i  Directory with *_sorted_dedup.bam + .bai (required, repeatable)
  -o  Output directory (required). BAM in OUTPUT_DIR/bam, tables in OUTPUT_DIR/tables,
      logs in OUTPUT_DIR/logs
  -j  Maximum number of queued jobs (default: ${MAX_JOBS} = MaxJobsPU)
  -t  Threads per job (default: ${THREADS})
  -h  Show this help

Examples:
  conda activate gatk
  $(basename "$0") -r /lustre/home/iata/aaguilar/paper/ref/nuclear.fasta \\
      -k /lustre/home/iata/aaguilar/paper/vc_bp/vcf/cohort.filtered.PASS.vcf.gz \\
      -i /lustre/home/iata/aaguilar/paper/mapped/bam -o /lustre/home/iata/aaguilar/paper/bqsr

Samples whose recalibrated BAM and index already exist are skipped when re-launching.
USAGE
}

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
while getopts "r:k:i:o:j:t:h" opt; do
    case ${opt} in
        r) REF=${OPTARG} ;;
        k) KNOWN=${OPTARG} ;;
        i) INPUT_DIRS+=("${OPTARG}") ;;
        o) OUTPUT_DIR=${OPTARG} ;;
        j) MAX_JOBS=${OPTARG} ;;
        t) THREADS=${OPTARG} ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

if [[ -z "${REF}" || -z "${KNOWN}" || -z "${OUTPUT_DIR}" || ${#INPUT_DIRS[@]} -eq 0 ]]; then
    echo "ERROR: missing required arguments (-r, -k, -i, -o)."
    usage
    exit 1
fi
if ! [[ "${MAX_JOBS}" =~ ^[0-9]+$ && "${THREADS}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: -j and -t must be integers."
    exit 1
fi
[[ -f "${REF}" ]] || { echo "ERROR: reference ${REF} not found"; exit 1; }
[[ -f "${KNOWN}" ]] || { echo "ERROR: known-sites VCF ${KNOWN} not found"; exit 1; }
[[ -f "${KNOWN}.tbi" ]] || { echo "ERROR: index ${KNOWN}.tbi not found"; exit 1; }
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
command -v samtools > /dev/null || { echo "ERROR: 'samtools' is not in the 'gatk' environment."; exit 1; }

# ---------------------------------------------------------------------------
# Directories
# ---------------------------------------------------------------------------
REF=$(realpath "${REF}")
KNOWN=$(realpath "${KNOWN}")
OUTPUT_DIR=$(realpath -m "${OUTPUT_DIR}")
BAM_DIR="${OUTPUT_DIR}/bam"
TABLE_DIR="${OUTPUT_DIR}/tables"
mkdir -p "${BAM_DIR}" "${TABLE_DIR}" "${OUTPUT_DIR}/logs" "${OUTPUT_DIR}/tmp"

for f in "${REF}.fai" "${REF%.*}.dict"; do
    [[ -f "${f}" ]] || { echo "ERROR: ${f} not found. Create it with step 1 of run_bwa_gatk.sh."; exit 1; }
done
N_KNOWN=$(zgrep -vc '^#' "${KNOWN}" || true)
[[ ${N_KNOWN} -gt 0 ]] || { echo "ERROR: ${KNOWN} contains no variants."; exit 1; }

# ---------------------------------------------------------------------------
# Sample list: SAMPLE<TAB>BAM
# ---------------------------------------------------------------------------
LIST_FILE="${OUTPUT_DIR}/samples.tsv"
: > "${LIST_FILE}"
for DIR in "${INPUT_DIRS[@]}"; do
    DIR=$(realpath "${DIR}")
    while read -r BAM; do
        [[ -z "${BAM}" ]] && continue
        SAMPLE=$(basename "${BAM}" _sorted_dedup.bam)
        if [[ ! -f "${BAM}.bai" ]]; then
            echo "WARNING: no index for ${BAM}, skipped."
            continue
        fi
        if cut -f1 "${LIST_FILE}" | grep -qxF "${SAMPLE}"; then
            echo "ERROR: sample duplicated across directories: ${SAMPLE}"
            exit 1
        fi
        printf '%s\t%s\n' "${SAMPLE}" "${BAM}" >> "${LIST_FILE}"
    done <<< "$(find "${DIR}" -maxdepth 1 -type f -name "*_sorted_dedup.bam" | sort)"
done

TOTAL=$(wc -l < "${LIST_FILE}")
[[ ${TOTAL} -gt 0 ]] || { echo "ERROR: no *_sorted_dedup.bam found in ${INPUT_DIRS[*]}"; exit 1; }

NUM_JOBS=${MAX_JOBS}
[[ ${NUM_JOBS} -gt ${TOTAL} ]] && NUM_JOBS=${TOTAL}
[[ ${NUM_JOBS} -lt 1 ]] && { echo "ERROR: -j demasiado bajo."; exit 1; }

# ---------------------------------------------------------------------------
# SLURM script — BQSR per sample (array)
# ---------------------------------------------------------------------------
SCRIPT1="${OUTPUT_DIR}/bqsr.slurm"
cat > "${SCRIPT1}" << 'EOF'
#!/bin/bash
#SBATCH --job-name=bqsr
#SBATCH --partition=generic
#SBATCH --array=0-__LAST_TASK__
#SBATCH --cpus-per-task=__THREADS__
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --output=__OUTPUT_DIR__/logs/bqsr_%A_%a.out
#SBATCH --error=__OUTPUT_DIR__/logs/bqsr_%A_%a.err

set -eo pipefail
module purge
module load rama0.4 GCCcore/13.2.0 GATK/4.6.0.0-Java-17
eval "$(conda shell.bash hook)"
conda activate gatk
set -u

REF="__REF__"
KNOWN="__KNOWN__"
LIST_FILE="__LIST_FILE__"
BAM_DIR="__BAM_DIR__"
TABLE_DIR="__TABLE_DIR__"
TMP_DIR="__OUTPUT_DIR__/tmp"
THREADS=__THREADS__
NUM_JOBS=__NUM_JOBS__
TOTAL=$(wc -l < "${LIST_FILE}")

echo "=== BQSR (known sites: initial call set) ==="
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID}, Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Reference: ${REF}"
echo "Known sites: ${KNOWN}"
echo "Start: $(date)"
echo "========================"

ITEMS_PER_JOB=$(( (TOTAL + NUM_JOBS - 1) / NUM_JOBS ))
START_LINE=$(( SLURM_ARRAY_TASK_ID * ITEMS_PER_JOB + 1 ))
END_LINE=$(( START_LINE + ITEMS_PER_JOB - 1 ))
[[ ${END_LINE} -gt ${TOTAL} ]] && END_LINE=${TOTAL}
[[ ${START_LINE} -gt ${TOTAL} ]] && { echo "No items assigned"; exit 0; }

# AnalyzeCovariates only writes the PDF if R with ggplot2/gplots/gsalib is available; the CSV always
PLOTS_OK=0
Rscript -e 'for (p in c("ggplot2","gplots","gsalib")) stopifnot(requireNamespace(p, quietly=TRUE))' \
    > /dev/null 2>&1 && PLOTS_OK=1
[[ ${PLOTS_OK} -eq 1 ]] || echo "WARNING: no R/ggplot2/gsalib, AnalyzeCovariates will only write the CSV."

OK=0
FAILED=0
N=0
N_BATCH=$(( END_LINE - START_LINE + 1 ))
while IFS=$'\t' read -r SAMPLE BAM; do
    N=$(( N + 1 ))
    RECAL="${BAM_DIR}/${SAMPLE}_recal.bam"
    BEFORE="${TABLE_DIR}/${SAMPLE}_bqsr_before.table"
    AFTER="${TABLE_DIR}/${SAMPLE}_bqsr_after.table"
    COVAR="${TABLE_DIR}/${SAMPLE}_bqsr_covariates"
    SAMPLE_TMP="${TMP_DIR}/${SAMPLE}"
    echo
    echo "[${N}/${N_BATCH}] Processing: ${SAMPLE}"
    echo "  BAM: ${BAM}"
    echo "  Time: $(date +%H:%M:%S)"

    if [[ -s "${RECAL}" && -s "${RECAL}.bai" && -s "${AFTER}" ]]; then
        echo "  Already recalibrated, skipped."
        OK=$(( OK + 1 ))
        continue
    fi
    mkdir -p "${SAMPLE_TMP}"
    PLOT_ARGS=()
    [[ ${PLOTS_OK} -eq 1 ]] && PLOT_ARGS=(-plots "${COVAR}.pdf")

    # set -e does not act inside an 'if': steps are chained with && to stop at the first failure
    if echo "  [1/5] BaseRecalibrator (before)..." \
        && gatk --java-options "-Xmx8g" BaseRecalibrator \
            -R "${REF}" -I "${BAM}" --known-sites "${KNOWN}" \
            -O "${BEFORE}" --tmp-dir "${SAMPLE_TMP}" \
        && echo "  [2/5] ApplyBQSR..." \
        && gatk --java-options "-Xmx8g" ApplyBQSR \
            -R "${REF}" -I "${BAM}" --bqsr-recal-file "${BEFORE}" \
            -O "${RECAL}" --create-output-bam-index false --tmp-dir "${SAMPLE_TMP}" \
        && echo "  [3/5] samtools index..." \
        && samtools index "${RECAL}" \
        && echo "  [4/5] BaseRecalibrator (after)..." \
        && gatk --java-options "-Xmx8g" BaseRecalibrator \
            -R "${REF}" -I "${RECAL}" --known-sites "${KNOWN}" \
            -O "${AFTER}" --tmp-dir "${SAMPLE_TMP}" \
        && echo "  [5/5] AnalyzeCovariates..." \
        && gatk AnalyzeCovariates \
            -before "${BEFORE}" -after "${AFTER}" \
            -csv "${COVAR}.csv" ${PLOT_ARGS[@]+"${PLOT_ARGS[@]}"} \
        && rm -rf "${SAMPLE_TMP:?}"; then
        echo "  OK: ${SAMPLE}"
        OK=$(( OK + 1 ))
    else
        echo "  ERROR: ${SAMPLE}" >&2
        rm -f "${RECAL}" "${RECAL}.bai" "${AFTER}"
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

sed -i \
    -e "s|__LAST_TASK__|$(( NUM_JOBS - 1 ))|g" \
    -e "s|__THREADS__|${THREADS}|g" \
    -e "s|__OUTPUT_DIR__|${OUTPUT_DIR}|g" \
    -e "s|__REF__|${REF}|g" \
    -e "s|__KNOWN__|${KNOWN}|g" \
    -e "s|__LIST_FILE__|${LIST_FILE}|g" \
    -e "s|__BAM_DIR__|${BAM_DIR}|g" \
    -e "s|__TABLE_DIR__|${TABLE_DIR}|g" \
    -e "s|__NUM_JOBS__|${NUM_JOBS}|g" \
    "${SCRIPT1}"

# ---------------------------------------------------------------------------
# Submission
# ---------------------------------------------------------------------------
OUTPUT=$(sbatch "${SCRIPT1}")
JOB1_ID=$(echo "${OUTPUT}" | awk '{print $NF}')

echo "============================================="
echo " Summary"
echo "============================================="
echo " Reference:    ${REF}"
echo " Known sites:  ${KNOWN} (${N_KNOWN} variants)"
for DIR in "${INPUT_DIRS[@]}"; do
echo " Input:        $(realpath "${DIR}")"
done
echo " Samples:      ${TOTAL} (${LIST_FILE})"
echo " Output BAM:   ${BAM_DIR}"
echo " Tables:       ${TABLE_DIR}"
echo " Threads/job:  ${THREADS}"
echo "---------------------------------------------"
echo " Job 1 (${JOB1_ID}): BQSR, array 0-$(( NUM_JOBS - 1 ))"
echo "---------------------------------------------"
echo " Total jobs: ${NUM_JOBS} / ${MAX_JOBS}"
echo "============================================="
echo
echo "Useful commands:"
echo "  squeue -u \$(whoami)"
echo "  tail -F ${OUTPUT_DIR}/logs/bqsr_*.out"
echo "  grep -h 'OK:\|ERROR' ${OUTPUT_DIR}/logs/bqsr_*.out"
echo "  grep -l 'ERROR' ${OUTPUT_DIR}/logs/*.err"
echo
echo "Next step (final call set on the recalibrated BAMs):"
echo "  $(dirname "$0")/run_variant_calling.sh -r ${REF} -i ${BAM_DIR} -o <OUTPUT>"
