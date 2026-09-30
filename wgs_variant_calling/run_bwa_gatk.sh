#!/bin/bash
# =============================================================================
# run_bwa_gatk.sh — Mapping with bwa-mem2 + MarkDuplicates (HPC Drago, SLURM)
#
# Step 1: reference indexes (bwa-mem2, samtools faidx, .dict) (single job)
# Step 2: bwa-mem2 mem | samtools sort + MarkDuplicates + flagstat (array)
#
# Output per sample: SAMPLE_sorted_dedup.bam(.bai), SAMPLE_dup_metrics.txt,
#                     SAMPLE_flagstat.txt
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
REF=""
INPUT_DIRS=()
OUTPUT_DIR=""
MAX_JOBS=60
THREADS=8
START_STEP=1

usage() {
    cat << USAGE
Usage: $(basename "$0") -r REF.fa -i DIR_FASTQ [-i DIR2 ...] -o OUTPUT_DIR [-j MAX_JOBS] [-t THREADS] [-s STEP] [-h]

Maps paired-end reads with bwa-mem2, sorts them, marks duplicates with GATK MarkDuplicates
and indexes the BAM. The read group uses the sample name (ID=SM=LB=SAMPLE).
All work is submitted to SLURM; this script only writes and submits the jobs.

Flags:
  -r  Reference FASTA (required). Indexes are created next to it if missing
  -i  Directory with FASTQ files (required, repeatable). Searched in this order:
        *_R1.trimmed.fastq.gz / *_1.trimmed.fastq.gz   (output of run_trim_galore.sh)
        *_R1.fastq.gz / *_1.fastq.gz                   (untrimmed, with a warning)
  -o  Output directory (required). BAM in OUTPUT_DIR/bam, logs in OUTPUT_DIR/logs
  -j  Maximum number of queued jobs (default: ${MAX_JOBS} = MaxJobsPU)
  -t  Threads per job (default: ${THREADS})
  -s  Step to start from (default: 1)
        1 = reference indexes
        2 = mapping + MarkDuplicates
  -h  Show this help

Examples:
  conda activate gatk
  $(basename "$0") -r /lustre/home/iata/aaguilar/variant_GATK/nuclear.fasta \\
      -i /lustre/home/iata/aaguilar/paper/data/fastq -o /lustre/home/iata/aaguilar/paper/mapping
  $(basename "$0") -r ref.fa -i fastq_1 -i fastq_2 -o mapping -t 16
  $(basename "$0") -r ref.fa -i fastq -o mapping -s 2      # indexes already built

Samples whose deduplicated BAM and index already exist are skipped when re-launching.
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
if ! [[ "${MAX_JOBS}" =~ ^[0-9]+$ && "${THREADS}" =~ ^[0-9]+$ && "${START_STEP}" =~ ^[12]$ ]]; then
    echo "ERROR: -j and -t must be integers and -s must be 1 or 2."
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
for tool in bwa-mem2 samtools; do
    command -v "${tool}" > /dev/null || { echo "ERROR: '${tool}' is not in the 'gatk' environment."; exit 1; }
done

# ---------------------------------------------------------------------------
# Directories
# ---------------------------------------------------------------------------
REF=$(realpath "${REF}")
OUTPUT_DIR=$(realpath -m "${OUTPUT_DIR}")
BAM_DIR="${OUTPUT_DIR}/bam"
mkdir -p "${BAM_DIR}" "${OUTPUT_DIR}/logs" "${OUTPUT_DIR}/tmp"

if [[ "${START_STEP}" -ge 2 ]]; then
    REF_DICT="${REF%.*}.dict"
    for f in "${REF}.fai" "${REF}.bwt.2bit.64" "${REF_DICT}"; do
        [[ -f "${f}" ]] || { echo "ERROR: ${f} not found. Start from step 1 (-s 1)."; exit 1; }
    done
fi

# ---------------------------------------------------------------------------
# Sample list: SAMPLE<TAB>R1<TAB>R2
# ---------------------------------------------------------------------------
LIST_FILE="${OUTPUT_DIR}/samples.tsv"
: > "${LIST_FILE}"
N_RAW=0
for DIR in "${INPUT_DIRS[@]}"; do
    DIR=$(realpath "${DIR}")
    R1_FILES=$(find "${DIR}" -type f \( -name "*_1.trimmed.fastq.gz" -o -name "*_R1.trimmed.fastq.gz" \) | sort)
    if [[ -z "${R1_FILES}" ]]; then
        R1_FILES=$(find "${DIR}" -type f \( -name "*_1.fastq.gz" -o -name "*_R1.fastq.gz" \) | sort)
        [[ -n "${R1_FILES}" ]] && N_RAW=$(( N_RAW + $(echo "${R1_FILES}" | wc -l) ))
    fi
    while read -r R1; do
        [[ -z "${R1}" ]] && continue
        BASE=$(basename "${R1}")
        case "${BASE}" in
            *_R1.trimmed.fastq.gz) SAMPLE=${BASE%_R1.trimmed.fastq.gz}; R2=${R1%_R1.trimmed.fastq.gz}_R2.trimmed.fastq.gz ;;
            *_1.trimmed.fastq.gz)  SAMPLE=${BASE%_1.trimmed.fastq.gz};  R2=${R1%_1.trimmed.fastq.gz}_2.trimmed.fastq.gz ;;
            *_R1.fastq.gz)         SAMPLE=${BASE%_R1.fastq.gz};         R2=${R1%_R1.fastq.gz}_R2.fastq.gz ;;
            *_1.fastq.gz)          SAMPLE=${BASE%_1.fastq.gz};          R2=${R1%_1.fastq.gz}_2.fastq.gz ;;
        esac
        if [[ ! -f "${R2}" ]]; then
            echo "WARNING: no R2 for ${R1}, skipped."
            continue
        fi
        if cut -f1 "${LIST_FILE}" | grep -qxF "${SAMPLE}"; then
            echo "ERROR: sample duplicated across directories: ${SAMPLE}"
            exit 1
        fi
        printf '%s\t%s\t%s\n' "${SAMPLE}" "${R1}" "${R2}" >> "${LIST_FILE}"
    done <<< "${R1_FILES}"
done

TOTAL=$(wc -l < "${LIST_FILE}")
[[ ${TOTAL} -gt 0 ]] || { echo "ERROR: no R1/R2 pairs found in ${INPUT_DIRS[*]}"; exit 1; }
[[ ${N_RAW} -gt 0 ]] && echo "WARNING: ${N_RAW} untrimmed samples (*.fastq.gz without .trimmed). Was run_trim_galore.sh run?"

# Step 1 is a single job: it is subtracted from the array pool
NUM_SEQ_JOBS=0
[[ "${START_STEP}" -le 1 ]] && NUM_SEQ_JOBS=1
NUM_JOBS=$(( MAX_JOBS - NUM_SEQ_JOBS ))
[[ ${NUM_JOBS} -gt ${TOTAL} ]] && NUM_JOBS=${TOTAL}
[[ ${NUM_JOBS} -lt 1 ]] && { echo "ERROR: -j demasiado bajo."; exit 1; }

# ---------------------------------------------------------------------------
# SLURM script — Step 1: indexes (single job)
# ---------------------------------------------------------------------------
SCRIPT1="${OUTPUT_DIR}/index_ref.slurm"
cat > "${SCRIPT1}" << EOF
#!/bin/bash
#SBATCH --job-name=index_ref
#SBATCH --partition=generic
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=${OUTPUT_DIR}/logs/index_ref_%j.out
#SBATCH --error=${OUTPUT_DIR}/logs/index_ref_%j.err

set -eo pipefail
module purge
module load rama0.4 GCCcore/13.2.0 GATK/4.6.0.0-Java-17
eval "\$(conda shell.bash hook)"
conda activate gatk
set -u

REF="${REF}"
DICT="\${REF%.*}.dict"

echo "=== Reference indexes ==="
echo "Job ID: \${SLURM_JOB_ID}"
echo "Reference: \${REF}"
echo "Start: \$(date)"
echo "========================"

echo "  [1/3] samtools faidx..."
[[ -f "\${REF}.fai" ]] && echo "    already exists" || samtools faidx "\${REF}"
echo "  [2/3] CreateSequenceDictionary..."
[[ -f "\${DICT}" ]] && echo "    already exists" || gatk CreateSequenceDictionary -R "\${REF}" -O "\${DICT}"
echo "  [3/3] bwa-mem2 index..."
[[ -f "\${REF}.bwt.2bit.64" ]] && echo "    already exists" || bwa-mem2 index "\${REF}"

echo
echo "========================"
echo "  End: \$(date)"
echo "========================"
EOF

# ---------------------------------------------------------------------------
# SLURM script — Step 2: mapping + MarkDuplicates (array)
# ---------------------------------------------------------------------------
SCRIPT2="${OUTPUT_DIR}/bwa_markdup.slurm"
cat > "${SCRIPT2}" << 'EOF'
#!/bin/bash
#SBATCH --job-name=bwa_markdup
#SBATCH --partition=generic
#SBATCH --array=0-__LAST_TASK__
#SBATCH --cpus-per-task=__THREADS__
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=__OUTPUT_DIR__/logs/bwa_markdup_%A_%a.out
#SBATCH --error=__OUTPUT_DIR__/logs/bwa_markdup_%A_%a.err

set -eo pipefail
module purge
module load rama0.4 GCCcore/13.2.0 GATK/4.6.0.0-Java-17
eval "$(conda shell.bash hook)"
conda activate gatk
set -u

REF="__REF__"
LIST_FILE="__LIST_FILE__"
BAM_DIR="__BAM_DIR__"
TMP_DIR="__OUTPUT_DIR__/tmp"
THREADS=__THREADS__
NUM_JOBS=__NUM_JOBS__
TOTAL=$(wc -l < "${LIST_FILE}")

echo "=== bwa-mem2 mapping + MarkDuplicates ==="
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
while IFS=$'\t' read -r SAMPLE R1 R2; do
    N=$(( N + 1 ))
    PREFIX="${BAM_DIR}/${SAMPLE}"
    SORTED="${PREFIX}_sorted_tmp.bam"
    DEDUP="${PREFIX}_sorted_dedup.bam"
    SAMPLE_TMP="${TMP_DIR}/${SAMPLE}"
    echo
    echo "[${N}/${N_BATCH}] Processing: ${SAMPLE}"
    echo "  R1: ${R1}"
    echo "  R2: ${R2}"
    echo "  Time: $(date +%H:%M:%S)"

    if [[ -s "${DEDUP}" && -s "${DEDUP}.bai" ]]; then
        echo "  Already mapped, skipped."
        OK=$(( OK + 1 ))
        continue
    fi
    mkdir -p "${SAMPLE_TMP}"
    RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA\tPU:${SAMPLE}"

    # bwa-mem2 exits with code 0 even if a FASTQ is empty or truncated:
    # reads are counted before mapping and compared with the BAM afterwards
    N_R1=$(( $(zcat "${R1}" | wc -l) / 4 ))
    N_R2=$(( $(zcat "${R2}" | wc -l) / 4 ))
    echo "  Reads: R1=${N_R1} R2=${N_R2}"
    if [[ ${N_R1} -eq 0 || ${N_R1} -ne ${N_R2} ]]; then
        echo "  ERROR: ${SAMPLE}: R1 and R2 empty or with different read numbers" >&2
        FAILED=$(( FAILED + 1 ))
        continue
    fi

    # set -e does not act inside an 'if': steps are chained with && to stop at the first failure
    if echo "  [1/5] bwa-mem2 mem | samtools sort..." \
        && bwa-mem2 mem -t "${THREADS}" -R "${RG}" "${REF}" "${R1}" "${R2}" \
            | samtools sort -@ "${THREADS}" -T "${SAMPLE_TMP}/sort" -o "${SORTED}" - \
        && echo "  [2/5] samtools index..." \
        && samtools index "${SORTED}" \
        && echo "  [3/5] MarkDuplicates..." \
        && gatk MarkDuplicates \
            -I "${SORTED}" \
            -O "${DEDUP}" \
            -M "${PREFIX}_dup_metrics.txt" \
            --READ_NAME_REGEX null \
            --ASSUME_SORT_ORDER coordinate \
            --TMP_DIR "${SAMPLE_TMP}" \
        && echo "  [4/5] samtools index + flagstat..." \
        && samtools index "${DEDUP}" \
        && samtools flagstat -@ "${THREADS}" "${DEDUP}" > "${PREFIX}_flagstat.txt" \
        && N_PRIMARY=$(awk '/ primary$/ {print $1}' "${PREFIX}_flagstat.txt") \
        && { [[ "${N_PRIMARY}" -eq $(( N_R1 + N_R2 )) ]] \
             || { echo "  ERROR: the BAM has ${N_PRIMARY} primary reads, $(( N_R1 + N_R2 )) expected" >&2; false; }; } \
        && echo "  [5/5] removing intermediate files..." \
        && rm -f "${SORTED}" "${SORTED}.bai" \
        && rm -rf "${SAMPLE_TMP:?}"; then
        echo "  Mapped: $(awk '/primary mapped/ {print $1, $6}' "${PREFIX}_flagstat.txt" | tr -d '(')"
        echo "  OK: ${SAMPLE}"
        OK=$(( OK + 1 ))
    else
        echo "  ERROR: ${SAMPLE}" >&2
        rm -f "${DEDUP}" "${DEDUP}.bai" "${SORTED}" "${SORTED}.bai"
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
    -e "s|__LIST_FILE__|${LIST_FILE}|g" \
    -e "s|__BAM_DIR__|${BAM_DIR}|g" \
    -e "s|__NUM_JOBS__|${NUM_JOBS}|g" \
    "${SCRIPT2}"

# ---------------------------------------------------------------------------
# Submission
# ---------------------------------------------------------------------------
PREV_DEP=""
JOB_LINES=()
N_SUBMITTED=0
if [[ "${START_STEP}" -le 1 ]]; then
    OUTPUT=$(sbatch ${PREV_DEP} "${SCRIPT1}")
    JOB1_ID=$(echo "${OUTPUT}" | awk '{print $NF}')
    PREV_DEP="--dependency=afterok:${JOB1_ID}"
    JOB_LINES+=(" Job 1 (${JOB1_ID}): reference indexes")
    N_SUBMITTED=$(( N_SUBMITTED + 1 ))
fi
if [[ "${START_STEP}" -le 2 ]]; then
    OUTPUT=$(sbatch ${PREV_DEP} "${SCRIPT2}")
    JOB2_ID=$(echo "${OUTPUT}" | awk '{print $NF}')
    JOB_LINES+=(" Job 2 (${JOB2_ID}): bwa-mem2 + MarkDuplicates, array 0-$(( NUM_JOBS - 1 ))${PREV_DEP:+ (${PREV_DEP#--dependency=})}")
    N_SUBMITTED=$(( N_SUBMITTED + NUM_JOBS ))
fi

echo "============================================="
echo " Summary"
echo "============================================="
echo " Reference:    ${REF}"
for DIR in "${INPUT_DIRS[@]}"; do
echo " Input:        $(realpath "${DIR}")"
done
echo " Samples:      ${TOTAL} (${LIST_FILE})"
echo " Output BAM:   ${BAM_DIR}"
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
echo "  tail -F ${OUTPUT_DIR}/logs/bwa_markdup_*.out"
echo "  grep -l 'ERROR' ${OUTPUT_DIR}/logs/*.err"
echo "  grep -H 'primary mapped' ${BAM_DIR}/*_flagstat.txt"
