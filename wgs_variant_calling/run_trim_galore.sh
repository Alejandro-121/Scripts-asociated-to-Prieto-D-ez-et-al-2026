#!/bin/bash
# =============================================================================
# run_trim_galore.sh — Read trimming with Trim Galore (HPC Drago, SLURM)
#
# Step 1: Trim Galore per sample (array)
#           --paired --nextera           Nextera adapters (Illumina DNA Prep)
#           --clip_R1/--clip_R2 N        hard trim of N bp at the 5' end (default 15)
#           --nextseq 20                 quality trimming that ignores signal-less G calls (2-colour)
#           --length 36                  discards pairs with a read < 36 bp
# Step 2: input/output read summary per sample (single job)
#
# Output per sample: trimmed/SAMPLE_R1.trimmed.fastq.gz / SAMPLE_R2.trimmed.fastq.gz
#                    reports/SAMPLE_R{1,2}.fastq.gz_trimming_report.txt
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
INPUT_DIRS=()
OUTPUT_DIR=""
MAX_JOBS=60
THREADS=8
START_STEP=1
CLIP=15
MODULES="rama0.4 GCCcore/12.3.0 Trim_Galore/0.6.10"   # module spider Trim_Galore/0.6.10

usage() {
    cat << USAGE
Usage: $(basename "$0") -i DIR_FASTQ [-i DIR2 ...] -o OUTPUT_DIR [-c CLIP] [-j MAX_JOBS] [-t THREADS] [-s STEP] [-h]

Trims paired-end reads with Trim Galore: Nextera adapters, hard trim of CLIP bp at the
5' end of R1 and R2, NextSeq quality trimming (--nextseq 20) and minimum length 36.
All work is submitted to SLURM; this script only writes and submits the jobs.

Flags:
  -i  Directory with untrimmed FASTQ files (required, repeatable):
        *_R1.fastq.gz / *_R2.fastq.gz   or   *_1.fastq.gz / *_2.fastq.gz
  -o  Output directory (required). FASTQ in OUTPUT_DIR/trimmed, logs in OUTPUT_DIR/logs
  -c  Base pairs to clip from the 5' end of R1 and R2 (default: ${CLIP})
  -j  Maximum number of queued jobs (default: ${MAX_JOBS} = MaxJobsPU)
  -t  Threads per job (default: ${THREADS}); Trim Galore uses --cores = THREADS/4
  -s  Step to start from (default: 1)
        1 = Trim Galore
        2 = read summary
  -h  Show this help

Environment: module load ${MODULES}

Examples:
  $(basename "$0") -i /lustre/home/iata/aaguilar/paper/data/fastq -o /lustre/home/iata/aaguilar/paper/trim
  $(basename "$0") -i fastq -o trim -c 10 -t 16
  $(basename "$0") -i fastq -o trim -s 2      # summary only

Samples that are already trimmed are skipped when re-launching.
USAGE
}

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
while getopts "i:o:c:j:t:s:h" opt; do
    case ${opt} in
        i) INPUT_DIRS+=("${OPTARG}") ;;
        o) OUTPUT_DIR=${OPTARG} ;;
        c) CLIP=${OPTARG} ;;
        j) MAX_JOBS=${OPTARG} ;;
        t) THREADS=${OPTARG} ;;
        s) START_STEP=${OPTARG} ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

if [[ -z "${OUTPUT_DIR}" || ${#INPUT_DIRS[@]} -eq 0 ]]; then
    echo "ERROR: missing required arguments (-i, -o)."
    usage
    exit 1
fi
if ! [[ "${MAX_JOBS}" =~ ^[0-9]+$ && "${THREADS}" =~ ^[0-9]+$ && "${CLIP}" =~ ^[0-9]+$ && "${START_STEP}" =~ ^[12]$ ]]; then
    echo "ERROR: -j, -t and -c must be integers and -s must be 1 or 2."
    exit 1
fi
for DIR in "${INPUT_DIRS[@]}"; do
    [[ -d "${DIR}" ]] || { echo "ERROR: directory ${DIR} not found"; exit 1; }
done

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
if ! type module > /dev/null 2>&1; then
    echo "ERROR: the 'module' command is not available in this session."
    exit 1
fi
# Checked in a subshell so that the modules of the session are not changed
if ! ( module purge && module load ${MODULES} && command -v trim_galore && command -v cutadapt ) > /dev/null 2>&1; then
    echo "ERROR: could not load trim_galore/cutadapt with: module load ${MODULES}"
    echo "  Diagnosis: module purge; module load ${MODULES}; which trim_galore cutadapt"
    echo "  Requisitos:  module spider ${MODULES##* }"
    exit 1
fi

# ---------------------------------------------------------------------------
# Directories
# ---------------------------------------------------------------------------
OUTPUT_DIR=$(realpath -m "${OUTPUT_DIR}")
TRIM_DIR="${OUTPUT_DIR}/trimmed"
mkdir -p "${TRIM_DIR}" "${OUTPUT_DIR}/reports" "${OUTPUT_DIR}/stats" "${OUTPUT_DIR}/logs" "${OUTPUT_DIR}/tmp"

CORES=$(( THREADS / 4 ))
[[ ${CORES} -lt 1 ]] && CORES=1

# ---------------------------------------------------------------------------
# Sample list: SAMPLE<TAB>R1<TAB>R2
# ---------------------------------------------------------------------------
LIST_FILE="${OUTPUT_DIR}/samples.tsv"
: > "${LIST_FILE}"
for DIR in "${INPUT_DIRS[@]}"; do
    DIR=$(realpath "${DIR}")
    while read -r R1; do
        [[ -z "${R1}" ]] && continue
        BASE=$(basename "${R1}")
        case "${BASE}" in
            *_R1.fastq.gz) SAMPLE=${BASE%_R1.fastq.gz}; R2=${R1%_R1.fastq.gz}_R2.fastq.gz ;;
            *_1.fastq.gz)  SAMPLE=${BASE%_1.fastq.gz};  R2=${R1%_1.fastq.gz}_2.fastq.gz ;;
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
    done < <(find "${DIR}" -maxdepth 1 \( -type f -o -type l \) \( -name "*_R1.fastq.gz" -o -name "*_1.fastq.gz" \) \
                 ! -name "*.trimmed.fastq.gz" | sort)
done

TOTAL=$(wc -l < "${LIST_FILE}")
[[ ${TOTAL} -gt 0 ]] || { echo "ERROR: no R1/R2 pairs found in ${INPUT_DIRS[*]}"; exit 1; }

# Step 2 is a single job: it is subtracted from the array pool
NUM_SEQ_JOBS=1
NUM_JOBS=$(( MAX_JOBS - NUM_SEQ_JOBS ))
[[ ${NUM_JOBS} -gt ${TOTAL} ]] && NUM_JOBS=${TOTAL}
[[ ${NUM_JOBS} -lt 1 ]] && { echo "ERROR: -j too low (minimum 2)."; exit 1; }

# ---------------------------------------------------------------------------
# SLURM script — Step 1: Trim Galore (array)
# ---------------------------------------------------------------------------
SCRIPT1="${OUTPUT_DIR}/trim_galore.slurm"
cat > "${SCRIPT1}" << 'EOF'
#!/bin/bash
#SBATCH --job-name=trim_galore
#SBATCH --partition=generic
#SBATCH --array=0-__LAST_TASK__
#SBATCH --cpus-per-task=__THREADS__
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH --output=__OUTPUT_DIR__/logs/trim_galore_%A_%a.out
#SBATCH --error=__OUTPUT_DIR__/logs/trim_galore_%A_%a.err

set -eo pipefail
module purge
module load __MODULES__
set -u

LIST_FILE="__LIST_FILE__"
TRIM_DIR="__TRIM_DIR__"
REPORT_DIR="__OUTPUT_DIR__/reports"
STATS_DIR="__OUTPUT_DIR__/stats"
TMP_DIR="__OUTPUT_DIR__/tmp"
CORES=__CORES__
CLIP=__CLIP__
NUM_JOBS=__NUM_JOBS__
TOTAL=$(wc -l < "${LIST_FILE}")

echo "=== Trim Galore ==="
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID}, Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "$(trim_galore --version | grep -i version | xargs) | cutadapt $(cutadapt --version 2>/dev/null || echo '-')"
echo "Hard trim 5': ${CLIP} pb | --nextera --nextseq 20 --length 36 | cores: ${CORES}"
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
    OUT_R1="${TRIM_DIR}/${SAMPLE}_R1.trimmed.fastq.gz"
    OUT_R2="${TRIM_DIR}/${SAMPLE}_R2.trimmed.fastq.gz"
    SAMPLE_TMP="${TMP_DIR}/${SAMPLE}"
    echo
    echo "[${N}/${N_BATCH}] Processing: ${SAMPLE}"
    echo "  R1: ${R1}"
    echo "  R2: ${R2}"
    echo "  Time: $(date +%H:%M:%S)"

    if [[ -s "${OUT_R1}" && -s "${OUT_R2}" && -s "${STATS_DIR}/${SAMPLE}.tsv" ]]; then
        echo "  Already trimmed, skipped."
        OK=$(( OK + 1 ))
        continue
    fi
    rm -rf "${SAMPLE_TMP}"
    mkdir -p "${SAMPLE_TMP}"

    # set -e does not act inside an 'if': steps are chained with && to stop at the first failure
    if echo "  [1/3] trim_galore..." \
        && trim_galore --paired --nextera \
            --clip_R1 "${CLIP}" --clip_R2 "${CLIP}" \
            --nextseq 20 --length 36 \
            --cores "${CORES}" \
            --basename "${SAMPLE}" \
            -o "${SAMPLE_TMP}" \
            "${R1}" "${R2}" \
        && echo "  [2/3] counting reads..." \
        && N_IN=$(( $(zcat "${R1}" | wc -l) / 4 )) \
        && N_OUT1=$(( $(zcat "${SAMPLE_TMP}/${SAMPLE}_val_1.fq.gz" | wc -l) / 4 )) \
        && N_OUT2=$(( $(zcat "${SAMPLE_TMP}/${SAMPLE}_val_2.fq.gz" | wc -l) / 4 )) \
        && { [[ ${N_OUT1} -gt 0 && ${N_OUT1} -eq ${N_OUT2} ]] \
             || { echo "  ERROR: empty output or unbalanced R1/R2 (${N_OUT1}/${N_OUT2})" >&2; false; }; } \
        && echo "  [3/3] moving results..." \
        && mv "${SAMPLE_TMP}/${SAMPLE}_val_1.fq.gz" "${OUT_R1}" \
        && mv "${SAMPLE_TMP}/${SAMPLE}_val_2.fq.gz" "${OUT_R2}" \
        && mv "${SAMPLE_TMP}"/*_trimming_report.* "${REPORT_DIR}/" \
        && printf '%s\t%s\t%s\n' "${SAMPLE}" "${N_IN}" "${N_OUT1}" > "${STATS_DIR}/${SAMPLE}.tsv" \
        && rm -rf "${SAMPLE_TMP:?}"; then
        echo "  Pairs: ${N_IN} -> ${N_OUT1} ($(awk -v a="${N_OUT1}" -v b="${N_IN}" 'BEGIN{printf "%.2f", 100*a/b}') % retained)"
        echo "  OK: ${SAMPLE}"
        OK=$(( OK + 1 ))
    else
        echo "  ERROR: ${SAMPLE}" >&2
        rm -f "${OUT_R1}" "${OUT_R2}" "${STATS_DIR}/${SAMPLE}.tsv"
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
    -e "s|__MODULES__|${MODULES}|g" \
    -e "s|__LIST_FILE__|${LIST_FILE}|g" \
    -e "s|__TRIM_DIR__|${TRIM_DIR}|g" \
    -e "s|__CORES__|${CORES}|g" \
    -e "s|__CLIP__|${CLIP}|g" \
    -e "s|__NUM_JOBS__|${NUM_JOBS}|g" \
    "${SCRIPT1}"

# ---------------------------------------------------------------------------
# SLURM script — Step 2: summary (single job)
# ---------------------------------------------------------------------------
SCRIPT2="${OUTPUT_DIR}/trim_summary.slurm"
cat > "${SCRIPT2}" << EOF
#!/bin/bash
#SBATCH --job-name=trim_summary
#SBATCH --partition=generic
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH --output=${OUTPUT_DIR}/logs/trim_summary_%j.out
#SBATCH --error=${OUTPUT_DIR}/logs/trim_summary_%j.err

set -euo pipefail

SUMMARY="${OUTPUT_DIR}/trim_summary.tsv"

echo "=== Trim Galore summary ==="
echo "Job ID: \${SLURM_JOB_ID}"
echo "Start: \$(date)"
echo "========================"

echo -e "sample\tpairs_in\tpairs_out\tpct_retained\tstatus" > "\${SUMMARY}"
OK=0
FAILED=0
while IFS=\$'\t' read -r SAMPLE R1 R2; do
    STATS="${OUTPUT_DIR}/stats/\${SAMPLE}.tsv"
    if [[ ! -s "\${STATS}" ]]; then
        echo -e "\${SAMPLE}\tNA\tNA\tNA\tMISSING" >> "\${SUMMARY}"
        FAILED=\$(( FAILED + 1 ))
        continue
    fi
    read -r _ N_IN N_OUT < "\${STATS}"
    PCT=\$(awk -v a="\${N_OUT}" -v b="\${N_IN}" 'BEGIN{printf "%.2f", 100*a/b}')
    echo -e "\${SAMPLE}\t\${N_IN}\t\${N_OUT}\t\${PCT}\tOK" >> "\${SUMMARY}"
    echo "  \${SAMPLE}: \${N_IN} -> \${N_OUT} (\${PCT} %)"
    OK=\$(( OK + 1 ))
done < "${LIST_FILE}"

echo
echo "========================"
echo "Summary:"
echo "  OK:        \${OK}"
echo "  Failed:    \${FAILED}"
echo "  Table:     \${SUMMARY}"
echo "  End: \$(date)"
echo "========================"
[[ \${FAILED} -eq 0 ]]
EOF

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
    JOB_LINES+=(" Job 1 (${JOB1_ID}): Trim Galore, array 0-$(( NUM_JOBS - 1 ))")
    N_SUBMITTED=$(( N_SUBMITTED + NUM_JOBS ))
fi
if [[ "${START_STEP}" -le 2 ]]; then
    OUTPUT=$(sbatch ${PREV_DEP} "${SCRIPT2}")
    JOB2_ID=$(echo "${OUTPUT}" | awk '{print $NF}')
    JOB_LINES+=(" Job 2 (${JOB2_ID}): read summary${PREV_DEP:+ (${PREV_DEP#--dependency=})}")
    N_SUBMITTED=$(( N_SUBMITTED + 1 ))
fi

echo "============================================="
echo " Summary"
echo "============================================="
for DIR in "${INPUT_DIRS[@]}"; do
echo " Input:        $(realpath "${DIR}")"
done
echo " Samples:      ${TOTAL} (${LIST_FILE})"
echo " Output FASTQ: ${TRIM_DIR}"
echo " Hard trim 5': ${CLIP} bp (R1 and R2)"
echo " Options:      --paired --nextera --nextseq 20 --length 36"
echo " Threads/job:  ${THREADS} (--cores ${CORES})"
echo " Environment:  module load ${MODULES}"
echo " Start step:   ${START_STEP}"
echo "---------------------------------------------"
printf '%s\n' "${JOB_LINES[@]}"
echo "---------------------------------------------"
echo " Total jobs: ${N_SUBMITTED} / ${MAX_JOBS}"
echo "============================================="
echo
echo "Useful commands:"
echo "  squeue -u \$(whoami)"
echo "  tail -F ${OUTPUT_DIR}/logs/trim_galore_*.out"
echo "  grep -l 'ERROR' ${OUTPUT_DIR}/logs/*.err"
echo "  column -t ${OUTPUT_DIR}/trim_summary.tsv"
