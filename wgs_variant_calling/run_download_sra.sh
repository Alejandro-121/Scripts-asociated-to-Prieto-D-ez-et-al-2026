#!/bin/bash
# =============================================================================
# run_download_sra.sh — Download of WGS reads from SRA (HPC Drago, SLURM)
#
# Step 1: prefetch + fasterq-dump + pigz (SLURM array)
# Step 2: R1/R2 pairing check and read summary (single job)
#
# By default it downloads BioProject PRJNA1418127 (Prieto-Díez et al. 2026), naming the
# FASTQ files after the sample names used in the VCFs: SAMPLE_R1.fastq.gz / SAMPLE_R2.fastq.gz
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
OUTPUT_DIR=""
ACC_FILE=""
MAX_JOBS=60
THREADS=8
START_STEP=1
MODULES="rama0.4 SRA-Toolkit/3.2.0"   # prefetch, fasterq-dump
CONDA_ENV="download"                  # pigz

usage() {
    cat << USAGE
Usage: $(basename "$0") -o OUTPUT_DIR [-a ACCESSIONS.tsv] [-j MAX_JOBS] [-t THREADS] [-s STEP] [-h]

Downloads paired-end reads from SRA with prefetch + fasterq-dump and compresses them with pigz.
All work is submitted to SLURM; this script only writes and submits the jobs.

Flags:
  -o  Output directory (required). FASTQ in OUTPUT_DIR/fastq, logs in OUTPUT_DIR/logs
  -a  TSV with two columns: RUN<TAB>SAMPLE (optional).
      Default: the 11 samples of PRJNA1418127
  -j  Maximum number of queued jobs (default: ${MAX_JOBS} = MaxJobsPU)
  -t  Threads per job for fasterq-dump and pigz (default: ${THREADS})
  -s  Step to start from (default: 1)
        1 = download
        2 = R1/R2 check
  -h  Show this help

Examples:
  conda activate ${CONDA_ENV}
  $(basename "$0") -o /lustre/home/iata/aaguilar/eif5a/sra
  $(basename "$0") -o /lustre/home/iata/aaguilar/eif5a/sra -t 16
  $(basename "$0") -o /lustre/home/iata/aaguilar/eif5a/sra -s 2      # check only

Environment: modules '${MODULES}' (SRA Toolkit) + conda '${CONDA_ENV}' (pigz).
Note: compute nodes need internet access for prefetch.
USAGE
}

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
while getopts "o:a:j:t:s:h" opt; do
    case ${opt} in
        o) OUTPUT_DIR=${OPTARG} ;;
        a) ACC_FILE=${OPTARG} ;;
        j) MAX_JOBS=${OPTARG} ;;
        t) THREADS=${OPTARG} ;;
        s) START_STEP=${OPTARG} ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

if [[ -z "${OUTPUT_DIR}" ]]; then
    echo "ERROR: output directory missing (-o)."
    usage
    exit 1
fi
if ! [[ "${MAX_JOBS}" =~ ^[0-9]+$ && "${THREADS}" =~ ^[0-9]+$ && "${START_STEP}" =~ ^[12]$ ]]; then
    echo "ERROR: -j and -t must be integers and -s must be 1 or 2."
    exit 1
fi

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
if [[ "${CONDA_DEFAULT_ENV:-}" != "${CONDA_ENV}" ]]; then
    echo "WARNING: activate the conda environment '${CONDA_ENV}' before launching."
    echo "  conda activate ${CONDA_ENV}"
    exit 1
fi
if ! command -v pigz > /dev/null; then
    echo "ERROR: 'pigz' is not in the '${CONDA_ENV}' environment. Install it with:"
    echo "  conda install -n ${CONDA_ENV} -c conda-forge pigz"
    exit 1
fi
if ! type module > /dev/null 2>&1; then
    echo "ERROR: the 'module' command is not available in this session."
    exit 1
fi
# Checked in a subshell so that the modules of the session are not changed
if ! ( module purge && module load ${MODULES} && command -v prefetch && command -v fasterq-dump ) > /dev/null 2>&1; then
    echo "ERROR: could not load prefetch/fasterq-dump with: module load ${MODULES}"
    echo "  module spider SRA-Toolkit"
    exit 1
fi

# ---------------------------------------------------------------------------
# Directories and accession list
# ---------------------------------------------------------------------------
OUTPUT_DIR=$(realpath -m "${OUTPUT_DIR}")
FASTQ_DIR="${OUTPUT_DIR}/fastq"
mkdir -p "${FASTQ_DIR}" "${OUTPUT_DIR}/logs" "${OUTPUT_DIR}/tmp" "${OUTPUT_DIR}/sra_cache"

LIST_FILE="${OUTPUT_DIR}/accessions.tsv"
if [[ -n "${ACC_FILE}" ]]; then
    [[ -f "${ACC_FILE}" ]] || { echo "ERROR: ${ACC_FILE} not found"; exit 1; }
    grep -v -e '^#' -e '^[[:space:]]*$' "${ACC_FILE}" > "${LIST_FILE}"
else
    # PRJNA1418127 — RUN<TAB>SAMPLE (name used in the VCFs)
    cat > "${LIST_FILE}" << 'ACC'
SRR37083327	WT
SRR37083326	2-1
SRR37083324	2-3
SRR37083323	sup.1
SRR37083322	sup.2
SRR37083321	sup.22
SRR37083320	sup.23
SRR37083319	sup.11
SRR37083318	sup.15
SRR37083317	sup.25
SRR37083325	sup.27
ACC
fi

TOTAL=$(wc -l < "${LIST_FILE}")
[[ ${TOTAL} -gt 0 ]] || { echo "ERROR: the accession list is empty."; exit 1; }

# Step 2 is a single job: it is subtracted from the array pool
NUM_SEQ_JOBS=1
NUM_JOBS=$(( MAX_JOBS - NUM_SEQ_JOBS ))
[[ ${NUM_JOBS} -gt ${TOTAL} ]] && NUM_JOBS=${TOTAL}
[[ ${NUM_JOBS} -lt 1 ]] && { echo "ERROR: -j demasiado bajo."; exit 1; }

# ---------------------------------------------------------------------------
# SLURM script — Step 1: download (array)
# ---------------------------------------------------------------------------
SCRIPT1="${OUTPUT_DIR}/download_sra.slurm"
cat > "${SCRIPT1}" << 'EOF'
#!/bin/bash
#SBATCH --job-name=sra_dl
#SBATCH --partition=generic
#SBATCH --array=0-__LAST_TASK__
#SBATCH --cpus-per-task=__THREADS__
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH --output=__OUTPUT_DIR__/logs/sra_dl_%A_%a.out
#SBATCH --error=__OUTPUT_DIR__/logs/sra_dl_%A_%a.err

set -eo pipefail
module purge
module load __MODULES__
eval "$(conda shell.bash hook)"
conda activate __CONDA_ENV__
set -u

LIST_FILE="__LIST_FILE__"
FASTQ_DIR="__FASTQ_DIR__"
CACHE_DIR="__OUTPUT_DIR__/sra_cache"
TMP_DIR="__OUTPUT_DIR__/tmp"
THREADS=__THREADS__
NUM_JOBS=__NUM_JOBS__
TOTAL=$(wc -l < "${LIST_FILE}")

echo "=== SRA download ==="
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID}, Task ID: ${SLURM_ARRAY_TASK_ID}"
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
while IFS=$'\t' read -r RUN SAMPLE; do
    N=$(( N + 1 ))
    R1="${FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
    R2="${FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"
    echo
    echo "[${N}/${N_BATCH}] Processing: ${SAMPLE} (${RUN})"
    echo "  Time: $(date +%H:%M:%S)"

    if [[ -s "${R1}" && -s "${R2}" ]]; then
        echo "  Already downloaded, skipped."
        OK=$(( OK + 1 ))
        continue
    fi

    # set -e does not act inside an 'if': steps are chained with && to stop at the first failure
    if echo "  [1/4] prefetch..." \
        && prefetch "${RUN}" --output-directory "${CACHE_DIR}" --max-size 20G \
        && echo "  [2/4] fasterq-dump..." \
        && fasterq-dump "${CACHE_DIR}/${RUN}" --split-files --threads "${THREADS}" \
            --outdir "${TMP_DIR}/${RUN}" --temp "${TMP_DIR}/${RUN}" \
        && echo "  [3/4] pigz..." \
        && pigz -p "${THREADS}" "${TMP_DIR}/${RUN}/${RUN}_1.fastq" "${TMP_DIR}/${RUN}/${RUN}_2.fastq" \
        && echo "  [4/4] renaming and cleaning up..." \
        && mv "${TMP_DIR}/${RUN}/${RUN}_1.fastq.gz" "${R1}" \
        && mv "${TMP_DIR}/${RUN}/${RUN}_2.fastq.gz" "${R2}" \
        && rm -rf "${CACHE_DIR:?}/${RUN}" "${TMP_DIR:?}/${RUN}"; then
        echo "  OK: ${SAMPLE}"
        OK=$(( OK + 1 ))
    else
        echo "  ERROR: ${SAMPLE} (${RUN})" >&2
        rm -f "${R1}" "${R2}"
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
    -e "s|__LIST_FILE__|${LIST_FILE}|g" \
    -e "s|__FASTQ_DIR__|${FASTQ_DIR}|g" \
    -e "s|__NUM_JOBS__|${NUM_JOBS}|g" \
    -e "s|__CONDA_ENV__|${CONDA_ENV}|g" \
    -e "s|__MODULES__|${MODULES}|g" \
    "${SCRIPT1}"

# ---------------------------------------------------------------------------
# SLURM script — Step 2: R1/R2 check (single job)
# ---------------------------------------------------------------------------
SCRIPT2="${OUTPUT_DIR}/check_fastq.slurm"
cat > "${SCRIPT2}" << EOF
#!/bin/bash
#SBATCH --job-name=sra_check
#SBATCH --partition=generic
#SBATCH --cpus-per-task=${THREADS}
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH --output=${OUTPUT_DIR}/logs/sra_check_%j.out
#SBATCH --error=${OUTPUT_DIR}/logs/sra_check_%j.err

set -eo pipefail
eval "\$(conda shell.bash hook)"
conda activate ${CONDA_ENV}
set -u

SUMMARY="${OUTPUT_DIR}/fastq_summary.tsv"

echo "=== FASTQ check ==="
echo "Job ID: \${SLURM_JOB_ID}"
echo "Start: \$(date)"
echo "========================"

echo -e "run\tsample\treads_R1\treads_R2\tstatus" > "\${SUMMARY}"
OK=0
FAILED=0
while IFS=\$'\t' read -r RUN SAMPLE; do
    R1="${FASTQ_DIR}/\${SAMPLE}_R1.fastq.gz"
    R2="${FASTQ_DIR}/\${SAMPLE}_R2.fastq.gz"
    if [[ ! -s "\${R1}" || ! -s "\${R2}" ]]; then
        echo -e "\${RUN}\t\${SAMPLE}\tNA\tNA\tMISSING" >> "\${SUMMARY}"
        FAILED=\$(( FAILED + 1 ))
        continue
    fi
    N1=\$(( \$(pigz -dc -p ${THREADS} "\${R1}" | wc -l) / 4 ))
    N2=\$(( \$(pigz -dc -p ${THREADS} "\${R2}" | wc -l) / 4 ))
    if [[ \${N1} -eq \${N2} ]]; then
        STATUS=OK; OK=\$(( OK + 1 ))
    else
        STATUS=MISMATCH; FAILED=\$(( FAILED + 1 ))
    fi
    echo -e "\${RUN}\t\${SAMPLE}\t\${N1}\t\${N2}\t\${STATUS}" >> "\${SUMMARY}"
    echo "  \${SAMPLE}: R1=\${N1} R2=\${N2} \${STATUS}"
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
    JOB_LINES+=(" Job 1 (${JOB1_ID}): SRA download, array 0-$(( NUM_JOBS - 1 ))")
    N_SUBMITTED=$(( N_SUBMITTED + NUM_JOBS ))
fi
if [[ "${START_STEP}" -le 2 ]]; then
    OUTPUT=$(sbatch ${PREV_DEP} "${SCRIPT2}")
    JOB2_ID=$(echo "${OUTPUT}" | awk '{print $NF}')
    JOB_LINES+=(" Job 2 (${JOB2_ID}): R1/R2 check${PREV_DEP:+ (${PREV_DEP#--dependency=})}")
    N_SUBMITTED=$(( N_SUBMITTED + 1 ))
fi

echo "============================================="
echo " Summary"
echo "============================================="
echo " Accessions:   ${LIST_FILE} (${TOTAL} samples)"
echo " Output FASTQ: ${FASTQ_DIR}"
echo " Threads/job:  ${THREADS}"
echo " Start step:   ${START_STEP}"
echo " Environment:  module load ${MODULES} + conda ${CONDA_ENV}"
echo "---------------------------------------------"
printf '%s\n' "${JOB_LINES[@]}"
echo "---------------------------------------------"
echo " Total jobs: ${N_SUBMITTED} / ${MAX_JOBS}"
echo "============================================="
echo
echo "Useful commands:"
echo "  squeue -u \$(whoami)"
echo "  tail -f ${OUTPUT_DIR}/logs/sra_dl_*.out"
echo "  grep -l 'ERROR' ${OUTPUT_DIR}/logs/*.err"
echo "  column -t ${OUTPUT_DIR}/fastq_summary.tsv"
