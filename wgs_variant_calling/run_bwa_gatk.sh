#!/bin/bash
# =============================================================================
# run_bwa_gatk.sh — Mapeo con bwa-mem2 + MarkDuplicates (HPC Drago, SLURM)
#
# Paso 1: índices de la referencia (bwa-mem2, samtools faidx, .dict) (job único)
# Paso 2: bwa-mem2 mem | samtools sort + MarkDuplicates + flagstat (array)
#
# Salida por muestra: SAMPLE_sorted_dedup.bam(.bai), SAMPLE_dup_metrics.txt,
#                     SAMPLE_flagstat.txt
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Valores por defecto
# ---------------------------------------------------------------------------
REF=""
INPUT_DIRS=()
OUTPUT_DIR=""
MAX_JOBS=60
THREADS=8
START_STEP=1

usage() {
    cat << USAGE
Uso: $(basename "$0") -r REF.fa -i DIR_FASTQ [-i DIR2 ...] -o OUTPUT_DIR [-j MAX_JOBS] [-t THREADS] [-s PASO] [-h]

Mapea reads paired-end con bwa-mem2, ordena, marca duplicados con GATK MarkDuplicates
e indexa. El read group usa el nombre de la muestra (ID=SM=LB=SAMPLE).
Todo el trabajo se lanza a SLURM; este script solo genera y envía los jobs.

Flags:
  -r  Referencia FASTA (obligatorio). Los índices se crean junto a ella si faltan
  -i  Directorio con FASTQ (obligatorio, repetible). Busca por orden de preferencia:
        *_R1.trimmed.fastq.gz / *_1.trimmed.fastq.gz   (salida de run_trim_galore.sh)
        *_R1.fastq.gz / *_1.fastq.gz                   (sin recortar, con aviso)
  -o  Directorio de salida (obligatorio). BAM en OUTPUT_DIR/bam, logs en OUTPUT_DIR/logs
  -j  Límite total de jobs en cola (default: ${MAX_JOBS} = MaxJobsPU)
  -t  Hilos por job (default: ${THREADS})
  -s  Paso desde el que empezar (default: 1)
        1 = índices de la referencia
        2 = mapeo + MarkDuplicates
  -h  Muestra esta ayuda

Ejemplos:
  conda activate gatk
  $(basename "$0") -r /lustre/home/iata/aaguilar/variant_GATK/nuclear.fasta \\
      -i /lustre/home/iata/aaguilar/paper/data/fastq -o /lustre/home/iata/aaguilar/paper/mapping
  $(basename "$0") -r ref.fa -i fastq_1 -i fastq_2 -o mapping -t 16
  $(basename "$0") -r ref.fa -i fastq -o mapping -s 2      # índices ya creados

Las muestras con BAM dedup e índice ya presentes se omiten al relanzar.
USAGE
}

# ---------------------------------------------------------------------------
# Argumentos
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
    echo "ERROR: faltan argumentos obligatorios (-r, -i, -o)."
    usage
    exit 1
fi
if ! [[ "${MAX_JOBS}" =~ ^[0-9]+$ && "${THREADS}" =~ ^[0-9]+$ && "${START_STEP}" =~ ^[12]$ ]]; then
    echo "ERROR: -j y -t deben ser enteros y -s debe ser 1 o 2."
    exit 1
fi
[[ -f "${REF}" ]] || { echo "ERROR: no existe la referencia ${REF}"; exit 1; }
for DIR in "${INPUT_DIRS[@]}"; do
    [[ -d "${DIR}" ]] || { echo "ERROR: no existe el directorio ${DIR}"; exit 1; }
done

# ---------------------------------------------------------------------------
# Entorno
# ---------------------------------------------------------------------------
if [[ "${CONDA_DEFAULT_ENV:-}" != "gatk" ]]; then
    echo "AVISO: Activa el entorno conda 'gatk' antes de lanzar."
    echo "  conda activate gatk"
    exit 1
fi
for tool in bwa-mem2 samtools; do
    command -v "${tool}" > /dev/null || { echo "ERROR: '${tool}' no está en el entorno 'gatk'."; exit 1; }
done

# ---------------------------------------------------------------------------
# Directorios
# ---------------------------------------------------------------------------
REF=$(realpath "${REF}")
OUTPUT_DIR=$(realpath -m "${OUTPUT_DIR}")
BAM_DIR="${OUTPUT_DIR}/bam"
mkdir -p "${BAM_DIR}" "${OUTPUT_DIR}/logs" "${OUTPUT_DIR}/tmp"

if [[ "${START_STEP}" -ge 2 ]]; then
    REF_DICT="${REF%.*}.dict"
    for f in "${REF}.fai" "${REF}.bwt.2bit.64" "${REF_DICT}"; do
        [[ -f "${f}" ]] || { echo "ERROR: falta ${f}. Lanza desde el paso 1 (-s 1)."; exit 1; }
    done
fi

# ---------------------------------------------------------------------------
# Lista de muestras: SAMPLE<TAB>R1<TAB>R2
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
            echo "AVISO: sin R2 para ${R1}, se omite."
            continue
        fi
        if cut -f1 "${LIST_FILE}" | grep -qxF "${SAMPLE}"; then
            echo "ERROR: muestra duplicada entre directorios: ${SAMPLE}"
            exit 1
        fi
        printf '%s\t%s\t%s\n' "${SAMPLE}" "${R1}" "${R2}" >> "${LIST_FILE}"
    done <<< "${R1_FILES}"
done

TOTAL=$(wc -l < "${LIST_FILE}")
[[ ${TOTAL} -gt 0 ]] || { echo "ERROR: no se encontraron pares R1/R2 en ${INPUT_DIRS[*]}"; exit 1; }
[[ ${N_RAW} -gt 0 ]] && echo "AVISO: ${N_RAW} muestras sin recortar (*.fastq.gz sin .trimmed). ¿Falta run_trim_galore.sh?"

# El paso 1 es un job secuencial: se resta del pool de arrays
NUM_SEQ_JOBS=0
[[ "${START_STEP}" -le 1 ]] && NUM_SEQ_JOBS=1
NUM_JOBS=$(( MAX_JOBS - NUM_SEQ_JOBS ))
[[ ${NUM_JOBS} -gt ${TOTAL} ]] && NUM_JOBS=${TOTAL}
[[ ${NUM_JOBS} -lt 1 ]] && { echo "ERROR: -j demasiado bajo."; exit 1; }

# ---------------------------------------------------------------------------
# Script SLURM — Paso 1: índices (job único)
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

echo "=== Índices de la referencia ==="
echo "Job ID: \${SLURM_JOB_ID}"
echo "Referencia: \${REF}"
echo "Inicio: \$(date)"
echo "========================"

echo "  [1/3] samtools faidx..."
[[ -f "\${REF}.fai" ]] && echo "    ya existe" || samtools faidx "\${REF}"
echo "  [2/3] CreateSequenceDictionary..."
[[ -f "\${DICT}" ]] && echo "    ya existe" || gatk CreateSequenceDictionary -R "\${REF}" -O "\${DICT}"
echo "  [3/3] bwa-mem2 index..."
[[ -f "\${REF}.bwt.2bit.64" ]] && echo "    ya existe" || bwa-mem2 index "\${REF}"

echo
echo "========================"
echo "  Fin: \$(date)"
echo "========================"
EOF

# ---------------------------------------------------------------------------
# Script SLURM — Paso 2: mapeo + MarkDuplicates (array)
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

echo "=== Mapeo bwa-mem2 + MarkDuplicates ==="
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID}, Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Referencia: ${REF}"
echo "Inicio: $(date)"
echo "========================"

ITEMS_PER_JOB=$(( (TOTAL + NUM_JOBS - 1) / NUM_JOBS ))
START_LINE=$(( SLURM_ARRAY_TASK_ID * ITEMS_PER_JOB + 1 ))
END_LINE=$(( START_LINE + ITEMS_PER_JOB - 1 ))
[[ ${END_LINE} -gt ${TOTAL} ]] && END_LINE=${TOTAL}
[[ ${START_LINE} -gt ${TOTAL} ]] && { echo "Sin items asignados"; exit 0; }

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
    echo "[${N}/${N_BATCH}] Procesando: ${SAMPLE}"
    echo "  R1: ${R1}"
    echo "  R2: ${R2}"
    echo "  Hora: $(date +%H:%M:%S)"

    if [[ -s "${DEDUP}" && -s "${DEDUP}.bai" ]]; then
        echo "  Ya mapeado, se omite."
        OK=$(( OK + 1 ))
        continue
    fi
    mkdir -p "${SAMPLE_TMP}"
    RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA\tPU:${SAMPLE}"

    # bwa-mem2 sale con código 0 aunque un FASTQ esté vacío o truncado:
    # se cuentan las lecturas antes y se comparan con las del BAM después
    N_R1=$(( $(zcat "${R1}" | wc -l) / 4 ))
    N_R2=$(( $(zcat "${R2}" | wc -l) / 4 ))
    echo "  Lecturas: R1=${N_R1} R2=${N_R2}"
    if [[ ${N_R1} -eq 0 || ${N_R1} -ne ${N_R2} ]]; then
        echo "  ERROR: ${SAMPLE}: R1 y R2 vacíos o con distinto número de lecturas" >&2
        FAILED=$(( FAILED + 1 ))
        continue
    fi

    # set -e no actúa dentro de un 'if': cada paso se encadena con && para cortar al primer fallo
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
             || { echo "  ERROR: el BAM tiene ${N_PRIMARY} lecturas primarias y se esperaban $(( N_R1 + N_R2 ))" >&2; false; }; } \
        && echo "  [5/5] limpiando intermedios..." \
        && rm -f "${SORTED}" "${SORTED}.bai" \
        && rm -rf "${SAMPLE_TMP:?}"; then
        echo "  Mapeadas: $(awk '/primary mapped/ {print $1, $6}' "${PREFIX}_flagstat.txt" | tr -d '(')"
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
echo "Resumen tarea ${SLURM_ARRAY_TASK_ID}:"
echo "  Procesados: ${OK}"
echo "  Fallidos:   ${FAILED}"
echo "  Fin: $(date)"
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
# Lanzamiento
# ---------------------------------------------------------------------------
PREV_DEP=""
JOB_LINES=()
N_SUBMITTED=0
if [[ "${START_STEP}" -le 1 ]]; then
    OUTPUT=$(sbatch ${PREV_DEP} "${SCRIPT1}")
    JOB1_ID=$(echo "${OUTPUT}" | awk '{print $NF}')
    PREV_DEP="--dependency=afterok:${JOB1_ID}"
    JOB_LINES+=(" Job 1 (${JOB1_ID}): índices de la referencia")
    N_SUBMITTED=$(( N_SUBMITTED + 1 ))
fi
if [[ "${START_STEP}" -le 2 ]]; then
    OUTPUT=$(sbatch ${PREV_DEP} "${SCRIPT2}")
    JOB2_ID=$(echo "${OUTPUT}" | awk '{print $NF}')
    JOB_LINES+=(" Job 2 (${JOB2_ID}): bwa-mem2 + MarkDuplicates, array 0-$(( NUM_JOBS - 1 ))${PREV_DEP:+ (${PREV_DEP#--dependency=})}")
    N_SUBMITTED=$(( N_SUBMITTED + NUM_JOBS ))
fi

echo "============================================="
echo " Resumen"
echo "============================================="
echo " Referencia:   ${REF}"
for DIR in "${INPUT_DIRS[@]}"; do
echo " Entrada:      $(realpath "${DIR}")"
done
echo " Muestras:     ${TOTAL} (${LIST_FILE})"
echo " Salida BAM:   ${BAM_DIR}"
echo " Hilos/job:    ${THREADS}"
echo " Paso inicial: ${START_STEP}"
echo "---------------------------------------------"
printf '%s\n' "${JOB_LINES[@]}"
echo "---------------------------------------------"
echo " Total jobs: ${N_SUBMITTED} / ${MAX_JOBS}"
echo "============================================="
echo
echo "Comandos útiles:"
echo "  squeue -u \$(whoami)"
echo "  tail -F ${OUTPUT_DIR}/logs/bwa_markdup_*.out"
echo "  grep -l 'ERROR' ${OUTPUT_DIR}/logs/*.err"
echo "  grep -H 'primary mapped' ${BAM_DIR}/*_flagstat.txt"
