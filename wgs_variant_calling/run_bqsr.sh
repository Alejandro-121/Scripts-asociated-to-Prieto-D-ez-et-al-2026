#!/bin/bash
# =============================================================================
# run_bqsr.sh — Recalibrado de calidades de base (BQSR) con GATK (HPC Drago, SLURM)
#
# Paso único (array), por muestra:
#   BaseRecalibrator (antes) → ApplyBQSR → BaseRecalibrator (después) → AnalyzeCovariates
#
# Sin catálogo de variantes conocidas para la cepa: se usan como sitios conocidos las
# variantes PASS de una primera llamada (bootstrapping, cohort.filtered.PASS.vcf.gz
# de run_variant_calling.sh sobre los BAM sin recalibrar).
#
# Salida por muestra: SAMPLE_recal.bam(.bai), SAMPLE_bqsr_before.table,
#                     SAMPLE_bqsr_after.table, SAMPLE_bqsr_covariates.csv (y .pdf si hay R)
# Después: run_variant_calling.sh -i OUTPUT_DIR/bam para la llamada definitiva.
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Valores por defecto
# ---------------------------------------------------------------------------
REF=""
KNOWN=""
INPUT_DIRS=()
OUTPUT_DIR=""
MAX_JOBS=60
THREADS=2

usage() {
    cat << USAGE
Uso: $(basename "$0") -r REF.fa -k KNOWN.vcf.gz -i DIR_BAM [-i DIR2 ...] -o OUTPUT_DIR [-j MAX_JOBS] [-t THREADS] [-h]

Recalibra las calidades de base de BAM deduplicados (salida de run_bwa_gatk.sh) con
GATK BaseRecalibrator + ApplyBQSR, usando como sitios conocidos las variantes de alta
confianza de una primera llamada. Evalúa el recalibrado con una segunda pasada de
BaseRecalibrator y AnalyzeCovariates.
Todo el trabajo se lanza a SLURM; este script solo genera y envía los jobs.

Flags:
  -r  Referencia FASTA (obligatorio). Debe tener .fai y .dict (los crea run_bwa_gatk.sh)
  -k  VCF de sitios conocidos, bgzip + .tbi (obligatorio). Normalmente
        <salida de run_variant_calling.sh>/vcf/cohort.filtered.PASS.vcf.gz
  -i  Directorio con *_sorted_dedup.bam + .bai (obligatorio, repetible)
  -o  Directorio de salida (obligatorio). BAM en OUTPUT_DIR/bam, tablas en OUTPUT_DIR/tables,
      logs en OUTPUT_DIR/logs
  -j  Límite total de jobs en cola (default: ${MAX_JOBS} = MaxJobsPU)
  -t  Hilos por job (default: ${THREADS})
  -h  Muestra esta ayuda

Ejemplos:
  conda activate gatk
  $(basename "$0") -r /lustre/home/iata/aaguilar/paper/ref/nuclear.fasta \\
      -k /lustre/home/iata/aaguilar/paper/vc_bp/vcf/cohort.filtered.PASS.vcf.gz \\
      -i /lustre/home/iata/aaguilar/paper/mapped/bam -o /lustre/home/iata/aaguilar/paper/bqsr

Las muestras con BAM recalibrado e índice ya presentes se omiten al relanzar.
USAGE
}

# ---------------------------------------------------------------------------
# Argumentos
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
    echo "ERROR: faltan argumentos obligatorios (-r, -k, -i, -o)."
    usage
    exit 1
fi
if ! [[ "${MAX_JOBS}" =~ ^[0-9]+$ && "${THREADS}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: -j y -t deben ser enteros."
    exit 1
fi
[[ -f "${REF}" ]] || { echo "ERROR: no existe la referencia ${REF}"; exit 1; }
[[ -f "${KNOWN}" ]] || { echo "ERROR: no existe el VCF de sitios conocidos ${KNOWN}"; exit 1; }
[[ -f "${KNOWN}.tbi" ]] || { echo "ERROR: falta el índice ${KNOWN}.tbi"; exit 1; }
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
command -v samtools > /dev/null || { echo "ERROR: 'samtools' no está en el entorno 'gatk'."; exit 1; }

# ---------------------------------------------------------------------------
# Directorios
# ---------------------------------------------------------------------------
REF=$(realpath "${REF}")
KNOWN=$(realpath "${KNOWN}")
OUTPUT_DIR=$(realpath -m "${OUTPUT_DIR}")
BAM_DIR="${OUTPUT_DIR}/bam"
TABLE_DIR="${OUTPUT_DIR}/tables"
mkdir -p "${BAM_DIR}" "${TABLE_DIR}" "${OUTPUT_DIR}/logs" "${OUTPUT_DIR}/tmp"

for f in "${REF}.fai" "${REF%.*}.dict"; do
    [[ -f "${f}" ]] || { echo "ERROR: falta ${f}. Créalo con el paso 1 de run_bwa_gatk.sh."; exit 1; }
done
N_KNOWN=$(zgrep -vc '^#' "${KNOWN}" || true)
[[ ${N_KNOWN} -gt 0 ]] || { echo "ERROR: ${KNOWN} no tiene variantes."; exit 1; }

# ---------------------------------------------------------------------------
# Lista de muestras: SAMPLE<TAB>BAM
# ---------------------------------------------------------------------------
LIST_FILE="${OUTPUT_DIR}/samples.tsv"
: > "${LIST_FILE}"
for DIR in "${INPUT_DIRS[@]}"; do
    DIR=$(realpath "${DIR}")
    while read -r BAM; do
        [[ -z "${BAM}" ]] && continue
        SAMPLE=$(basename "${BAM}" _sorted_dedup.bam)
        if [[ ! -f "${BAM}.bai" ]]; then
            echo "AVISO: sin índice para ${BAM}, se omite."
            continue
        fi
        if cut -f1 "${LIST_FILE}" | grep -qxF "${SAMPLE}"; then
            echo "ERROR: muestra duplicada entre directorios: ${SAMPLE}"
            exit 1
        fi
        printf '%s\t%s\n' "${SAMPLE}" "${BAM}" >> "${LIST_FILE}"
    done <<< "$(find "${DIR}" -maxdepth 1 -type f -name "*_sorted_dedup.bam" | sort)"
done

TOTAL=$(wc -l < "${LIST_FILE}")
[[ ${TOTAL} -gt 0 ]] || { echo "ERROR: no se encontraron *_sorted_dedup.bam en ${INPUT_DIRS[*]}"; exit 1; }

NUM_JOBS=${MAX_JOBS}
[[ ${NUM_JOBS} -gt ${TOTAL} ]] && NUM_JOBS=${TOTAL}
[[ ${NUM_JOBS} -lt 1 ]] && { echo "ERROR: -j demasiado bajo."; exit 1; }

# ---------------------------------------------------------------------------
# Script SLURM — BQSR por muestra (array)
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

echo "=== BQSR (sitios conocidos: primera llamada) ==="
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID}, Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Referencia: ${REF}"
echo "Sitios conocidos: ${KNOWN}"
echo "Inicio: $(date)"
echo "========================"

ITEMS_PER_JOB=$(( (TOTAL + NUM_JOBS - 1) / NUM_JOBS ))
START_LINE=$(( SLURM_ARRAY_TASK_ID * ITEMS_PER_JOB + 1 ))
END_LINE=$(( START_LINE + ITEMS_PER_JOB - 1 ))
[[ ${END_LINE} -gt ${TOTAL} ]] && END_LINE=${TOTAL}
[[ ${START_LINE} -gt ${TOTAL} ]] && { echo "Sin items asignados"; exit 0; }

# AnalyzeCovariates solo genera el PDF si hay R con ggplot2/gplots/gsalib; el CSV siempre
PLOTS_OK=0
Rscript -e 'for (p in c("ggplot2","gplots","gsalib")) stopifnot(requireNamespace(p, quietly=TRUE))' \
    > /dev/null 2>&1 && PLOTS_OK=1
[[ ${PLOTS_OK} -eq 1 ]] || echo "AVISO: sin R/ggplot2/gsalib, AnalyzeCovariates solo generará el CSV."

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
    echo "[${N}/${N_BATCH}] Procesando: ${SAMPLE}"
    echo "  BAM: ${BAM}"
    echo "  Hora: $(date +%H:%M:%S)"

    if [[ -s "${RECAL}" && -s "${RECAL}.bai" && -s "${AFTER}" ]]; then
        echo "  Ya recalibrado, se omite."
        OK=$(( OK + 1 ))
        continue
    fi
    mkdir -p "${SAMPLE_TMP}"
    PLOT_ARGS=()
    [[ ${PLOTS_OK} -eq 1 ]] && PLOT_ARGS=(-plots "${COVAR}.pdf")

    # set -e no actúa dentro de un 'if': cada paso se encadena con && para cortar al primer fallo
    if echo "  [1/5] BaseRecalibrator (antes)..." \
        && gatk --java-options "-Xmx8g" BaseRecalibrator \
            -R "${REF}" -I "${BAM}" --known-sites "${KNOWN}" \
            -O "${BEFORE}" --tmp-dir "${SAMPLE_TMP}" \
        && echo "  [2/5] ApplyBQSR..." \
        && gatk --java-options "-Xmx8g" ApplyBQSR \
            -R "${REF}" -I "${BAM}" --bqsr-recal-file "${BEFORE}" \
            -O "${RECAL}" --create-output-bam-index false --tmp-dir "${SAMPLE_TMP}" \
        && echo "  [3/5] samtools index..." \
        && samtools index "${RECAL}" \
        && echo "  [4/5] BaseRecalibrator (después)..." \
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
    -e "s|__KNOWN__|${KNOWN}|g" \
    -e "s|__LIST_FILE__|${LIST_FILE}|g" \
    -e "s|__BAM_DIR__|${BAM_DIR}|g" \
    -e "s|__TABLE_DIR__|${TABLE_DIR}|g" \
    -e "s|__NUM_JOBS__|${NUM_JOBS}|g" \
    "${SCRIPT1}"

# ---------------------------------------------------------------------------
# Lanzamiento
# ---------------------------------------------------------------------------
OUTPUT=$(sbatch "${SCRIPT1}")
JOB1_ID=$(echo "${OUTPUT}" | awk '{print $NF}')

echo "============================================="
echo " Resumen"
echo "============================================="
echo " Referencia:   ${REF}"
echo " Conocidos:    ${KNOWN} (${N_KNOWN} variantes)"
for DIR in "${INPUT_DIRS[@]}"; do
echo " Entrada:      $(realpath "${DIR}")"
done
echo " Muestras:     ${TOTAL} (${LIST_FILE})"
echo " Salida BAM:   ${BAM_DIR}"
echo " Tablas:       ${TABLE_DIR}"
echo " Hilos/job:    ${THREADS}"
echo "---------------------------------------------"
echo " Job 1 (${JOB1_ID}): BQSR, array 0-$(( NUM_JOBS - 1 ))"
echo "---------------------------------------------"
echo " Total jobs: ${NUM_JOBS} / ${MAX_JOBS}"
echo "============================================="
echo
echo "Comandos útiles:"
echo "  squeue -u \$(whoami)"
echo "  tail -F ${OUTPUT_DIR}/logs/bqsr_*.out"
echo "  grep -h 'OK:\|ERROR' ${OUTPUT_DIR}/logs/bqsr_*.out"
echo "  grep -l 'ERROR' ${OUTPUT_DIR}/logs/*.err"
echo
echo "Siguiente paso (llamada definitiva sobre los BAM recalibrados):"
echo "  $(dirname "$0")/run_variant_calling.sh -r ${REF} -i ${BAM_DIR} -o <SALIDA>"
