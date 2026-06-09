#!/bin/bash
#===============================================================================
# run_sppIDer_array.sh — Lanzador de sppIDer como array job en Drago
#
# Reparte N muestras entre M jobs (default: 60 = MaxJobsPU).
# Cada job procesa su lote secuencialmente.
# Las muestras con _run2/_run3/... se fusionan automáticamente.
# Detecta automáticamente muestras single-end (SE) y paired-end (PE).
#
# Uso:
#   bash run_sppIDer_array.sh
#   bash run_sppIDer_array.sh -i /path/reads -r /path/ref.fasta -o /path/out
#   bash run_sppIDer_array.sh -s _R1.trimmed.fastq.gz -S _R2.trimmed.fastq.gz
#===============================================================================

set -euo pipefail

# ── Valores por defecto ──
READS_DIR="/lustre/home/iata/aaguilar/somics/seub/fastp_out/trimmed"
REF="/lustre/home/iata/aaguilar/somics/ref_genomes/comb/comb.fasta"
OUTDIR="/lustre/home/iata/aaguilar/somics/seub/sppider_out"
SUFFIX1=""
SUFFIX2=""
NUM_JOBS=60
EXTRA_ARGS=""
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Sufijos candidatos R1→R2, en orden de preferencia
CANDIDATE_S1=(
    "_R1.trimmed.fastq.gz"
    "_1.trimmed.fastq.gz"
    "_R1.fastq.gz"
    "_1.fastq.gz"
)
CANDIDATE_S2=(
    "_R2.trimmed.fastq.gz"
    "_2.trimmed.fastq.gz"
    "_R2.fastq.gz"
    "_2.fastq.gz"
)

# Sufijos candidatos single-end, en orden de preferencia
CANDIDATE_SE=(
    ".trimmed.fastq.gz"
    ".fastq.gz"
)

# ── Parseo de argumentos ──
while getopts "i:r:o:j:s:S:bkh" opt; do
    case ${opt} in
        i) READS_DIR="${OPTARG}" ;;
        r) REF="${OPTARG}" ;;
        o) OUTDIR="${OPTARG}" ;;
        j) NUM_JOBS="${OPTARG}" ;;
        s) SUFFIX1="${OPTARG}" ;;
        S) SUFFIX2="${OPTARG}" ;;
        b) EXTRA_ARGS="${EXTRA_ARGS} --byGroup" ;;
        k) EXTRA_ARGS="${EXTRA_ARGS} --keep-intermediates" ;;
        h)
            echo "Uso: bash $0 [-i reads_dir] [-r ref.fasta] [-o outdir] [-j max_jobs]"
            echo "            [-s SUFFIX_R1] [-S SUFFIX_R2] [-b] [-k]"
            echo ""
            echo "Sufijos PE candidatos autodetectados (en orden de prioridad):"
            for i in "${!CANDIDATE_S1[@]}"; do
                echo "  $((i+1)). ${CANDIDATE_S1[$i]} / ${CANDIDATE_S2[$i]}"
            done
            echo ""
            echo "Sufijos SE candidatos autodetectados:"
            for s in "${CANDIDATE_SE[@]}"; do
                echo "  - ${s}"
            done
            exit 0 ;;
        *) exit 1 ;;
    esac
done

READS_DIR=$(realpath "${READS_DIR}")
REF=$(realpath "${REF}")

[[ ! -d "$READS_DIR" ]] && { echo "ERROR: reads-dir no existe: $READS_DIR"; exit 1; }
[[ ! -f "$REF" ]]       && { echo "ERROR: ref no existe: $REF"; exit 1; }

# ── Autodetección de sufijos PE ──
if [[ -z "${SUFFIX1}" ]]; then
    echo "Sufijo PE no especificado, autodetectando en ${READS_DIR}..."
    for i in "${!CANDIDATE_S1[@]}"; do
        S1="${CANDIDATE_S1[$i]}"
        S2="${CANDIDATE_S2[$i]}"
        COUNT=$(find "${READS_DIR}" -name "*${S1}" | wc -l)
        if [[ "${COUNT}" -gt 0 ]]; then
            SUFFIX1="${S1}"
            SUFFIX2="${S2}"
            echo "  → Sufijo PE detectado: R1='${SUFFIX1}' / R2='${SUFFIX2}' (${COUNT} archivos R1)"
            break
        fi
    done
    if [[ -z "${SUFFIX1}" ]]; then
        echo "AVISO: No se detectó ningún sufijo R1 PE conocido. Solo se procesarán muestras SE."
    fi
else
    if [[ -z "${SUFFIX2}" ]]; then
        SUFFIX2="${SUFFIX1/R1/R2}"
        SUFFIX2="${SUFFIX2/_1./_2.}"
        echo "AVISO: -S no especificado, infiriendo SUFFIX2='${SUFFIX2}'"
    fi
    echo "Sufijos PE manuales: R1='${SUFFIX1}' / R2='${SUFFIX2}'"
fi

mkdir -p "${OUTDIR}/logs"

# ── Generar lista de muestras ──
# Formato TSV: SAMPLE_BASE <TAB> MODE(PE|SE) <TAB> R1_1:R1_2:... <TAB> R2_1:R2_2:... (R2 vacío para SE)
SAMPLES_FILE="${OUTDIR}/sppIDer_samples.tsv"
> "${SAMPLES_FILE}"

declare -A R1_MAP
declare -A R2_MAP
declare -A MODE_MAP

# ── Recolectar muestras PE ──
if [[ -n "${SUFFIX1}" ]]; then
    while IFS= read -r r1; do
        raw_name="$(basename "${r1}" "${SUFFIX1}")"
        r1_dir="$(dirname "${r1}")"
        r2="${r1_dir}/${raw_name}${SUFFIX2}"

        if [[ ! -f "$r2" ]]; then
            echo "AVISO: R2 no encontrado para '${raw_name}' (esperado: ${r2}), tratando como SE" >&2
            # Tratar como SE si no hay R2
            sample_base=$(echo "${raw_name}" | sed 's/_run[0-9]\+$//')
            if [[ -z "${R1_MAP[$sample_base]+x}" ]]; then
                R1_MAP["$sample_base"]="${r1}"
                R2_MAP["$sample_base"]=""
                MODE_MAP["$sample_base"]="SE"
            else
                R1_MAP["$sample_base"]="${R1_MAP[$sample_base]}:${r1}"
                MODE_MAP["$sample_base"]="SE"
            fi
            continue
        fi

        sample_base=$(echo "${raw_name}" | sed 's/_run[0-9]\+$//')

        if [[ -z "${R1_MAP[$sample_base]+x}" ]]; then
            R1_MAP["$sample_base"]="${r1}"
            R2_MAP["$sample_base"]="${r2}"
            MODE_MAP["$sample_base"]="PE"
        else
            R1_MAP["$sample_base"]="${R1_MAP[$sample_base]}:${r1}"
            R2_MAP["$sample_base"]="${R2_MAP[$sample_base]}:${r2}"
            # Mantener PE si ya era PE
        fi

    done < <(find "${READS_DIR}" -name "*${SUFFIX1}" | sort)
fi

# ── Recolectar muestras SE (archivos que no tienen sufijo R1/R2) ──
# Buscar archivos que coincidan con sufijos SE pero que NO sean R1/R2 ya registrados
for SE_SUFFIX in "${CANDIDATE_SE[@]}"; do
    # Saltar si el sufijo SE es más genérico y ya está cubierto por el sufijo PE
    # (evitar capturar _R1.trimmed.fastq.gz con .trimmed.fastq.gz)
    while IFS= read -r se_file; do
        fname="$(basename "${se_file}")"

        # Excluir archivos que ya son R1 o R2 de PE
        IS_PE=false
        if [[ -n "${SUFFIX1}" ]] && [[ "${fname}" == *"${SUFFIX1}" ]]; then
            IS_PE=true
        fi
        if [[ -n "${SUFFIX2}" ]] && [[ "${fname}" == *"${SUFFIX2}" ]]; then
            IS_PE=true
        fi
        # Excluir también patrones genéricos de R1/R2
        if [[ "${fname}" =~ _R[12]\. ]] || [[ "${fname}" =~ _[12]\. ]]; then
            IS_PE=true
        fi

        [[ "${IS_PE}" == "true" ]] && continue

        raw_name="$(basename "${se_file}" "${SE_SUFFIX}")"
        sample_base=$(echo "${raw_name}" | sed 's/_run[0-9]\+$//')

        # No sobrescribir si ya fue registrado como PE
        if [[ -n "${R1_MAP[$sample_base]+x}" ]]; then
            continue
        fi

        if [[ -z "${R1_MAP[$sample_base]+x}" ]]; then
            R1_MAP["$sample_base"]="${se_file}"
            R2_MAP["$sample_base"]=""
            MODE_MAP["$sample_base"]="SE"
        else
            R1_MAP["$sample_base"]="${R1_MAP[$sample_base]}:${se_file}"
            MODE_MAP["$sample_base"]="SE"
        fi

    done < <(find "${READS_DIR}" -name "*${SE_SUFFIX}" | sort)

    # Solo usar el primer sufijo SE que encuentre archivos nuevos
    # (para no duplicar con sufijos más específicos)
    break
done

# ── Volcar al TSV ordenado ──
for sample_base in $(echo "${!R1_MAP[@]}" | tr ' ' '\n' | sort); do
    echo -e "${sample_base}\t${MODE_MAP[$sample_base]}\t${R1_MAP[$sample_base]}\t${R2_MAP[$sample_base]}" >> "${SAMPLES_FILE}"
done

N_SAMPLES=$(wc -l < "${SAMPLES_FILE}")

if [[ "$N_SAMPLES" -eq 0 ]]; then
    echo "ERROR: No se encontraron reads en ${READS_DIR}"
    exit 1
fi

# ── Informe de muestras detectadas ──
echo ""
echo "Muestras detectadas:"
N_PE=0
N_SE=0
while IFS=$'\t' read -r sample mode r1s r2s; do
    N_RUNS=$(echo "${r1s}" | tr ':' '\n' | wc -l)
    if [[ "${mode}" == "PE" ]]; then
        N_PE=$((N_PE + 1))
        if [[ "${N_RUNS}" -gt 1 ]]; then
            echo "  [PE][${N_RUNS} runs] ${sample}"
        else
            echo "  [PE][1 run ] ${sample}"
        fi
    else
        N_SE=$((N_SE + 1))
        if [[ "${N_RUNS}" -gt 1 ]]; then
            echo "  [SE][${N_RUNS} runs] ${sample}"
        else
            echo "  [SE][1 run ] ${sample}"
        fi
    fi
done < "${SAMPLES_FILE}"
echo ""
echo "  Total: ${N_SAMPLES} muestras (${N_PE} PE, ${N_SE} SE)"
echo ""

# ── Ajustar NUM_JOBS ──
if [[ "${N_SAMPLES}" -lt "${NUM_JOBS}" ]]; then
    NUM_JOBS="${N_SAMPLES}"
    echo "Ajustando a ${NUM_JOBS} jobs (menos muestras que jobs solicitados)"
fi

SAMPLES_PER_JOB=$(( (N_SAMPLES + NUM_JOBS - 1) / NUM_JOBS ))
ARRAY_MAX=$(( NUM_JOBS - 1 ))

echo "Muestras: ${N_SAMPLES} | Jobs: ${NUM_JOBS} | ~${SAMPLES_PER_JOB} muestras/job"

# ── Crear el script SLURM worker ──
SLURM_SCRIPT="${OUTDIR}/sppIDer_worker.slurm"

# Parte 1: heredoc sin expansión (variables SLURM)
cat > "${SLURM_SCRIPT}" << 'SLURM_EOF'
#!/bin/bash
#SBATCH --job-name=sppIDer
#SBATCH --partition=generic
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=OUTDIR_PLACEHOLDER/logs/sppIDer_%A_%a.out
#SBATCH --error=OUTDIR_PLACEHOLDER/logs/sppIDer_%A_%a.err

set -euo pipefail

module load rama0.4
module load GCCcore/13.3.0
module load BWA/0.7.18
module load GCC/13.3.0
module load SAMtools/1.21
module load BEDTools/2.31.1

export PS1=""
source /dragofs/sw/foss/0.2/software/Miniconda3/4.9.2/etc/profile.d/conda.sh
conda activate new_sppider

SLURM_EOF

# Parte 2: heredoc con expansión (variables del launcher)
cat >> "${SLURM_SCRIPT}" << SLURM_EOF2
# ── Configuración (expandida en tiempo de submit) ──
SAMPLES_FILE="${SAMPLES_FILE}"
REF="${REF}"
OUTDIR="${OUTDIR}"
SCRIPT_DIR="${SCRIPT_DIR}"
EXTRA_ARGS="${EXTRA_ARGS}"
N_SAMPLES=${N_SAMPLES}
NUM_JOBS=${NUM_JOBS}
SLURM_EOF2

# Parte 3: heredoc sin expansión (lógica del worker)
cat >> "${SLURM_SCRIPT}" << 'SLURM_EOF3'

# ── Calcular rango de muestras para este task (0-based) ──
SAMPLES_PER_JOB=$(( (N_SAMPLES + NUM_JOBS - 1) / NUM_JOBS ))
START_LINE=$(( SLURM_ARRAY_TASK_ID * SAMPLES_PER_JOB + 1 ))
END_LINE=$(( START_LINE + SAMPLES_PER_JOB - 1 ))

[[ ${END_LINE} -gt ${N_SAMPLES} ]] && END_LINE=${N_SAMPLES}

BATCH=$(sed -n "${START_LINE},${END_LINE}p" "${SAMPLES_FILE}")
BATCH_SIZE=$(echo "${BATCH}" | grep -c . || true)

echo "════════════════════════════════════════════════════════════"
echo " Task      : ${SLURM_ARRAY_TASK_ID}"
echo " Job       : ${SLURM_JOB_ID}"
echo " Nodo      : $(hostname)"
echo " Muestras  : ${START_LINE}-${END_LINE} (${BATCH_SIZE} muestras)"
echo " Inicio    : $(date)"
echo "════════════════════════════════════════════════════════════"

PROCESSED=0
FAILED=0

while IFS=$'\t' read -r SAMPLE MODE R1S R2S; do
    [[ -z "${SAMPLE}" ]] && continue

    PROCESSED=$((PROCESSED + 1))
    SAMPLE_OUTDIR="${OUTDIR}/${SAMPLE}"
    mkdir -p "${SAMPLE_OUTDIR}"

    # Convertir listas ":" en arrays
    IFS=':' read -ra R1_LIST <<< "${R1S}"
    N_RUNS=${#R1_LIST[@]}

    echo ""
    echo "[${PROCESSED}/${BATCH_SIZE}] Procesando: ${SAMPLE} (${MODE}, ${N_RUNS} run/s)"
    for ((i=0; i<N_RUNS; i++)); do
        echo "  R1[run$((i+1))]: ${R1_LIST[$i]}"
    done
    if [[ "${MODE}" == "PE" ]]; then
        IFS=':' read -ra R2_LIST <<< "${R2S}"
        for ((i=0; i<N_RUNS; i++)); do
            echo "  R2[run$((i+1))]: ${R2_LIST[$i]}"
        done
    fi
    echo "  Hora: $(date '+%H:%M:%S')"

    # ── Fusionar runs si hay más de uno ──
    if [[ ${N_RUNS} -eq 1 ]]; then
        FINAL_R1="${R1_LIST[0]}"
        [[ "${MODE}" == "PE" ]] && FINAL_R2="${R2_LIST[0]}" || FINAL_R2=""
        MERGED=false
    else
        TMP_R1="${SAMPLE_OUTDIR}/${SAMPLE}_merged_R1.fastq.gz"
        echo "  Fusionando ${N_RUNS} runs → $(basename ${TMP_R1})"
        cat "${R1_LIST[@]}" > "${TMP_R1}"
        FINAL_R1="${TMP_R1}"
        MERGED=true

        if [[ "${MODE}" == "PE" ]]; then
            TMP_R2="${SAMPLE_OUTDIR}/${SAMPLE}_merged_R2.fastq.gz"
            echo "  Fusionando ${N_RUNS} runs → $(basename ${TMP_R2})"
            cat "${R2_LIST[@]}" > "${TMP_R2}"
            FINAL_R2="${TMP_R2}"
        else
            FINAL_R2=""
        fi
    fi

    # ── Llamar a sppIDer (PE o SE) ──
    if [[ "${MODE}" == "PE" ]]; then
        python3 "${SCRIPT_DIR}/sppIDer.py" \
            --sample  "${SAMPLE}" \
            --ref     "${REF}" \
            --r1      "${FINAL_R1}" \
            --r2      "${FINAL_R2}" \
            --outdir  "${SAMPLE_OUTDIR}" \
            --threads "${SLURM_CPUS_PER_TASK}" \
            ${EXTRA_ARGS}
    else
        python3 "${SCRIPT_DIR}/sppIDer.py" \
            --sample  "${SAMPLE}" \
            --ref     "${REF}" \
            --r1      "${FINAL_R1}" \
            --outdir  "${SAMPLE_OUTDIR}" \
            --threads "${SLURM_CPUS_PER_TASK}" \
            ${EXTRA_ARGS}
    fi

    EXIT_CODE=$?

    # ── Limpiar temporales de fusión ──
    if [[ "${MERGED}" == "true" ]] && [[ "${EXTRA_ARGS}" != *"--keep-intermediates"* ]]; then
        rm -f "${SAMPLE_OUTDIR}/${SAMPLE}_merged_R1.fastq.gz"
        rm -f "${SAMPLE_OUTDIR}/${SAMPLE}_merged_R2.fastq.gz"
        echo "  Temporales eliminados."
    fi

    if [[ ${EXIT_CODE} -ne 0 ]]; then
        echo "  ERROR: sppIDer falló para ${SAMPLE}"
        FAILED=$((FAILED + 1))
    else
        echo "  OK: ${SAMPLE} [${MODE}]"
    fi

done <<< "${BATCH}"

echo ""
echo "════════════════════════════════════════════════════════════"
echo " Resumen task ${SLURM_ARRAY_TASK_ID}:"
echo "   Procesados: ${PROCESSED}"
echo "   Fallidos:   ${FAILED}"
echo "   Fin: $(date)"
echo "════════════════════════════════════════════════════════════"

[[ ${FAILED} -gt 0 ]] && exit 1
exit 0
SLURM_EOF3

# Sustituir placeholder del OUTDIR en las directivas #SBATCH
sed -i "s|OUTDIR_PLACEHOLDER|${OUTDIR}|g" "${SLURM_SCRIPT}"

# ── Lanzar el array job ──
echo "Lanzando array job..."
JOB_OUTPUT=$(sbatch --array="0-${ARRAY_MAX}" "${SLURM_SCRIPT}")
echo "${JOB_OUTPUT}"
JOB_ID=$(echo "${JOB_OUTPUT}" | awk '{print $NF}')

echo ""
echo "============================================="
echo " Resumen"
echo "============================================="
echo " Reads dir:       ${READS_DIR}"
echo " Ref:             ${REF}"
echo " Outdir:          ${OUTDIR}"
echo " Sufijo R1:       ${SUFFIX1:-N/A}"
echo " Sufijo R2:       ${SUFFIX2:-N/A}"
echo " Total muestras:  ${N_SAMPLES} (${N_PE} PE, ${N_SE} SE)"
echo " Jobs:            ${NUM_JOBS} (array 0-${ARRAY_MAX})"
echo " Muestras/job:    ~${SAMPLES_PER_JOB}"
echo " Extra args:      ${EXTRA_ARGS:-ninguno}"
echo " Samples file:    ${SAMPLES_FILE}"
echo " Worker script:   ${SLURM_SCRIPT}"
echo " Logs:            ${OUTDIR}/logs/"
echo "============================================="
echo ""
echo "Comandos útiles:"
echo "  squeue -j ${JOB_ID}"
echo "  scancel ${JOB_ID}"
echo "  tail -f ${OUTDIR}/logs/sppIDer_${JOB_ID}_0.out"
echo "  grep -l 'ERROR' ${OUTDIR}/logs/sppIDer_${JOB_ID}_*.err"
