#!/bin/bash
# =============================================================================
# run_variant_calling.sh — Variant calling conjunto con GATK (HPC Drago, SLURM)
#
# Paso 1: HaplotypeCaller -ploidy 1 -ERC BP_RESOLUTION por muestra (array)
# Paso 2: CombineGVCFs + GenotypeGVCFs -ploidy 1 (job único)
# Paso 3: filtros duros (SNP / no-SNP), MergeVcfs y VariantsToTable (job único)
#
# Salida: gvcf/SAMPLE.g.vcf.gz, vcf/cohort.raw.vcf.gz, vcf/cohort.filtered.vcf.gz,
#         vcf/cohort.filtered.PASS.vcf.gz, vcf/cohort.filtered.table.tsv
#
# El genotipado conjunto da genotipo a todas las muestras en cada sitio variante:
# 0 = referencia con cobertura, ./. = sin cobertura (nunca se asume referencia).
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Valores por defecto
# ---------------------------------------------------------------------------
REF=""
INPUT_DIRS=()
OUTPUT_DIR=""
MAX_JOBS=60
THREADS=4
START_STEP=1

usage() {
    cat << USAGE
Uso: $(basename "$0") -r REF.fa -i DIR_BAM [-i DIR2 ...] -o OUTPUT_DIR [-j MAX_JOBS] [-t THREADS] [-s PASO] [-h]

Llama variantes en BAM deduplicados (run_bwa_gatk.sh) o recalibrados (run_bqsr.sh) con HaplotypeCaller
en modo GVCF por base (BP_RESOLUTION, haploide), genotipa todas las muestras juntas y aplica los filtros duros
recomendados por GATK. Los sitios que no pasan se marcan en FILTER, no se eliminan.
Todo el trabajo se lanza a SLURM; este script solo genera y envía los jobs.

Flags:
  -r  Referencia FASTA (obligatorio). Debe tener .fai y .dict (los crea run_bwa_gatk.sh)
  -i  Directorio con *_recal.bam o *_sorted_dedup.bam + .bai (obligatorio, repetible).
      Si hay de los dos tipos en un directorio, se usan los *_recal.bam
  -o  Directorio de salida (obligatorio). GVCF en OUTPUT_DIR/gvcf, VCF en OUTPUT_DIR/vcf,
      logs en OUTPUT_DIR/logs
  -j  Límite total de jobs en cola (default: ${MAX_JOBS} = MaxJobsPU)
  -t  Hilos por job de HaplotypeCaller (default: ${THREADS})
  -s  Paso desde el que empezar (default: 1)
        1 = HaplotypeCaller por muestra
        2 = CombineGVCFs + GenotypeGVCFs
        3 = filtros + tabla
  -h  Muestra esta ayuda

Ejemplos:
  conda activate gatk
  $(basename "$0") -r /lustre/home/iata/aaguilar/paper/ref/nuclear.fasta \\
      -i /lustre/home/iata/aaguilar/paper/mapping/bam -o /lustre/home/iata/aaguilar/paper/calling
  $(basename "$0") -r ref.fa -i bam -o calling -s 3      # solo rehacer filtros

Las muestras con GVCF e índice ya presentes se omiten al relanzar.
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
if ! [[ "${MAX_JOBS}" =~ ^[0-9]+$ && "${THREADS}" =~ ^[0-9]+$ && "${START_STEP}" =~ ^[123]$ ]]; then
    echo "ERROR: -j y -t deben ser enteros y -s debe ser 1, 2 o 3."
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

# ---------------------------------------------------------------------------
# Directorios
# ---------------------------------------------------------------------------
REF=$(realpath "${REF}")
OUTPUT_DIR=$(realpath -m "${OUTPUT_DIR}")
GVCF_DIR="${OUTPUT_DIR}/gvcf"
VCF_DIR="${OUTPUT_DIR}/vcf"
mkdir -p "${GVCF_DIR}" "${VCF_DIR}" "${OUTPUT_DIR}/logs" "${OUTPUT_DIR}/tmp"

for f in "${REF}.fai" "${REF%.*}.dict"; do
    [[ -f "${f}" ]] || { echo "ERROR: falta ${f}. Créalo con el paso 1 de run_bwa_gatk.sh."; exit 1; }
done

# ---------------------------------------------------------------------------
# Lista de muestras: SAMPLE<TAB>BAM
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
            echo "AVISO: sin índice para ${BAM}, se omite."
            continue
        fi
        if cut -f1 "${LIST_FILE}" | grep -qxF "${SAMPLE}"; then
            echo "ERROR: muestra duplicada entre directorios: ${SAMPLE}"
            exit 1
        fi
        printf '%s\t%s\n' "${SAMPLE}" "${BAM}" >> "${LIST_FILE}"
    done <<< "${BAMS}"
done

TOTAL=$(wc -l < "${LIST_FILE}")
[[ ${TOTAL} -gt 0 ]] || { echo "ERROR: no se encontraron *_recal.bam ni *_sorted_dedup.bam en ${INPUT_DIRS[*]}"; exit 1; }

if [[ "${START_STEP}" -eq 2 ]]; then
    while IFS=$'\t' read -r SAMPLE BAM; do
        [[ -s "${GVCF_DIR}/${SAMPLE}.g.vcf.gz.tbi" ]] \
            || { echo "ERROR: falta ${GVCF_DIR}/${SAMPLE}.g.vcf.gz(.tbi). Lanza desde el paso 1."; exit 1; }
    done < "${LIST_FILE}"
fi
if [[ "${START_STEP}" -eq 3 && ! -s "${VCF_DIR}/cohort.raw.vcf.gz.tbi" ]]; then
    echo "ERROR: falta ${VCF_DIR}/cohort.raw.vcf.gz(.tbi). Lanza desde el paso 2."
    exit 1
fi

# Los pasos 2 y 3 son jobs secuenciales: se restan del pool de arrays
NUM_SEQ_JOBS=$(( 4 - START_STEP ))
[[ ${NUM_SEQ_JOBS} -gt 2 ]] && NUM_SEQ_JOBS=2
NUM_JOBS=$(( MAX_JOBS - NUM_SEQ_JOBS ))
[[ ${NUM_JOBS} -gt ${TOTAL} ]] && NUM_JOBS=${TOTAL}
[[ ${NUM_JOBS} -lt 1 ]] && { echo "ERROR: -j demasiado bajo."; exit 1; }

# ---------------------------------------------------------------------------
# Script SLURM — Paso 1: HaplotypeCaller por muestra (array)
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

echo "=== HaplotypeCaller (GVCF por base, ploidía 1) ==="
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
while IFS=$'\t' read -r SAMPLE BAM; do
    N=$(( N + 1 ))
    GVCF="${GVCF_DIR}/${SAMPLE}.g.vcf.gz"
    SAMPLE_TMP="${TMP_DIR}/${SAMPLE}"
    echo
    echo "[${N}/${N_BATCH}] Procesando: ${SAMPLE}"
    echo "  BAM: ${BAM}"
    echo "  Hora: $(date +%H:%M:%S)"

    if [[ -s "${GVCF}" && -s "${GVCF}.tbi" ]]; then
        echo "  Ya llamado, se omite."
        OK=$(( OK + 1 ))
        continue
    fi
    mkdir -p "${SAMPLE_TMP}"

    # BP_RESOLUTION: una línea por base, sin bloques de referencia. Así el DP de las
    # muestras sin la variante es el de esa posición y no el mínimo de un bloque
    # set -e no actúa dentro de un 'if': cada paso se encadena con && para cortar al primer fallo
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
echo "Resumen tarea ${SLURM_ARRAY_TASK_ID}:"
echo "  Procesados: ${OK}"
echo "  Fallidos:   ${FAILED}"
echo "  Fin: $(date)"
echo "========================"
[[ ${FAILED} -eq 0 ]]
EOF

# ---------------------------------------------------------------------------
# Script SLURM — Paso 2: CombineGVCFs + GenotypeGVCFs (job único)
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
echo "Referencia: ${REF}"
echo "Muestras: $(wc -l < "${LIST_FILE}")"
echo "Inicio: $(date)"
echo "========================"

VARIANT_ARGS=()
while IFS=$'\t' read -r SAMPLE BAM; do
    GVCF="${GVCF_DIR}/${SAMPLE}.g.vcf.gz"
    [[ -s "${GVCF}.tbi" ]] || { echo "ERROR: falta ${GVCF}(.tbi)" >&2; exit 1; }
    VARIANT_ARGS+=(-V "${GVCF}")
done < "${LIST_FILE}"

echo "  [1/2] CombineGVCFs..."
gatk --java-options "-Xmx12g" CombineGVCFs \
    -R "${REF}" \
    "${VARIANT_ARGS[@]}" \
    -O "${COMBINED}" \
    --tmp-dir "${TMP_DIR}"

echo "  [2/2] GenotypeGVCFs (ploidía 1)..."
gatk --java-options "-Xmx12g" GenotypeGVCFs \
    -R "${REF}" \
    -V "${COMBINED}" \
    -O "${RAW}" \
    -ploidy 1 \
    --tmp-dir "${TMP_DIR}"

echo "  Variantes en bruto: $(zgrep -vc '^#' "${RAW}")"

echo
echo "========================"
echo "  Fin: $(date)"
echo "========================"
EOF

# ---------------------------------------------------------------------------
# Script SLURM — Paso 3: filtros duros + tabla (job único)
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

echo "=== Filtros duros GATK ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "VCF: ${RAW}"
echo "Inicio: $(date)"
echo "========================"

# Se separa SNP / resto (INDEL + MIXED) para no perder los sitios mixtos
echo "  [1/5] SelectVariants SNP / no-SNP..."
gatk SelectVariants -R "${REF}" -V "${RAW}" --select-type-to-include SNP -O "${P}.snps.raw.vcf.gz"
gatk SelectVariants -R "${REF}" -V "${RAW}" --select-type-to-exclude SNP -O "${P}.indels.raw.vcf.gz"

# Umbrales recomendados por GATK para filtros duros. LowDP marca el genotipo (FT), no el sitio
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
echo "  Sitios totales: $(zgrep -vc '^#' "${P}.filtered.vcf.gz")"
echo "  Sitios PASS:    $(zgrep -vc '^#' "${P}.filtered.PASS.vcf.gz")"
echo "========================"
echo "  Fin: $(date)"
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
# Lanzamiento
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
JOB_LINES+=(" Job 3 (${JOB3_ID}): filtros + tabla${PREV_DEP:+ (${PREV_DEP#--dependency=})}")
N_SUBMITTED=$(( N_SUBMITTED + 1 ))

echo "============================================="
echo " Resumen"
echo "============================================="
echo " Referencia:   ${REF}"
for DIR in "${INPUT_DIRS[@]}"; do
echo " Entrada:      $(realpath "${DIR}")"
done
echo " Muestras:     ${TOTAL} (${LIST_FILE})"
echo " Salida GVCF:  ${GVCF_DIR}"
echo " Salida VCF:   ${VCF_DIR}"
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
echo "  tail -F ${OUTPUT_DIR}/logs/hc_gvcf_*.out"
echo "  grep -h 'OK:\|ERROR' ${OUTPUT_DIR}/logs/hc_gvcf_*.out"
echo "  grep -l 'ERROR' ${OUTPUT_DIR}/logs/*.err"
echo "  tail ${OUTPUT_DIR}/logs/filter_vcf_*.out"
