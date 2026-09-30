#!/bin/bash
# =============================================================================
# annotate_snpeff.sh — Anotación funcional con SnpEff (base de datos propia R64-1-1)
#
# Paso 1: construye la base de datos SnpEff "R64-1-1" a partir de la referencia usada
#         en el mapeo y del GFF de SGD R64-1-1 (2011), comprobando CDS y proteínas
#         contra orf_coding_all / orf_trans_all de SGD
# Paso 2: anota el VCF conjunto (cohort.filtered.vcf.gz de run_variant_calling.sh)
#
# Salida: VCF.ann.vcf.gz(.tbi), snpeff_stats.csv, snpeff_summary.html, snpeff.log
# Entorno: conda create -n snpeff -c conda-forge -c bioconda snpeff bcftools htslib
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Valores por defecto
# ---------------------------------------------------------------------------
REF=""
SGD_DIR=""
VCF=""
SNPEFF_DIR=""
GENOME="R64-1-1"
START_STEP=1

usage() {
    cat << USAGE
Uso: $(basename "$0") -r REF.fa -g DIR_SGD -v VCF.gz -d DIR_SNPEFF [-s PASO] [-h]

Construye una base de datos SnpEff para S. cerevisiae S288C R64-1-1 y anota un VCF.
Los cromosomas de la referencia deben llamarse chrI..chrXVI, como en el GFF de SGD.

Flags:
  -r  Referencia FASTA nuclear usada en el mapeo (obligatorio)
  -g  Directorio S288C_reference_genome_R64-1-1_20110203 de SGD (obligatorio en el paso 1)
  -v  VCF bgzip a anotar (obligatorio en el paso 2)
  -d  Directorio de SnpEff: snpEff.config y data/${GENOME}/ (obligatorio)
  -s  Paso desde el que empezar (default: 1)
        1 = construir la base de datos
        2 = anotar el VCF
  -h  Muestra esta ayuda

Ejemplo:
  conda activate snpeff
  $(basename "$0") -r reference/nuclear.fasta \\
      -g reference/S288C_reference_genome_R64-1-1_20110203 \\
      -v vc_final/vcf/cohort.filtered.vcf.gz -d snpeff
USAGE
}

# ---------------------------------------------------------------------------
# Argumentos
# ---------------------------------------------------------------------------
while getopts "r:g:v:d:s:h" opt; do
    case ${opt} in
        r) REF=${OPTARG} ;;
        g) SGD_DIR=${OPTARG} ;;
        v) VCF=${OPTARG} ;;
        d) SNPEFF_DIR=${OPTARG} ;;
        s) START_STEP=${OPTARG} ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

if [[ -z "${SNPEFF_DIR}" ]] || ! [[ "${START_STEP}" =~ ^[12]$ ]]; then
    echo "ERROR: falta -d o -s no es 1 o 2."
    usage
    exit 1
fi
if [[ "${START_STEP}" -le 1 ]]; then
    [[ -f "${REF}" ]] || { echo "ERROR: no existe la referencia '${REF}' (-r)"; exit 1; }
    [[ -d "${SGD_DIR}" ]] || { echo "ERROR: no existe el directorio de SGD '${SGD_DIR}' (-g)"; exit 1; }
fi
[[ -f "${VCF}" ]] || { echo "ERROR: no existe el VCF '${VCF}' (-v)"; exit 1; }

# ---------------------------------------------------------------------------
# Entorno
# ---------------------------------------------------------------------------
if [[ "${CONDA_DEFAULT_ENV:-}" != "snpeff" ]]; then
    echo "AVISO: Activa el entorno conda 'snpeff' antes de lanzar."
    echo "  conda activate snpeff"
    exit 1
fi
for tool in snpEff bgzip tabix; do
    command -v "${tool}" > /dev/null || { echo "ERROR: '${tool}' no está en el entorno 'snpeff'."; exit 1; }
done

SNPEFF_DIR=$(realpath -m "${SNPEFF_DIR}")
DB_DIR="${SNPEFF_DIR}/data/${GENOME}"
CONFIG="${SNPEFF_DIR}/snpEff.config"

# ---------------------------------------------------------------------------
# Paso 1: base de datos
# ---------------------------------------------------------------------------
if [[ "${START_STEP}" -le 1 ]]; then
    GFF=$(find "${SGD_DIR}" -maxdepth 1 -name "saccharomyces_cerevisiae_R64-1-1_*.gff" | head -1)
    CDS=$(find "${SGD_DIR}" -maxdepth 1 -name "orf_coding_all_R64-1-1_*.fasta" | head -1)
    PROT=$(find "${SGD_DIR}" -maxdepth 1 -name "orf_trans_all_R64-1-1_*.fasta" | head -1)
    for f in "${GFF}" "${CDS}" "${PROT}"; do
        [[ -n "${f}" && -f "${f}" ]] || { echo "ERROR: falta GFF, orf_coding_all u orf_trans_all en ${SGD_DIR}"; exit 1; }
    done
    mkdir -p "${DB_DIR}"

    echo "[1/4] genes.gff: cromosomas nucleares, sin la sección ##FASTA, Name = nombre estándar..."
    # El GFF de 2011 lleva el nombre sistemático en Name y el estándar en gene=:
    # SnpEff toma Name como nombre del gen, así que se sustituye cuando hay gene=
    awk -F'\t' -v OFS='\t' '
        { sub(/\r$/, "") }
        /^##FASTA/ { exit }
        /^#/ { print; next }
        NF == 9 && $1 ~ /^chr[IVX]+$/ {
            if ($3 != "CDS" && match($9, /(^|;)gene=[^;]+/)) {
                std = substr($9, RSTART, RLENGTH); sub(/^;?gene=/, "", std)
                sub(/(^|;)Name=[^;]+/, (substr($9, 1, 5) == "Name=" ? "" : ";") "Name=" std, $9)
            }
            print
        }' "${GFF}" > "${DB_DIR}/genes.gff"

    echo "[2/4] sequences.fa: referencia del mapeo..."
    cp "${REF}" "${DB_DIR}/sequences.fa"

    # SnpEff crea un tránscrito TRANSCRIPT_<ORF> por gen (el GFF no tiene mRNA)
    echo "[3/4] cds.fa / protein.fa para la comprobación..."
    sed -E 's/^>([^ ]+).*/>TRANSCRIPT_\1/' "${CDS}" > "${DB_DIR}/cds.fa"
    sed -E 's/^>([^ ]+).*/>TRANSCRIPT_\1/' "${PROT}" > "${DB_DIR}/protein.fa"

    cat > "${CONFIG}" << EOF
# SnpEff: base de datos propia para S. cerevisiae S288C R64-1-1 (SGD, 2011)
data.dir = ./data/
${GENOME}.genome : Saccharomyces cerevisiae S288C R64-1-1 (nuclear)
EOF

    echo "[4/4] snpEff build..."
    snpEff build -gff3 -c "${CONFIG}" "${GENOME}" > "${SNPEFF_DIR}/build.log" 2>&1
    grep -E "CDS check:|Protein check:" "${SNPEFF_DIR}/build.log"
    if grep -E "CDS check:|Protein check:" "${SNPEFF_DIR}/build.log" | grep -qv "Errors: 0"; then
        echo "ERROR: la comprobación de CDS/proteínas tiene errores. Revisa ${SNPEFF_DIR}/build.log"
        exit 1
    fi
fi

# ---------------------------------------------------------------------------
# Paso 2: anotación
# ---------------------------------------------------------------------------
[[ -f "${CONFIG}" ]] || { echo "ERROR: falta ${CONFIG}. Lanza desde el paso 1."; exit 1; }
OUT_DIR=$(dirname "$(realpath "${VCF}")")
OUT="${VCF%.vcf.gz}.ann.vcf.gz"

echo "snpEff ann ${GENOME}: ${VCF}"
snpEff ann -c "${CONFIG}" \
    -csvStats "${OUT_DIR}/snpeff_stats.csv" -s "${OUT_DIR}/snpeff_summary.html" \
    "${GENOME}" "${VCF}" 2> "${OUT_DIR}/snpeff.log" | bgzip > "${OUT}"
tabix -f -p vcf "${OUT}"

N_IN=$(zgrep -vc '^#' "${VCF}" || true)
N_ANN=$(zgrep -v '^#' "${OUT}" | grep -c 'ANN=' || true)

echo "============================================="
echo " Resumen"
echo "============================================="
echo " Base de datos: ${DB_DIR}"
echo " VCF entrada:   ${VCF} (${N_IN} sitios)"
echo " VCF anotado:   ${OUT} (${N_ANN} con ANN)"
echo " Estadísticas:  ${OUT_DIR}/snpeff_summary.html"
echo "============================================="
