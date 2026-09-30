#!/bin/bash
# =============================================================================
# annotate_snpeff.sh — Functional annotation with SnpEff (custom R64-1-1 database)
#
# Step 1: builds the SnpEff database "R64-1-1" from the reference used for mapping
#         and the SGD R64-1-1 GFF (2011), checking CDS and protein sequences
#         against SGD orf_coding_all / orf_trans_all
# Step 2: annotates the joint VCF (cohort.filtered.vcf.gz from run_variant_calling.sh)
#
# Output: VCF.ann.vcf.gz(.tbi), snpeff_stats.csv, snpeff_summary.html, snpeff.log
# Environment: conda create -n snpeff -c conda-forge -c bioconda snpeff bcftools htslib
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
REF=""
SGD_DIR=""
VCF=""
SNPEFF_DIR=""
GENOME="R64-1-1"
START_STEP=1

usage() {
    cat << USAGE
Usage: $(basename "$0") -r REF.fa -g SGD_DIR -v VCF.gz -d SNPEFF_DIR [-s STEP] [-h]

Builds a SnpEff database for S. cerevisiae S288C R64-1-1 and annotates a VCF.
Reference chromosomes must be named chrI..chrXVI, as in the SGD GFF.

Flags:
  -r  Nuclear reference FASTA used for mapping (required)
  -g  SGD directory S288C_reference_genome_R64-1-1_20110203 (required for step 1)
  -v  bgzipped VCF to annotate (required for step 2)
  -d  SnpEff directory: snpEff.config and data/${GENOME}/ (required)
  -s  Step to start from (default: 1)
        1 = build the database
        2 = annotate the VCF
  -h  Show this help

Example:
  conda activate snpeff
  $(basename "$0") -r reference/nuclear.fasta \\
      -g reference/S288C_reference_genome_R64-1-1_20110203 \\
      -v vc_final/vcf/cohort.filtered.vcf.gz -d snpeff
USAGE
}

# ---------------------------------------------------------------------------
# Arguments
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
    echo "ERROR: -d is missing or -s is not 1 or 2."
    usage
    exit 1
fi
if [[ "${START_STEP}" -le 1 ]]; then
    [[ -f "${REF}" ]] || { echo "ERROR: reference '${REF}' not found (-r)"; exit 1; }
    [[ -d "${SGD_DIR}" ]] || { echo "ERROR: SGD directory '${SGD_DIR}' not found (-g)"; exit 1; }
fi
[[ -f "${VCF}" ]] || { echo "ERROR: VCF '${VCF}' not found (-v)"; exit 1; }

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
if [[ "${CONDA_DEFAULT_ENV:-}" != "snpeff" ]]; then
    echo "WARNING: activate the conda environment 'snpeff' before launching."
    echo "  conda activate snpeff"
    exit 1
fi
for tool in snpEff bgzip tabix; do
    command -v "${tool}" > /dev/null || { echo "ERROR: '${tool}' is not in the 'snpeff' environment."; exit 1; }
done

SNPEFF_DIR=$(realpath -m "${SNPEFF_DIR}")
DB_DIR="${SNPEFF_DIR}/data/${GENOME}"
CONFIG="${SNPEFF_DIR}/snpEff.config"

# ---------------------------------------------------------------------------
# Step 1: database
# ---------------------------------------------------------------------------
if [[ "${START_STEP}" -le 1 ]]; then
    GFF=$(find "${SGD_DIR}" -maxdepth 1 -name "saccharomyces_cerevisiae_R64-1-1_*.gff" | head -1)
    CDS=$(find "${SGD_DIR}" -maxdepth 1 -name "orf_coding_all_R64-1-1_*.fasta" | head -1)
    PROT=$(find "${SGD_DIR}" -maxdepth 1 -name "orf_trans_all_R64-1-1_*.fasta" | head -1)
    for f in "${GFF}" "${CDS}" "${PROT}"; do
        [[ -n "${f}" && -f "${f}" ]] || { echo "ERROR: GFF, orf_coding_all or orf_trans_all missing in ${SGD_DIR}"; exit 1; }
    done
    mkdir -p "${DB_DIR}"

    echo "[1/4] genes.gff: nuclear chromosomes, without the ##FASTA section, Name = standard name..."
    # The 2011 GFF has the systematic name in Name and the standard name in gene=:
    # SnpEff uses Name as the gene name, so it is replaced when gene= is present
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

    echo "[2/4] sequences.fa: mapping reference..."
    cp "${REF}" "${DB_DIR}/sequences.fa"

    # SnpEff creates one transcript TRANSCRIPT_<ORF> per gene (the GFF has no mRNA features)
    echo "[3/4] cds.fa / protein.fa for the sequence check..."
    sed -E 's/^>([^ ]+).*/>TRANSCRIPT_\1/' "${CDS}" > "${DB_DIR}/cds.fa"
    sed -E 's/^>([^ ]+).*/>TRANSCRIPT_\1/' "${PROT}" > "${DB_DIR}/protein.fa"

    cat > "${CONFIG}" << EOF
# SnpEff: custom database for S. cerevisiae S288C R64-1-1 (SGD, 2011)
data.dir = ./data/
${GENOME}.genome : Saccharomyces cerevisiae S288C R64-1-1 (nuclear)
EOF

    echo "[4/4] snpEff build..."
    snpEff build -gff3 -c "${CONFIG}" "${GENOME}" > "${SNPEFF_DIR}/build.log" 2>&1
    grep -E "CDS check:|Protein check:" "${SNPEFF_DIR}/build.log"
    if grep -E "CDS check:|Protein check:" "${SNPEFF_DIR}/build.log" | grep -qv "Errors: 0"; then
        echo "ERROR: the CDS/protein check reported errors. See ${SNPEFF_DIR}/build.log"
        exit 1
    fi
fi

# ---------------------------------------------------------------------------
# Step 2: annotation
# ---------------------------------------------------------------------------
[[ -f "${CONFIG}" ]] || { echo "ERROR: ${CONFIG} not found. Start from step 1."; exit 1; }
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
echo " Summary"
echo "============================================="
echo " Database:      ${DB_DIR}"
echo " Input VCF:     ${VCF} (${N_IN} sites)"
echo " Annotated VCF: ${OUT} (${N_ANN} with ANN)"
echo " Statistics:    ${OUT_DIR}/snpeff_summary.html"
echo "============================================="
