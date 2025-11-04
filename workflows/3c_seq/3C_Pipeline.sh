#!/usr/bin/env bash
#
# juicer_wrapper.sh — Run Juicer and convert output to .mcool format
#
# Example:
#   bash juicer_wrapper.sh \
#     -j /mnt/data0/noah/software/juicer/CPU/juicer.sh \
#     -d /mnt/data0/noah/processing/032524_hic/large_Pre_B_MboI \
#     -s MboI \
#     -p /mnt/data0/noah/references/Mus_musculus/GRCm38/Annotation/ChromSizes_ordered.tsv \
#     -z /mnt/data0/noah/references/Mus_musculus/GRCm38/Sequence/BWAIndex/genome.fa \
#     -y /mnt/data0/noah/references/Mus_musculus/GRCm38/GRCm38_MboI.txt \
#     -D /mnt/data0/noah/software/juicer \
#     -t 12 \
#     -c myCoolFile
#
# Author: Noah Burget

set -euo pipefail

########################################
# Default values
########################################
THREADS=8
COOL_PREFIX="output"

########################################
# Help message
########################################
usage() {
    echo "Usage: $0 -j JUICER -d TOP_DIR -s RESTRICTION_ENZYME -p CHROMSIZES -z BWA_GENOME -y RESTRICTION_SITES -D JUICER_DIR [-t THREADS] [-c COOL_PREFIX]"
    echo
    echo "Required arguments:"
    echo "  -j  Path to juicer.sh executable"
    echo "  -d  Top-level output directory (Juicer working dir)"
    echo "  -s  Restriction enzyme (e.g. MboI, HindIII)"
    echo "  -p  Chromosome sizes file"
    echo "  -z  BWA genome FASTA"
    echo "  -y  Restriction sites file"
    echo "  -D  Juicer installation directory"
    echo
    echo "Optional arguments:"
    echo "  -t  Number of threads to use (default: 8)"
    echo "  -c  Prefix for generated .mcool file (default: output)"
    echo "  -h  Show this help message and exit"
    echo
    echo "Example:"
    echo "  bash $0 -j /path/to/juicer.sh -d /data/exp -s MboI -p chrom.sizes -z genome.fa -y MboI.txt -D /software/juicer -t 12 -c myCoolFile"
    exit 1
}

########################################
# Parse arguments
########################################
while getopts ":j:d:s:p:z:y:D:t:c:h" opt; do
    case "${opt}" in
        j) JUICER=${OPTARG} ;;
        d) TOP_DIR=${OPTARG} ;;
        s) RESTRICTION_ENZYME=${OPTARG} ;;
        p) CHROMSIZES=${OPTARG} ;;
        z) BWA_GENOME=${OPTARG} ;;
        y) RESTRICTION_SITES=${OPTARG} ;;
        D) JUICER_DIR=${OPTARG} ;;
        t) THREADS=${OPTARG} ;;
        c) COOL_PREFIX=${OPTARG} ;;
        h) usage ;;
        *) echo "Invalid option: -${OPTARG}" >&2; usage ;;
    esac
done

########################################
# Check required args
########################################
for var in JUICER TOP_DIR RESTRICTION_ENZYME CHROMSIZES BWA_GENOME RESTRICTION_SITES JUICER_DIR; do
    if [[ -z "${!var:-}" ]]; then
        echo "Error: Missing required argument: $var" >&2
        usage
    fi
done

########################################
# Run Juicer
########################################
echo "[INFO] Running Juicer..."
echo "Command:"
echo "  bash ${JUICER} -d ${TOP_DIR} -s ${RESTRICTION_ENZYME} -p ${CHROMSIZES} -z ${BWA_GENOME} -y ${RESTRICTION_SITES} -D ${JUICER_DIR} -t ${THREADS}"
echo

bash "${JUICER}" \
  -d "${TOP_DIR}" \
  -s "${RESTRICTION_ENZYME}" \
  -p "${CHROMSIZES}" \
  -z "${BWA_GENOME}" \
  -y "${RESTRICTION_SITES}" \
  -D "${JUICER_DIR}" \
  -t "${THREADS}"

########################################
# Convert to .mcool
########################################
HIC_FILE="${TOP_DIR}/aligned/inter_30.hic"
OUT_FILE="${COOL_PREFIX}.mcool"

if [[ -f "${HIC_FILE}" ]]; then
    echo "[INFO] Converting ${HIC_FILE} -> ${OUT_FILE}"
    hic2cool convert -r 0 "${HIC_FILE}" "${OUT_FILE}"
    echo "[INFO] Conversion complete: ${OUT_FILE}"
else
    echo "[ERROR] Expected Hi-C file not found at ${HIC_FILE}" >&2
    exit 1
fi

echo "[INFO] Done."

