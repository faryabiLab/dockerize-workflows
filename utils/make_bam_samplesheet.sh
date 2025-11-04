#!/bin/bash
# Noah Burget — BAM samplesheet generator
# Generates samplesheet.tsv with <sample_name>\t<path_to_bam>

USAGE="USAGE: $0 -d|--directory <bam directory>"

if [[ $# -lt 2 ]]; then
    echo "${USAGE}"
    exit 1
fi

# Parse args
while [[ $# -gt 0 ]]; do
    case "$1" in
        -d | --directory )
            dir="$2"
            shift 2
            ;;
        -* | --* )
            echo "${USAGE}"
            exit 1
            ;;
    esac
done

# Validate directory
if [[ ! -d "${dir}" ]]; then
    echo "Error: directory '${dir}' not found."
    exit 1
fi

# Remove existing samplesheet if present
if [[ -f "samplesheet.tsv" ]]; then
    rm samplesheet.tsv
fi

touch samplesheet.tsv

echo "Generating samplesheet from BAM files in ${dir}..."

# Loop through all BAM files
for bam in "${dir}"/*.bam; do
    # Skip if no BAMs found
    [[ -e "$bam" ]] || { echo "No BAM files found in ${dir}"; exit 1; }

    sample=$(basename "$bam")
    sample=${sample%%.*}   # remove extension
    echo -e "${sample}\t$(realpath "$bam")" >> samplesheet.tsv
done

echo "Done. Wrote $(wc -l < samplesheet.tsv) entries to samplesheet.tsv"

