#!/bin/bash
if [[ -f "samplesheet.tsv" ]]; then
	rm samplesheet.tsv
fi

USAGE="USAGE:  $0 <-p|--paired> -s|--pair-id-style <R|N> -d <fastq directory>"

MODE="single"
STYLE="R"

if [[ $# -lt 2 ]]; then
	echo ${USAGE}
	exit 1
fi

while [[ $# -gt 0 ]]; do
	case "$1" in
		-p | --paired ) 
			MODE="paired"
			shift 
			;;
		-d | --directory ) 
			dir="$2"
			shift 2
			;;
		-s | --pair-id-style )
			STYLE="$2"
			shift 2
			;;
		-* | --* ) 
			echo ${USAGE}
			exit 1
			;;
	esac
done

if [[ ${STYLE} == "N" ]]; then
	STYLE=""
fi

if [[ ! -f "samplesheet.tsv" ]]; then
	touch samplesheet.tsv
fi
#echo "samples" >> samplesheet.csv
if [[ $MODE == "single" ]]; then
	echo "Making single-end samplesheet"
	for i in $(ls ${dir}); do
		j="${i##*/}"
		j="${j%%.*}"
		echo -e "$j\t$(realpath ${i})\n" >> samplesheet.tsv
	done

elif [[ $MODE == "paired" ]]; then
    for R1 in "${dir}"/*_"${STYLE}"1*.fastq.gz; do
        # Derive basename without _R1/1
        sample=$(basename "$R1")
        sample=${sample%%_"${STYLE}"1*}

        # Matching R2
        R2="${dir}/${sample}_${STYLE}2.fastq.gz"

        # Check that R2 exists
        if [[ -f "$R2" ]]; then
            echo -e "${sample}\t$(realpath "$R1")\t$(realpath "$R2")" >> samplesheet.tsv
        else
            echo "Warning: no R2 found for ${sample}" >&2
        fi
    done
fi
