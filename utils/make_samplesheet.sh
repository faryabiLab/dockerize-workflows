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
		echo "$j\t$(realpath ${i})" >> samplesheet.tsv
	done
elif [[ $MODE == "paired" ]]; then
	echo "Making paired-end samplesheet"
	for i in $(ls ${dir}); do
		i="${i##*/}"
		echo $i
		if [[ $i =~ "_${STYLE}2" ]]; then
			R2="${i}"
			continue
		fi
		R1="${i}"
		i=${i%%_${STYLE}1*}
		echo "$i\t$(realpath ${R1})\t$(realpath ${R2})" >> samplesheet.tsv
	done
fi
