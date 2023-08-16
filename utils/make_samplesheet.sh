#!/bin/bash
MODE="single"
while getopts ":p" FLAG; do
	case $FLAG in
		p)
			echo "Paired-end mode"
			MODE="paired"
			;;
	esac
done
if [[ ! -f "samplesheet.csv" ]]; then
	touch samplesheet.csv
fi
echo "samples" >> samplesheet.csv
if [[ $MODE == "single" ]]; then
	echo "Making single-end samplesheet"
	for i in $(ls "./00.fastq/"); do
		i="${i##*/}"
		i="${i%%.*}"
		echo $i >> samplesheet.csv
	done
elif [[ $MODE == "paired" ]]; then
	echo "Making paired-end samplesheet"
	for i in $(ls "./00.fastq/"); do
		i="${i##*/}"
		if [[ $i =~ "R2" ]]; then
			continue
		fi
		i="${i%%_R1*}"
		echo $i >> samplesheet.csv
	done
fi
