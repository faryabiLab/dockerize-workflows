#!/bin/bash
while true; do
        case "$1" in
                -d | --directory ) dir="$2"; shift ;;
                * ) break ;;
        esac
done

if [[ ! -d $dir ]]; then
	echo "Directory $dir not found, exiting."
	exit 1
fi

for file in $dir/*_1*; do
	mv $file "${file/_1/_R1}"
done

for file in $dir/*_2*; do
        mv $file "${file/_2/_R2}"
done
