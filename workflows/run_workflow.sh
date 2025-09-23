#!/bin/bash
# Execute a WDL workflow with Cromwell

USAGE="USAGE:  $0 -j|--jar <cromwell jar, if not in PATH> -c|--config <cromwell config> -i <input JSON> -o <options file> -p <zipped imports> -w <wdl workflow file>"

MODE="single"
STYLE="R"

if [[ $# -lt 5 ]]; then
	echo ${USAGE}
	exit 1
fi

while [[ $# -gt 0 ]]; do
	case "$1" in
		-j | --jar ) 
			JAR=$2
			shift 2
			;;
		-c | --config ) 
			CONF="$2"
			shift 2
			;;
		-i )
			INPUT="$2"
			shift 2
			;;
		-o ) 
			OPT="$2"
			shift 2
			;;
		-p ) 
			IMPORTS="$2"
			shift 2
			;;
		-w )
			WDL="$2"
			shift 2
			;;
		-* | --* )
			echo ${USAGE}
			exit 1
			;;
	esac
done

# Run cromwell command
java -Dconfig.file=${CONF} -jar ${JAR} run -i ${INPUT} -o ${OPT} -p ${IMPORTS} ${WDL}
