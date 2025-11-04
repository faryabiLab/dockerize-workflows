#!/usr/bin/env bash
set -euo pipefail

USAGE="Usage: $0 -j <cromwell jar> -c <cromwell config> -h <http://ip:port> -p <imports.zip> -o <options.json> -i <inputs.json> -w <workflow.wdl>"

while getopts ":j:c:h:p:o:i:w:" opt; do
  case $opt in
    j) JAR="$OPTARG";;
    c) CONFIG="$OPTARG";;
    h) HOST="$OPTARG";;
    p) IMPORTS="$OPTARG";;
    o) OPTIONS="$OPTARG";;
    i) INPUTS="$OPTARG";;
    w) WORKFLOW="$OPTARG";;
    *) echo $USAGE; exit 1;;
  esac
done

if [[ -z "${JAR:-}" || -z "${CONFIG:-}" || -z "${HOST:-}" || -z "${WORKFLOW:-}" ]]; then
  echo -e "Missing required arguments.\n$USAGE"; exit 1
fi

if ! curl -s --max-time 5 "$HOST/api/workflows/v1/version" | grep -q "cromwell"; then
  echo "Cromwell server not responding at $HOST"; exit 1
fi

java -Dconfig.file="$CONFIG" -jar "$JAR" submit -h "$HOST" -w "$WORKFLOW" \
  ${IMPORTS:+-p "$IMPORTS"} ${OPTIONS:+-o "$OPTIONS"} ${INPUTS:+-i "$INPUTS"}

