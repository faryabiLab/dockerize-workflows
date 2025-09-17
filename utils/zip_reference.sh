#!/bin/bash
# Zip a reference directory for use in the Dockerized Pipelines

dir=""
USAGE="Usage: $0 -d <index directory>"

function print_usage_exit()
{
	if [[ ! -z $1 ]]; then
		echo -e "Unknown parameter '$1' \n${USAGE}"
		exit 1
	fi
	echo ${USAGE}
	exit 1
}

if [[ $# -eq 0 ]]; then
	echo $USAGE
	exit 1
fi

while [[ $# -gt 0 ]]; do
	case "$1" in
		-d | --dir )
			dir="$2"
			shift 2 # Shift past both argument and value
			;;
		-* | --* )
			print_usage_exit $1
			;;
	esac
done

ref_full=$(realpath "${ref}")

stamp=$(date +"%Y%m%d_%H%M%S")
base_name=$(basename $0)
base="${base_name%.*}"

script_name="${base}_${stamp}.sh"

# Write the command into a new script
cat > "$script_name" <<EOF
#!/bin/bash
tar -czhf index.tar.gz -C ${dir} .
EOF

chmod +x "$script_name"

echo "Wrote executable script: $script_name"

# Run the command immediately
./"$script_name"
