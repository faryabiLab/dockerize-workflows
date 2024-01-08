###############
# This script will take in an arbitrary amount of BED files and return the Union of all regions

while getopts ":p:n:" opt; do
  case $opt in
    p) peaks="$OPTARG";;
    n) name="$OPTARG";;
    \?) echo "Invalid option: -$OPTARG, Please use as such: <>.sh -p <path to dir of peak files>" >&2; exit 1;;
  esac
done

FILES=$(ls $peaks)

cat ${peaks}/* | cut -f 1-3 - | bedtools sort -i stdin | bedtools merge -i stdin -sorted > ${name}_consensusPeaks.bed
