#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# getting root directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# sourcing helper and mode-specific pre-processing scripts
source "$SCRIPT_DIR/scripts/helper.sh"
source "$SCRIPT_DIR/scripts/sample_preprocessing.sh"

# function to define preprocess specific command usage
usage() {
  cat << EOF
Usage: ${TOOL_NAME} preprocess [arguments] [options]

Arguments: 
    -s, --sample-name <string>      sample identifier [REQUIRED]
        --fastq1 <path>             forward FASTQ file [REQUIRED]
        --fastq2 <path>             reverse FASTQ file [REQUIRED]
    -r, --reference <path>          reference hg38/GRCh38 FASTA file [REQUIRED]
    -i, --intervals <path>          BED/VCF/.interval_list/.list/.intervals file specifying exon capture intervals [REQUIRED]
    -t, --threads <integer>         threads to use in programs that support multithreading [OPTIONAL] [1]    

Options: 
    -v, --verbose                   display tool outputs
    -h, --help                      show help message                        
EOF
}

preprocess_main() {
# initializing variables
VERBOSE=0
SAMPLE_NAME=""
FQ1=""
FQ2=""
THREADS=1
REFERENCE=""
INTERVAL_LIST=""

# parse arguments
OPTS=$(getopt -o s:t:i:r:vh --long sample-name:,fastq1:,fastq2:,threads:,intervals:,reference:,verbose,help -n "$0" -- "$@") || exit 1
eval set -- "$OPTS"

# setting variables based on inputted arguments
while true; do
  case "$1" in
    -s | --sample-name) SAMPLE_NAME="$2"; shift 2 ;;
    --fastq1) FQ1="$2"; shift 2 ;;
    --fastq2) FQ2="$2"; shift 2 ;;
    -r | --reference) REFERENCE="$2"; shift 2 ;;
    -t | --threads) THREADS="$2"; shift 2 ;;
    -i | --intervals) INTERVAL_LIST="$2"; shift 2 ;;
    -v | --verbose) VERBOSE=1; shift ;;
    -h | --help) usage; exit 0 ;;
    --) shift; break ;;
    *) echo "Invalid argument"; exit 1 ;;
  esac
done

#### Creating errors if required inputs are not present ####

# --sample-name 
if [[ -z "$SAMPLE_NAME" ]]; then 
  echo "[ERROR] --sample-name is required"
  exit 1
fi

# --reference
if [[ -z "$REFERENCE" ]]; then 
  echo "[ERROR] --reference is required"
  exit 1
fi

# --intervals
if [[ -z "$INTERVAL_LIST" ]]; then 
  echo "[ERROR] --intervals is required"
  exit 1
fi

# --fastq1 and --fastq2
if [[ -z "$FQ1" || -z "$FQ2" ]]; then
  echo "[ERROR] --fastq1 and --fastq2 are required"
  exit 1
fi

# setting start 
SECONDS=0

# running preprocessing
preprocess_workflow

# getting elapsed time
elapsed=$SECONDS
hours=$(( elapsed / 3600 ))
minutes=$(( (elapsed % 3600) / 60 ))
seconds=$(( elapsed % 60 ))

echo "Total runtime: ${hours}h ${minutes}m ${seconds}s"

}