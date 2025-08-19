#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# getting root directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# sourcing helper and mode-specific pre-processing scripts
source "$SCRIPT_DIR/scripts/helper.sh"
source "$SCRIPT_DIR/scripts/tumor_only_preprocessing.sh"
source "$SCRIPT_DIR/scripts/tumor_normal_preprocessing.sh"

# function to define preprocess specific command usage
usage() {
  cat << EOF
Usage: ${TOOL_NAME} preprocess [arguments] [options]

Arguments: 
    -s, --sample-name <string>      sample identifier [REQUIRED]
        --tumor-fq1 <path>          tumor forward FASTQ file [REQUIRED]
        --tumor-fq2 <path>          tumor reverse FASTQ file [REQUIRED]
        --normal-fq1 <path>         normal forward FASTQ file [REQUIRED if --mode tumor-normal]
        --normal-fq2 <path>         normal reverse FASTQ file [REQUIRED if --mode tumor-normal]
    -r, --reference <path>          reference hg38/GRCh38 FASTA file [REQUIRED]
    -m, --mode <string>             mode to run preprocessing. possible values: {tumor-only, tumor-normal} [REQUIRED]
    -i, --intervals <path>          BED/VCF/.interval_list/.list/.intervals file specifying exon capture intervals [REQUIRED]
    -t, --threads <integer>         threads to use in programs that support multithreading [OPTIONAL] [1]    

Options: 
-v, --verbose                       display tool outputs
-h, --help                          show help message                        
EOF
}

preprocess_main() {
# initializing variables
VERBOSE=0
SAMPLE_NAME=""
MODE=""
TUMOR_FQ1=""
TUMOR_FQ2=""
NORMAL_FQ1=""
NORMAL_FQ2=""
THREADS=1
REFERENCE=""
BQSR_DB_DIR="./bqsr_training"
BQSR_OUT_DIR="./bqsr_output"
MAPPED_READS_DIR="./mapped_reads"
INTERVAL_LIST=""

# parse arguments
OPTS=$(getopt -o s:t:i:m:r:vh --long sample-name:,tumor-fq1:,tumor-fq2:,normal-fq1:,normal-fq2:,threads:,intervals:,mode:,reference:,verbose,help -n "$0" -- "$@") || exit 1
eval set -- "$OPTS"

# setting variables based on inputted arguments
while true; do
  case "$1" in
    -s | --sample-name) SAMPLE_NAME="$2"; shift 2 ;;
    --tumor-fq1) TUMOR_FQ1="$2"; shift 2 ;;
    --tumor-fq2) TUMOR_FQ2="$2"; shift 2 ;;
    --normal-fq1) NORMAL_FQ1="$2"; shift 2 ;;
    --normal-fq2) NORMAL_FQ2="$2"; shift 2 ;;
    -r | --reference) REFERENCE="$2"; shift 2 ;;
    -m | --mode) MODE="$2"; shift 2 ;;
    -t | --threads) THREADS="$2"; shift 2 ;;
    -i | --intervals) INTERVAL_LIST="$2"; shift 2 ;;
    -v | --verbose) VERBOSE=1; shift ;;
    -h | --help) usage; exit 0 ;;
    --) shift; break ;;
    *) echo "Invalid argument"; exit 1 ;;
  esac
done

#### Creating errors if required inputs are not present ####

# --mode
if [[ "$MODE" != "tumor-only" && "$MODE" != "tumor-normal" ]]; then
    echo "[ERROR] must specify either --mode tumor-only or --mode tumor-normal"
    exit 1
fi

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

# --tumor-fq1 and --tumor-fq2
if [[ -z "$TUMOR_FQ1" || -z "$TUMOR_FQ2" ]]; then
  echo "[ERROR] --tumor-fq1 and --tumor-fq2 are required"
  exit 1
fi

# --normal-fq1 and --normal-fq2 for tumor-normal mode
if [[ "$MODE" == "tumor-normal" ]]; then
  if [[ -z "$NORMAL_FQ1" || -z "$NORMAL_FQ2" ]]; then 
    echo "[ERROR] --normal-fq1 and --normal-fq2 are required in tumor-normal mode"
    exit 1
  fi
fi

# setting start 
SECONDS=0

# running either tumor-only or tumor-normal preprocessing steps by calling functions from sources scripts
if [[ "$MODE" == "tumor-only" ]]; then
    tumor_only_preprocess_workflow
  else
    tumor_normal_preprocess_workflow
fi

# getting elapsed time
elapsed=$SECONDS
hours=$(( elapsed / 3600 ))
minutes=$(( (elapsed % 3600) / 60 ))
seconds=$(( elapsed % 60 ))

echo "Total runtime: ${hours}h ${minutes}m ${seconds}s"

}