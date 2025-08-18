#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# getting root directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# sourcing helper and mode-specific variant calling scripts
source "$SCRIPT_DIR/scripts/helper.sh"
source "$SCRIPT_DIR/scripts/tumor_only_calling.sh"
source "$SCRIPT_DIR/scripts/tumor_normal_calling.sh"

# function to define preprocess specific command usage
usage() {
  cat << EOF
Usage: ${TOOL_NAME} call [arguments] [options]

Arguments: 
    -s, --sample-name <string>      sample identifier [REQUIRED]
        --tumor-bam <path>          sorted, de-duplicated, recalibrated tumor BAM file [REQUIRED]
        --normal-bam <path>         sorted, de-duplicated, recalibrated normal BAM file [REQUIRED]
    -r, --reference <path>          reference FASTA file [REQUIRED]
    -m, --mode <string>             mode to run preprocessing. possible values: {tumor-only, tumor-normal} [REQUIRED]
    -i, --intervals <path>          BED/VCF/.interval_list/.list/.intervals file specifying exon capture intervals [REQUIRED]
    -t, --threads <integer>         threads to used in programs that support multithreading [OPTIONAL] [1]    

Options: 
-a --keep-all                       retain all variants, including those that fail filters (not recommended, may increase annotation time)
-v, --verbose                       display tool outputs
-h, --help                          show help message                        
EOF
}

call_main() {
# initializing variables
VERBOSE=0
SAMPLE_NAME=""
TUMOR_BAM=""
NORMAL_BAM=""
KEEP_ALL=0
MODE=""
GENOME_VERS="hg38"
THREADS=1
REFERENCE=""
MUTECT2_SUPPORTING_FILES_DIR=./mutect2_supporting_files
MUTECT2_OUTPUT_DIR=./mutect2_output
MUTECT2_FILTERING_DIR=./mutect2_filtering_data
MAPPED_READS_DIR="./mapped_reads"
INTERVAL_LIST=""

# parse arguments
OPTS=$(getopt -o s:t:i:m:r:avh --long sample-name:,tumor-bam:,normal-bam:,threads:,intervals:,mode:,reference:,keep-all,verbose,help -n "$0" -- "$@") || exit 1
eval set -- "$OPTS"

# setting variables based on inputted arguments
while true; do
  case "$1" in
    -s | --sample-name) SAMPLE_NAME="$2"; shift 2 ;;
    --tumor-bam) TUMOR_BAM="$2"; shift 2 ;;
    --normal-bam) NORMAL_BAM="$2"; shift 2 ;;
    -r | --reference) REFERENCE="$2"; shift 2 ;;
    -m | --mode) MODE="$2"; shift 2 ;;
    -t | --threads) THREADS="$2"; shift 2 ;;
    -i | --intervals) INTERVAL_LIST="$2"; shift 2 ;;
    -a | --keep-all) KEEP_ALL=1; shift;;
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

# --tumor-bam
if [[ -z "$TUMOR_BAM" ]]; then
  echo "[ERROR] --tumor-bam is required"
  exit 1
fi

# --normal-bam for tumor-normal mode
if [[ "$MODE" == "tumor-normal" ]]; then
  if [[ -z "$NORMAL_BAM" ]]; then 
    echo "[ERROR] --normal-bam required in tumor-normal mode"
    exit 1
  fi
fi

# setting start 
SECONDS=0

# running either tumor-only or tumor-normal variant calling by calling functions from sources scripts
if [[ "$MODE" == "tumor-only" ]]; then
    tumor_only_variant_calling_workflow
  else
    tumor_normal_variant_calling_workflow
fi

# getting elapsed time
elapsed=$SECONDS
hours=$(( elapsed / 3600 ))
minutes=$(( (elapsed % 3600) / 60 ))
seconds=$(( elapsed % 60 ))

echo "Total runtime: ${hours}h ${minutes}m ${seconds}s"

}