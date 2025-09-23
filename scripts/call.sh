#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# getting root directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# sourcing helper and mode-specific variant calling scripts
source "$SCRIPT_DIR/scripts/helper.sh"
source "$SCRIPT_DIR/scripts/tumor_only_calling_workflow.sh"
source "$SCRIPT_DIR/scripts/tumor_normal_calling_workflow.sh"
source "$SCRIPT_DIR/scripts/germline_calling_workflow.sh"

# function to define call specific command usage
usage() {
  cat << EOF
Usage: ${TOOL_NAME} call [arguments] [options]

Arguments: 
        --tumor-bam <path>          preprocessed tumor BAM file [REQUIRED]
        --normal-bam <path>         preprocessed normal BAM file [REQUIRED if --mode tumor-normal]
    -r, --reference <path>          reference FASTA file [REQUIRED]
    -m, --mode <string>             mode to run variant calling. possible values: {tumor-only, tumor-normal, germline} [REQUIRED]
    -i, --intervals <path>          BED/VCF/.interval_list/.list/.intervals file specifying exon capture intervals [REQUIRED]
    -t, --threads <integer>         threads to use in programs that support multithreading [OPTIONAL] [1]    

Options: 
    -v, --verbose                   display tool outputs
    -h, --help                      show help message                        
EOF
}

call_main() {
# initializing variables
VERBOSE=0
TUMOR_BAM=""
NORMAL_BAM=""
MODE=""
THREADS=1
REFERENCE=""
INTERVAL_LIST=""

# parse arguments
OPTS=$(getopt -o t:i:m:r:vh --long tumor-bam:,normal-bam:,threads:,intervals:,mode:,reference:,verbose,help -n "$0" -- "$@") || exit 1
eval set -- "$OPTS"

# setting variables based on inputted arguments
while true; do
  case "$1" in
    --tumor-bam) TUMOR_BAM="$2"; shift 2 ;;
    --normal-bam) NORMAL_BAM="$2"; shift 2 ;;
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
if [[ "$MODE" != "tumor-only" && "$MODE" != "tumor-normal" && "$MODE" != "germline" ]]; then
    echo "[ERROR] must specify either --mode tumor-only, --mode tumor-normal or --mode germline"
    exit 1
fi

# --reference
if [[ -z "$REFERENCE" ]]; then 
  echo "[ERROR] --reference is required"
  exit 1
fi

if [[ ! -f "$REFERENCE" ]]; then
  echo "[ERROR] Reference file does not exist: $REFERENCE"
  exit 1
fi

# --intervals
if [[ -z "$INTERVAL_LIST" ]]; then 
  echo "[ERROR] --intervals is required"
  exit 1
fi

if [[ ! -f "$INTERVAL_LIST" ]]; then
  echo "[ERROR] Interval file does not exist: $INTERVAL_LIST"
  exit 1
fi

# --tumor-bam
if [[ -z "$TUMOR_BAM" ]]; then
  echo "[ERROR] --tumor-bam is required"
  exit 1
fi

if [[ ! -f "$TUMOR_BAM" ]]; then
  echo "[ERROR] Tumor BAM file does not exist: $TUMOR_BAM"
  exit 1
fi

# --normal-bam for tumor-normal mode
if [[ "$MODE" == "tumor-normal" ]]; then
  if [[ -z "$NORMAL_BAM" ]]; then 
    echo "[ERROR] --normal-bam required in tumor-normal mode"
    exit 1
  fi
  if [[ ! -f "$NORMAL_BAM" ]]; then
    echo "[ERROR] Normal BAM file does not exist: $NORMAL_BAM"
    exit 1
  fi
fi

# setting start 
SECONDS=0

# running either tumor-only, tumor-normal, or germline variant calling by calling functions from sourced scripts
if [[ "$MODE" == "tumor-only" ]]; then
    tumor_only_variant_calling_workflow
elif [[ "$MODE" == "tumor-normal" ]]; then
    tumor_normal_variant_calling_workflow
elif [[ "$MODE" == "germline" ]]; then
    germline_variant_calling_workflow
fi

# getting elapsed time
elapsed=$SECONDS
hours=$(( elapsed / 3600 ))
minutes=$(( (elapsed % 3600) / 60 ))
seconds=$(( elapsed % 60 ))

echo "Total runtime: ${hours}h ${minutes}m ${seconds}s"

}