#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# getting root directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# sourcing helper and annotation scripts
source "$SCRIPT_DIR/scripts/helper.sh"
source "$SCRIPT_DIR/scripts/annotation_workflow.sh"

usage(){
    cat << EOF
Usage: ${TOOL_NAME} annotate [arguments] [options]

Arguments: 
    -s, --sample-name <string>      sample identifier [REQUIRED]
        --vcf <path>                input VCF file for annotation [REQUIRED]
    -m, --mode <string>             mode in which variant calling was run. possible values: {tumor-only, tumor-normal, germline} [REQUIRED]
    -c, --cache <path>              directory housing the unzipped VEP cache file [REQUIRED]
    -r, --reference <path>          reference hg38/GRCh38 FASTA file [REQUIRED]
    -a, --pop-af <float>            maximum population allele frequency (0-1) threshold for filtering [OPTIONAL] [0.01]

Options: 
    -v, --verbose                   display tool outputs
    -h, --help                      show help message  
EOF
}

annotate_main(){
# initializing variables
VERBOSE=0
SAMPLE_NAME=""
REFERENCE=""
CACHE=""
MODE=""
POP_AF=0.01
VCF=""

# parse arguments
OPTS=$(getopt -o s:m:c:r:a:vh --long sample-name:,mode:,vcf:,cache:,reference:,pop-af:,verbose,help -n "$0" -- "$@") || exit 1
eval set -- "$OPTS"

# setting variables based on inputted arguments
while true; do
  case "$1" in
    -s | --sample-name) SAMPLE_NAME="$2"; shift 2 ;;
    -m | --mode) MODE="$2"; shift 2 ;;
    -c | --cache) CACHE="$2"; shift 2 ;;
    --vcf) VCF="$2"; shift 2 ;;
    -r | --reference) REFERENCE="$2"; shift 2 ;;
    -a | --pop-af) POP_AF="$2"; shift 2 ;;
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

# --mode
if [[ "$MODE" != "tumor-only" && "$MODE" != "tumor-normal" && "$MODE" != "germline" ]]; then
    echo "[ERROR] must specify either --mode tumor-only, --mode tumor-normal or --mode germline"
    exit 1
fi

# --vcf 
if [[ -z "$VCF" ]]; then
    echo "[ERROR] --vcf is required"
    exit 1
fi

if [[ ! -f "$VCF" ]]; then
  echo "[ERROR] VCF file does not exist: $VCF"
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

# --cache
if [[ -z "$CACHE" ]]; then 
  echo "[ERROR] --cache is required"
  exit 1
fi

# --pop-af
if (( $(echo "$POP_AF > 1.0 || $POP_AF < 0.0" | bc) )); then 
    echo "[ERROR] --pop-af must be between 0.0 and 1.0 (got: $POP_AF)"
    exit 1
fi

# setting start 
SECONDS=0

# running annotation
annotation_workflow

# getting elapsed time
elapsed=$SECONDS
hours=$(( elapsed / 3600 ))
minutes=$(( (elapsed % 3600) / 60 ))
seconds=$(( elapsed % 60 ))

echo "Total runtime: ${hours}h ${minutes}m ${seconds}s"

}