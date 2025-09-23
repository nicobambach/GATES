#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# setting log file 
LOG_FILE="./$(date +%Y-%m-%d_%H-%M-%S)_gates_${GATES_COMMAND:-unknown}.log"

# run_cmd: silences tool outputs unless VERBOSE=1
run_cmd() {
echo "[$(date '+%F %T')] Running: $*" >> "$LOG_FILE"

  if [[ ${VERBOSE:-0} -eq 1 ]]; then
    "$@" 2>&1 | tee -a "$LOG_FILE"
  else
    "$@" >> "$LOG_FILE" 2>&1
  fi
}

# log: always prints timestamped messages
log() {
  echo "[$(date '+%F %T')] $*" >&2
  echo "[$(date '+%F %T')] $*" >> "$LOG_FILE"
}