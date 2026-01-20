#!/bin/bash
set -euo pipefail

WORKDIR="."
SLEEP_SEC=20

while [[ $# -gt 0 ]]; do
  case "$1" in
    --workdir) WORKDIR="$2"; shift 2;;
    --sleep)   SLEEP_SEC="$2"; shift 2;;
    *) echo "unknown arg: $1"; exit 1;;
  esac
done

WORKDIR="$(readlink -f "$WORKDIR")"
OUT_DIR="${WORKDIR}/Output"
BATCH_LIST_FILE="${OUT_DIR}/expected_batches.txt"

normalize_batch() {
  local b="$1"
  b="${b//$'\r'/}"      # strip CR if file has Windows line endings
  b="${b##*/}"          # strip any path
  b="${b#batch_}"       # strip prefix if present
  b="${b%.DONE}"        # strip suffix if present
  b="${b%.done}"
  b="${b%.lsf}"
  echo "$b"
}

has_done_batch() {
  local raw="$1"
  local b
  b="$(normalize_batch "$raw")"

  if [[ "$b" =~ ^[0-9]+$ ]]; then
    local b4
    b4="$(printf "%04d" "$b")"
    [[ -f "${OUT_DIR}/batch_${b}.DONE" || -f "${OUT_DIR}/batch_${b4}.DONE" ]]
  else
    [[ -f "${OUT_DIR}/batch_${b}.DONE" ]]
  fi
}

[[ -s "$BATCH_LIST_FILE" ]] || { echo "[03] Error: missing/empty $BATCH_LIST_FILE"; exit 1; }

mapfile -t BATCHES < "$BATCH_LIST_FILE"
N=${#BATCHES[@]}
echo "[03] will wait for ${N} batches to finish"

while true; do
  done_cnt=0
  for b in "${BATCHES[@]}"; do
    if has_done_batch "$b"; then
      ((++done_cnt))   # SAFE with set -e
    fi
  done

  echo "[03] ${done_cnt}/${N} batches done"

  if (( done_cnt == N )); then
    echo "[03] all batches done"
    break
  fi

  sleep "$SLEEP_SEC"
done

