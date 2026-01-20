#!/bin/bash
set -euo pipefail

WORKDIR="."
CONFIG_FILE=""
WHITELIST=""
BCLEN=16
UMI_LEN=12
SAMPLE_NAME="myGoT"
LSF_QUEUE="standard"
PERL_BIN="/path/to/IronThrone-GoT_only_preprocess"
PY_BIN="python"
PY_POST="got_post_process.py"
MAX_CHUNKS=0
RUNMODE="linear"
MISMATCH=0.2
POSTP=0.99
DUPCUT=1
KEEPOUTS=1
LSF_CORES=2
LSF_MEM=64000
LSF_WALL="04:00"
CHUNKS_PER_JOB=10

usage() {
  echo "Usage: $0 --workdir DIR --config cfg --whitelist wl.txt --sample NAME [--max-chunks N] [--chunks-per-job N]"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --workdir) WORKDIR="$2"; shift 2;;
    --config) CONFIG_FILE="$2"; shift 2;;
    --whitelist) WHITELIST="$2"; shift 2;;
    --bclen) BCLEN="$2"; shift 2;;
    --umilen) UMI_LEN="$2"; shift 2;;
    --sample) SAMPLE_NAME="$2"; shift 2;;
    --queue) LSF_QUEUE="$2"; shift 2;;
    --perl) PERL_BIN="$2"; shift 2;;
    --python-bin) PY_BIN="$2"; shift 2;;
    --python-post) PY_POST="$2"; shift 2;;
    --max-chunks) MAX_CHUNKS="$2"; shift 2;;
    --runmode) RUNMODE="$2"; shift 2;;
    --mismatch) MISMATCH="$2"; shift 2;;
    --postP) POSTP="$2"; shift 2;;
    --dupcut) DUPCUT="$2"; shift 2;;
    --keepouts) KEEPOUTS="$2"; shift 2;;
    --lsf-cores) LSF_CORES="$2"; shift 2;;
    --lsf-mem) LSF_MEM="$2"; shift 2;;
    --lsf-wall) LSF_WALL="$2"; shift 2;;
    --chunks-per-job) CHUNKS_PER_JOB="$2"; shift 2;;
    *) echo "Unknown arg: $1"; usage;;
  esac
done

WORKDIR="$(readlink -f "$WORKDIR")"
CONFIG_FILE="$(readlink -f "$CONFIG_FILE")"
WHITELIST="$(readlink -f "$WHITELIST")"

SPLIT_DIR="${WORKDIR}/shuffled_split"
BATCH_TEMPLATE="$(dirname "$0")/IronThrone_batch.lsf.tpl"

# discover chunk IDs from shuffled.R1<digits>.fastq[.gz]
mapfile -t CHUNK_IDS < <(
  ls -1 "${SPLIT_DIR}"/shuffled.R1[0-9]*.fastq* 2>/dev/null \
    | sed -E 's|.*/shuffled\.R1([0-9]+)\.fastq(\.gz)?$|\1|' \
    | sort -n
)
[[ ${#CHUNK_IDS[@]} -eq 0 ]] && { echo "[submit] no chunks found in ${SPLIT_DIR}"; exit 1; }

# ensure R2 exists
filtered=()
for id in "${CHUNK_IDS[@]}"; do
  if [[ -s "${SPLIT_DIR}/shuffled.R2${id}.fastq" || -s "${SPLIT_DIR}/shuffled.R2${id}.fastq.gz" ]]; then
    filtered+=( "$id" )
  fi
done
CHUNK_IDS=( "${filtered[@]}" )
[[ ${#CHUNK_IDS[@]} -eq 0 ]] && { echo "[submit] no paired chunks in ${SPLIT_DIR}"; exit 1; }

# cap
if [[ "${MAX_CHUNKS}" -gt 0 && "${#CHUNK_IDS[@]}" -gt "${MAX_CHUNKS}" ]]; then
  CHUNK_IDS=("${CHUNK_IDS[@]:0:${MAX_CHUNKS}}")
fi

echo "[submit] submitting ${#CHUNK_IDS[@]} chunks from ${SPLIT_DIR}"

# canonicalize tools
PERL_BIN="$(readlink -f "$PERL_BIN")"
PY_BIN="$(command -v "$PY_BIN")"
PY_POST="$(readlink -f "$PY_POST")"
[[ -x "$PERL_BIN" ]] || { echo "PERL_BIN not executable: $PERL_BIN"; exit 2; }
[[ -n "$PY_BIN" ]] || { echo "PY_BIN not found in PATH"; exit 2; }
[[ -f "$PY_POST" ]] || { echo "PY_POST not found: $PY_POST"; exit 2; }
[[ -f "$BATCH_TEMPLATE" ]] || { echo "Batch template missing: $BATCH_TEMPLATE"; exit 2; }

mkdir -p "${WORKDIR}/Output"
BATCH_LIST_FILE="${WORKDIR}/Output/expected_batches.txt"
: > "$BATCH_LIST_FILE"

# batching
total="${#CHUNK_IDS[@]}"
batches=$(( (total + CHUNKS_PER_JOB - 1) / CHUNKS_PER_JOB ))
echo "[submit] scheduling ${batches} jobs (${CHUNKS_PER_JOB} chunks/job), queue=${LSF_QUEUE}, wall=${LSF_WALL}"

for ((b=0; b<batches; b++)); do
  start=$(( b * CHUNKS_PER_JOB ))
  end=$(( start + CHUNKS_PER_JOB ))
  (( end > total )) && end=$total
  ids=( "${CHUNK_IDS[@]:start:end-start}" )

  CHUNK_LIST="${ids[*]}"
  BATCH_INDEX=$(printf "%04d" "$b")
  LSF_FILE="${WORKDIR}/Output/batch_${BATCH_INDEX}.lsf"

  sed \
    -e "s|{{QUEUE}}|${LSF_QUEUE}|g" \
    -e "s|{{SAMPLE}}|${SAMPLE_NAME}|g" \
    -e "s|{{BATCH_INDEX}}|${BATCH_INDEX}|g" \
    -e "s|{{CORES}}|${LSF_CORES}|g" \
    -e "s|{{MEM_MB}}|${LSF_MEM}|g" \
    -e "s|{{WALLTIME}}|${LSF_WALL}|g" \
    -e "s|{{WORKDIR}}|${WORKDIR}|g" \
    -e "s|{{SPLIT_DIR}}|${SPLIT_DIR}|g" \
    -e "s|{{PERL_BIN}}|${PERL_BIN}|g" \
    -e "s|{{PY_BIN}}|${PY_BIN}|g" \
    -e "s|{{PY_POST}}|${PY_POST}|g" \
    -e "s|{{CONFIG}}|${CONFIG_FILE}|g" \
    -e "s|{{WHITELIST}}|${WHITELIST}|g" \
    -e "s|{{BC_LEN}}|${BCLEN}|g" \
    -e "s|{{UMI_LEN}}|${UMI_LEN}|g" \
    -e "s|{{MMTCH}}|${MISMATCH}|g" \
    -e "s|{{POSTP}}|${POSTP}|g" \
    -e "s|{{DUPCUT}}|${DUPCUT}|g" \
    -e "s|{{KEEPOUTS}}|${KEEPOUTS}|g" \
    -e "s|{{RUNMODE}}|${RUNMODE}|g" \
    -e "s|{{CHUNK_LIST}}|${CHUNK_LIST}|g" \
    "$BATCH_TEMPLATE" > "$LSF_FILE"

  bsub < "$LSF_FILE"
  echo "$BATCH_INDEX" >> "$BATCH_LIST_FILE"
done

