#!/bin/bash
set -euo pipefail

# ------------------ user args ------------------
FASTQ_R1=""
FASTQ_R2=""
CONFIG=""
WHITELIST=""
SAMPLE="myGoT"
OUTDIR="./Output"
BCLEN=16
UMILEN=12
LEV_DIST=0.2
PCR_THRESH=0.75
DUPCUT=2
THREADS=8
TARGET_LINES=500000
KEEPOUTS=1
MAX_CHUNKS=0
QUEUE="standard"
PERL_BIN="/research/sharedresources/immunoinformatics/spoudel1/GoT_pipeline/pipeline/final_scripts_110225/IronThrone-GoT_only_preprocess"
PYTHON_BIN="python"
PYTHON_POST="$(dirname "$0")/got_post_process.py"
# ------------------------------------------------

usage() {
  echo "Usage: $0 --fastqR1 R1.fastq.gz --fastqR2 R2.fastq.gz --config cfg --whitelist wl.txt --sample S --outdir DIR [--python-bin python] [--python-post script.py]"
  exit 1
}

# ---- parse CLI ----
while [[ $# -gt 0 ]]; do
  case "$1" in
    --fastqR1) FASTQ_R1="$2"; shift 2;;
    --fastqR2) FASTQ_R2="$2"; shift 2;;
    --config) CONFIG="$2"; shift 2;;
    --whitelist) WHITELIST="$2"; shift 2;;
    --sample) SAMPLE="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --bclen) BCLEN="$2"; shift 2;;
    --umilen) UMILEN="$2"; shift 2;;
    --lev-dist|--ld) LEV_DIST="$2"; shift 2;;
    --pcr) PCR_THRESH="$2"; shift 2;;
    --dupcut) DUPCUT="$2"; shift 2;;
    --keepouts) KEEPOUTS="$2"; shift 2;;
    --threads|-t) THREADS="$2"; shift 2;;
    --target_lines) TARGET_LINES="$2"; shift 2;;
    --max_chunks) MAX_CHUNKS="$2"; shift 2;;
    --queue) QUEUE="$2"; shift 2;;
    --perl) PERL_BIN="$2"; shift 2;;
    --python-bin) PYTHON_BIN="$2"; shift 2;;
    --python-post) PYTHON_POST="$2"; shift 2;;
    *) echo "Unknown arg: $1"; usage;;
  esac
done

[[ -z "$FASTQ_R1" || -z "$FASTQ_R2" || -z "$CONFIG" || -z "$WHITELIST" ]] && usage

mkdir -p "$OUTDIR"
OUTDIR=$(readlink -f "$OUTDIR")

echo "[run] OUTDIR = $OUTDIR"

# 1) shuffle + split
 bash "$(dirname "$0")/01_shuffle_and_split.sh" \
   --r1 "$FASTQ_R1" \
   --r2 "$FASTQ_R2" \
   --target-lines "$TARGET_LINES" \
   --workdir "$OUTDIR" \
   --sample "$SAMPLE"

# 2) submit chunks (Perl pre + Python post in SAME job)
 bash "$(dirname "$0")/02_submit_chunks.sh" \
   --workdir "$OUTDIR" \
   --config "$CONFIG" \
   --whitelist "$WHITELIST" \
   --bclen "$BCLEN" \
   --umilen "$UMILEN" \
   --sample "$SAMPLE" \
   --queue "$QUEUE" \
   --perl "$PERL_BIN" \
   --python-bin "$PYTHON_BIN" \
   --python-post "$PYTHON_POST" \
   --max-chunks "$MAX_CHUNKS" \
   --dupcut "$DUPCUT"\
   --keepouts "$KEEPOUTS"

# # 3) wait for all chunks to finish
 bash "$(dirname "$0")/03_wait_for_chunks.sh" --workdir "$OUTDIR"


# 4) combine
$PYTHON_BIN "$(dirname "$0")/04_combine_chunks.py" \
  --outdir "$OUTDIR/Output" \
  --pcr-thresh "$PCR_THRESH" \
  --lev-dist "$LEV_DIST" \
  --dupcut "$DUPCUT" \
  --threads "$THREADS"




echo "[run] done."
