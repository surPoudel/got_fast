#!/bin/bash
set -euo pipefail

# ---- 1. Initialize Defaults (Matching Original Script Variable Names) ----
fastqR1=""
fastqR2=""
target_lines=500000
low_mem=0
WORKDIR="."

# ---- 2. Parse User Arguments ----
while [[ $# -gt 0 ]]; do
  case "$1" in
    --r1) fastqR1="$2"; shift 2;;
    --r2) fastqR2="$2"; shift 2;;
    --target-lines) target_lines="$2"; shift 2;;
    --workdir) WORKDIR="$2"; shift 2;;
    --sample) sample="$2"; shift 2;;
    --low-mem) low_mem=1; shift 1;;
    *) echo "Unknown arg $1"; exit 1;;
  esac
done

# ---- 3. Setup Environment ----
[[ -z "$fastqR1" || -z "$fastqR2" ]] && { echo "ERROR: R1/R2 FASTQs not found."; exit 1; }

mkdir -p "$WORKDIR"
fastqR1=$(readlink -f "$fastqR1")
fastqR2=$(readlink -f "$fastqR2")

cd "$WORKDIR"

# Convert desired line number of split fastq into number of reads
target_reads=$((target_lines / 4))

# ==================== SHUFFLE & SPLIT (EXACT COPY) ====================
# The code below is pasted verbatim from IronThroneParLSF.sh

# join gz or plain
if [[ "$fastqR1" == *.gz ]]; then
  paste <(zcat "$fastqR1") <(zcat "$fastqR2") > combined.fastq
else
  paste "$fastqR1" "$fastqR2" > combined.fastq
fi
echo "fastq files joined"

if (( low_mem == 1 )); then
  if (($(grep ";" combined.fastq | wc -l) == 0)); then
    awk '{printf("%s%s",$0,(NR%4==0)?"\n":";")}' combined.fastq | sort -R | tr ";" "\n" > combined_shuffled.fastq
  elif (($(grep "|" combined.fastq | wc -l) == 0)); then
    awk '{printf("%s%s",$0,(NR%4==0)?"\n":"|")}' combined.fastq | sort -R | tr "|" "\n" > combined_shuffled.fastq
  else
    echo "New awk-line character needed"; exit 1
  fi
else
  if (($(grep ";" combined.fastq | wc -l) == 0)); then
    awk '{printf("%s%s",$0,(NR%4==0)?"\n":";")}' combined.fastq | shuf | tr ";" "\n" > combined_shuffled.fastq
  elif (($(grep "|" combined.fastq | wc -l) == 0)); then
    awk '{printf("%s%s",$0,(NR%4==0)?"\n":"|")}' combined.fastq | shuf | tr "|" "\n" > combined_shuffled.fastq
  else
    echo "New awk-line character needed"; exit 1
  fi
fi
echo "fastq files shuffled"

cut -f1 -d$'\t' combined_shuffled.fastq > shuffled.R1.fastq
cut -f2 -d$'\t' combined_shuffled.fastq > shuffled.R2.fastq
echo "shuffled fastq files cut back into R1 and R2"

total_lines=$(wc -l < combined_shuffled.fastq)
total_reads=$((total_lines / 4))
reads_mod=$((total_reads % target_reads))
if (( reads_mod*10 < target_reads*9 && reads_mod != 0 )); then
  remainder_reads=$((total_reads % target_reads))
  total_files=$((total_reads / target_reads))
  to_add=$(((remainder_reads / total_files) + 1))
  target_reads=$((target_reads + to_add))
  target_lines=$((target_reads * 4))
fi

split -d -a 4 -l $target_lines shuffled.R1.fastq shuffled.R1
split -d -a 4 -l $target_lines shuffled.R2.fastq shuffled.R2

for file in $(ls | grep '.*R[0-9][0-9][0-9][0-9][0-9]'); do mv "$file" "$file.fastq"; done
mkdir -p shuffled_split preprocessing_fastqs
for file in $(ls | grep '.*R[0-9][0-9][0-9][0-9][0-9]'); do mv "$file" "./shuffled_split/"; done
mv combined.fastq combined_shuffled.fastq shuffled.R1.fastq shuffled.R2.fastq preprocessing_fastqs/

echo "[DONE] Shuffle and split complete in $WORKDIR"
