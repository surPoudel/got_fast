#BSUB -P GOT
#BSUB -q {{QUEUE}}
#BSUB -J got_{{SAMPLE}}_batch{{BATCH_INDEX}}
#BSUB -n {{CORES}}
#BSUB -R "rusage[mem={{MEM_MB}}]"
#BSUB -W {{WALLTIME}}
#BSUB -oo {{WORKDIR}}/batch_{{BATCH_INDEX}}.out
#BSUB -eo {{WORKDIR}}/batch_{{BATCH_INDEX}}.err

set -euo pipefail

echo "[`date`] starting batch {{BATCH_INDEX}}"

SPLIT_DIR="{{SPLIT_DIR}}"
WORKDIR="{{WORKDIR}}"
OUTROOT="${WORKDIR}/Output"
PERL_BIN="{{PERL_BIN}}"
PY_BIN="{{PY_BIN}}"
PY_POST="{{PY_POST}}"

CONFIG="{{CONFIG}}"
WHITELIST="{{WHITELIST}}"
SAMPLE="{{SAMPLE}}"

RUNMODE="{{RUNMODE}}"
BC_LEN="{{BC_LEN}}"
UMI_LEN="{{UMI_LEN}}"
MMTCH="{{MMTCH}}"
POSTP="{{POSTP}}"
DUPCUT="{{DUPCUT}}"
KEEPOUTS="{{KEEPOUTS}}"
CHUNK_LIST="{{CHUNK_LIST}}"

mkdir -p "${OUTROOT}"

# Iterate ids in this batch
for id in ${CHUNK_LIST}; do
  R1="${SPLIT_DIR}/shuffled.R1${id}.fastq"
  R2="${SPLIT_DIR}/shuffled.R2${id}.fastq"
  [[ -f "${R1}.gz" ]] && R1="${R1}.gz"
  [[ -f "${R2}.gz" ]] && R2="${R2}.gz"

  if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "[batch {{BATCH_INDEX}}] missing R1/R2 for id=${id}; skipping"
    continue
  fi

  OUTDIR="${OUTROOT}/${id}"
  mkdir -p "${OUTDIR}"

  echo "[`date`] batch {{BATCH_INDEX}} id=${id} → Perl pre-process"
  "${PERL_BIN}" \
    --run "${RUNMODE}" \
    --fastqR1 "${R1}" \
    --fastqR2 "${R2}" \
    --config "${CONFIG}" \
    --whitelist "${WHITELIST}" \
    --sample "${SAMPLE}" \
    --outdir "${OUTDIR}" \
    --bclen "${BC_LEN}" \
    --umilen "${UMI_LEN}" \
    --mmtch "${MMTCH}" \
    --postP "${POSTP}" \
    --dupcut "${DUPCUT}"

  LOOKED="${OUTDIR}/${SAMPLE}.looked"
  if [[ ! -s "${LOOKED}" ]]; then
    echo "[batch {{BATCH_INDEX}}] id=${id} → missing ${LOOKED}; skipping post-process"
    continue
  fi

  echo "[`date`] batch {{BATCH_INDEX}} id=${id} → Python post-process"
  "${PY_BIN}" "${PY_POST}" \
    --run "${RUNMODE}" \
    --looked "${LOOKED}" \
    --config "${CONFIG}" \
    --whitelist "${WHITELIST}" \
    --sample "${SAMPLE}" \
    --outdir "${OUTDIR}" \
    --mmtch "${MMTCH}" \
    --postP "${POSTP}" \
    --dupcut "${DUPCUT}" \
    --keepouts "${KEEPOUTS}" \
    --verbose 1

  touch "${OUTDIR}/DONE"
done

# Batch-level sentinel
touch "${OUTROOT}/batch_{{BATCH_INDEX}}.DONE"
echo "[`date`] finished batch {{BATCH_INDEX}}"

