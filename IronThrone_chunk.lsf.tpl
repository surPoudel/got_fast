#BSUB -P GOT
#BSUB -q {{QUEUE}}
#BSUB -J got_{{SAMPLE}}_{{CHUNK_ID}}
#BSUB -n {{CORES}}
#BSUB -R "rusage[mem={{MEM_MB}}]"
#BSUB -W {{WALLTIME}}
#BSUB -oo {{WORKDIR}}/chunk_{{CHUNK_ID}}.out
#BSUB -eo {{WORKDIR}}/chunk_{{CHUNK_ID}}.err

set -euo pipefail

echo "[`date`] starting chunk {{CHUNK_ID}}"

# 1) run Perl pre_process (your IronThrone-GoT with post_process commented out)
#    IMPORTANT: only pass the options your Perl actually understands
{{PERL_BIN}} \
  --run {{RUNMODE}} \
  --fastqR1 {{R1}} \
  --fastqR2 {{R2}} \
  --config {{CONFIG}} \
  --whitelist {{WHITELIST}} \
  --sample {{SAMPLE}} \
  --outdir {{CHUNK_DIR}} \
  --bclen {{BC_LEN}} \
  --umilen {{UMI_LEN}} \
  --mmtch {{MMTCH}} \
  --postP {{POSTP}} \
  --dupcut {{DUPCUT}}

# At this point Perl should have written:
#   {{CHUNK_DIR}}/{{SAMPLE}}.looked

# 2) run Python fast post-process on that .looked
{{PY_BIN}} {{PY_POST}} \
  --run {{RUNMODE}} \
  --looked {{CHUNK_DIR}}/{{SAMPLE}}.looked \
  --config {{CONFIG}} \
  --whitelist {{WHITELIST}} \
  --sample {{SAMPLE}} \
  --outdir {{CHUNK_DIR}} \
  --mmtch {{MMTCH}} \
  --postP {{POSTP}} \
  --dupcut {{DUPCUT}} \
  --keepouts 0 \
  --verbose 1

# 3) mark done
touch {{CHUNK_DIR}}/DONE

echo "[`date`] finished chunk {{CHUNK_ID}}"








