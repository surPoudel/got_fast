#!/bin/bash
# ============================================================
# 00_config.sh — Central configuration for GoT hybrid pipeline
# ============================================================
# Defines all input, runtime, and cluster parameters
# Used by: run_got_pipeline.sh, 01_shuffle_and_split.sh,
#          02_submit_chunks.sh, 03_wait_for_chunks.sh, 04_combine_chunks.py
# ============================================================


# -------------------------
# INPUT FASTQs + CONFIG
# -------------------------
FASTQ_R1="/research/sharedresources/immunoinformatics/spoudel1/GoT_pipeline/pipeline/data/3364681/3364681_JCC365_RPS19_Ery_GoTseq_S2_L006_R1_001.fastq.gz"
FASTQ_R2="/research/sharedresources/immunoinformatics/spoudel1/GoT_pipeline/pipeline/data/3364681/3364681_JCC365_RPS19_Ery_GoTseq_S2_L006_R2_001.fastq.gz"
CONFIG_FILE="JCC365_RPS19-250825-1.config"
WHITELIST="3M-5pgex-jan-2023.txt"


# -------------------------
# IRONTHRONE (Perl PREPROCESS)
# -------------------------
# This Perl binary must have post_process() DISABLED/commented out
PERL_BIN="/research/sharedresources/immunoinformatics/spoudel1/GoT_pipeline/pipeline/final_scripts_110225/IronThrone-GoT_only_preprocess"


# -------------------------
# PYTHON POSTPROCESS (FAST)
# -------------------------
PY_BIN="python"
PY_POST="/research/sharedresources/immunoinformatics/spoudel1/GoT_pipeline/pipeline/final_scripts_110225/got_post_process.py"


# -------------------------
# OUTPUT ROOT DIRECTORY
# -------------------------
# You can auto-timestamp for clean re-runs
DATESTAMP=$(date +%Y%m%d_%H%M)
OUTROOT="/research/sharedresources/immunoinformatics/spoudel1/GoT_pipeline/pipeline/final_scripts_110225/run_S2_${DATESTAMP}"
mkdir -p "$OUTROOT"


# -------------------------
# SPLIT / READ PARAMETERS
# -------------------------
READS_PER_CHUNK=500000      # 500k lines = 125k reads/chunk
BC_LEN=16
UMI_LEN=12
MISMATCH=0.2
POSTP=0.99
DUPCUT=1                    # 1 = strict; 2 = lenient
RUNMODE="linear"
SAMPLE_NAME="3364681_S2"


# -------------------------
# LSF (CLUSTER) SETTINGS
# -------------------------
# !! Confirm your queue allows this walltime !!
# For 4h runtime, use "standard" if "short" ≤1h
LSF_QUEUE="standard"        # or "short"
LSF_WALL="04:00"            # hours:minutes
LSF_CORES=4                 # requested CPU cores
LSF_MEM=64000               # MB per job (64 GB)


# -------------------------
# BATCHING CONTROL
# -------------------------
# How many chunks per LSF job (1 = per-chunk jobs)
CHUNKS_PER_JOB=10


# -------------------------
# COMBINE (PYTHON MERGE)
# -------------------------
PCR_THRESH=0.75
LEV_DIST=0.2
COMBINE_THREADS=8


# -------------------------
# JOB CAP (FOR TESTING)
# -------------------------
MAX_CHUNKS=10               # 0 = all chunks; >0 = limit for test runs


# -------------------------
# EXPORT ALL VARIABLES
# -------------------------
export FASTQ_R1 FASTQ_R2 CONFIG_FILE WHITELIST \
       PERL_BIN PY_BIN PY_POST OUTROOT DATESTAMP \
       READS_PER_CHUNK BC_LEN UMI_LEN MISMATCH POSTP DUPCUT \
       RUNMODE SAMPLE_NAME LSF_QUEUE LSF_WALL LSF_CORES LSF_MEM \
       CHUNKS_PER_JOB PCR_THRESH LEV_DIST COMBINE_THREADS MAX_CHUNKS

# ============================================================
# END OF CONFIG
# ============================================================

