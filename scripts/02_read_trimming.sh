#!/usr/bin/env bash

###############################################################################
#
# Apple RNA-seq Analysis Pipeline
#
# Module: 02 - Read Trimming
#
# Description:
# Performs adapter removal and quality trimming of paired-end RNA-Seq reads
# using fastp.
#
# Author: Jean Resende
#
###############################################################################

set -euo pipefail

########################################
# Directories
########################################

RAW_DIR="../../data/raw"

TRIM_DIR="../data/trimmed"

REPORT_DIR="../results/trimming"

########################################
# Start timer
########################################

START_TIME=$(date +%s)

########################################
# Detect CPUs
########################################

TOTAL_THREADS=$(nproc)

if [ "$TOTAL_THREADS" -le 4 ]; then
    THREADS=$TOTAL_THREADS
elif [ "$TOTAL_THREADS" -le 16 ]; then
    THREADS=$((TOTAL_THREADS - 1))
else
    THREADS=$((TOTAL_THREADS - 2))
fi

########################################
# Create directories
########################################

mkdir -p "$TRIM_DIR"
mkdir -p "$REPORT_DIR"

########################################
# Locate R1 files
########################################

R1_FILES=("$RAW_DIR"/*_R1_*.fastq.gz)

if [ ! -e "${R1_FILES[0]}" ]; then
    echo "ERROR: No paired-end FASTQ files found."
    exit 1
fi

N_SAMPLES=${#R1_FILES[@]}

########################################
# Header
########################################

echo "============================================================"
echo "           Apple RNA-seq Analysis Pipeline"
echo "              Module 02 - Read Trimming"
echo "============================================================"
echo
echo "Input directory  : ${RAW_DIR}"
echo "Output directory : ${TRIM_DIR}"
echo "Reports          : ${REPORT_DIR}"
echo
echo "Samples detected : ${N_SAMPLES}"
echo "CPUs detected    : ${TOTAL_THREADS}"
echo "Threads in use   : ${THREADS}"
echo
echo "Start time       : $(date)"
echo
echo "============================================================"

########################################
# Trimming
########################################

COUNT=1

for R1 in "${R1_FILES[@]}"
do

    R2="${R1/_R1_/_R2_}"

    SAMPLE=$(basename "$R1")
    SAMPLE=${SAMPLE%%_R1_*}

    echo
    echo "------------------------------------------------------------"
    echo "[$COUNT/$N_SAMPLES] Processing sample: ${SAMPLE}"
    echo "------------------------------------------------------------"

    fastp \
        --in1 "$R1" \
        --in2 "$R2" \
        --out1 "${TRIM_DIR}/${SAMPLE}_R1.trimmed.fastq.gz" \
        --out2 "${TRIM_DIR}/${SAMPLE}_R2.trimmed.fastq.gz" \
        --detect_adapter_for_pe \
        --cut_front \
        --cut_right \
        --cut_window_size 4 \
        --cut_mean_quality 20 \
        --length_required 50 \
        --thread "$THREADS" \
        --html "${REPORT_DIR}/${SAMPLE}.fastp.html" \
        --json "${REPORT_DIR}/${SAMPLE}.fastp.json"

    ((COUNT++))

done

########################################
# Finish timer
########################################

END_TIME=$(date +%s)

ELAPSED=$((END_TIME - START_TIME))

HOURS=$((ELAPSED / 3600))
MINUTES=$(((ELAPSED % 3600) / 60))
SECONDS=$((ELAPSED % 60))

########################################
# Summary
########################################

echo
echo "============================================================"
echo " Read trimming completed successfully!"
echo "============================================================"
echo
echo "Samples processed : ${N_SAMPLES}"
echo "Threads used      : ${THREADS}"

printf "Execution time    : %02d:%02d:%02d\n" \
       "$HOURS" "$MINUTES" "$SECONDS"

echo
echo "Trimmed reads:"
echo "  ${TRIM_DIR}"

echo
echo "Reports:"
echo "  ${REPORT_DIR}"

echo
echo "Finished at: $(date)"
echo "============================================================"