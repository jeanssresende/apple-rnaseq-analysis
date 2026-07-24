#!/usr/bin/env bash

###############################################################################
#
# Apple RNA-seq Analysis Pipeline
#
# Module: 01 - Quality Control
#
# Description:
# Performs quality assessment of raw FASTQ files using FastQC and
# summarizes the results with MultiQC.
#
# Author: Jean Resende
#
###############################################################################

set -euo pipefail

########################################
# Directories
########################################

RAW_DIR="../../data/raw"
FASTQC_DIR="../results/qc/fastqc"
MULTIQC_DIR="../results/qc/multiqc"

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
# Create output directories
########################################

mkdir -p "$FASTQC_DIR"
mkdir -p "$MULTIQC_DIR"

########################################
# Locate FASTQ files
########################################

FASTQ_FILES=("$RAW_DIR"/*.fastq.gz)

if [ ! -e "${FASTQ_FILES[0]}" ]; then
    echo "ERROR: No FASTQ files found in ${RAW_DIR}"
    exit 1
fi

N_SAMPLES=${#FASTQ_FILES[@]}

########################################
# Header
########################################

echo "============================================================"
echo "           Apple RNA-seq Analysis Pipeline"
echo "             Module 01 - Quality Control"
echo "============================================================"
echo
echo "Input directory  : ${RAW_DIR}"
echo "FastQC output    : ${FASTQC_DIR}"
echo "MultiQC output   : ${MULTIQC_DIR}"
echo
echo "Samples detected : ${N_SAMPLES}"
echo "CPUs detected    : ${TOTAL_THREADS}"
echo "Threads in use   : ${THREADS}"
echo
echo "Start time       : $(date)"
echo
echo "============================================================"

########################################
# Run FastQC
########################################

echo
echo "[1/2] Running FastQC..."
echo

fastqc \
    --threads "$THREADS" \
    --outdir "$FASTQC_DIR" \
    "${FASTQ_FILES[@]}"

########################################
# Run MultiQC
########################################

echo
echo "[2/2] Running MultiQC..."
echo

multiqc \
    "$FASTQC_DIR" \
    --outdir "$MULTIQC_DIR" \
    --force

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
echo " Quality Control completed successfully!"
echo "============================================================"
echo
echo "Samples processed : ${N_SAMPLES}"
echo "Threads used      : ${THREADS}"
printf "Execution time    : %02d:%02d:%02d\n" \
       "$HOURS" "$MINUTES" "$SECONDS"
echo
echo "Results:"
echo "  FastQC  -> ${FASTQC_DIR}"
echo "  MultiQC -> ${MULTIQC_DIR}"
echo
echo "Finished at: $(date)"
echo "============================================================"