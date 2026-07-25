#!/usr/bin/env bash

###############################################################################
#
# Apple RNA-seq Analysis Pipeline
#
# Module: 03 - Quality Control After Trimming
#
# Description:
# Performs quality assessment of trimmed paired-end RNA-Seq reads
# using FastQC and MultiQC.
#
# Author: Jean Resende
#
###############################################################################

set -euo pipefail

########################################
# Directories
########################################

INPUT_DIR="../data/trimmed"

FASTQC_DIR="../results/quality_control_trimmed/fastqc"

MULTIQC_DIR="../results/quality_control_trimmed/multiqc"

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

mkdir -p "$FASTQC_DIR"
mkdir -p "$MULTIQC_DIR"

########################################
# Locate FASTQ files
########################################

FASTQ_FILES=("$INPUT_DIR"/*.trimmed.fastq.gz)

if [ ! -e "${FASTQ_FILES[0]}" ]; then
    echo "ERROR: No trimmed FASTQ files found."
    exit 1
fi

N_FILES=${#FASTQ_FILES[@]}

########################################
# Header
########################################

echo "============================================================"
echo "        Apple RNA-seq Analysis Pipeline"
echo " Module 03 - Quality Control After Trimming"
echo "============================================================"

echo
echo "Input directory  : ${INPUT_DIR}"
echo "FastQC output    : ${FASTQC_DIR}"
echo "MultiQC output   : ${MULTIQC_DIR}"

echo
echo "FASTQ files      : ${N_FILES}"
echo "Threads          : ${THREADS}"

echo
echo "Start time       : $(date)"
echo

########################################
# FastQC
########################################

echo "Running FastQC..."

fastqc \
    --threads "$THREADS" \
    --outdir "$FASTQC_DIR" \
    "${FASTQ_FILES[@]}"

########################################
# MultiQC
########################################

echo
echo "Running MultiQC..."

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
echo " Quality control completed successfully!"
echo "============================================================"

echo
echo "FASTQ files processed : ${N_FILES}"
echo "Threads used          : ${THREADS}"

printf "Execution time        : %02d:%02d:%02d\n" \
       "$HOURS" "$MINUTES" "$SECONDS"

echo
echo "Results:"
echo "  FastQC : ${FASTQC_DIR}"
echo "  MultiQC: ${MULTIQC_DIR}"

echo
echo "Finished at: $(date)"
echo "============================================================"