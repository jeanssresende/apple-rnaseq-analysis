#!/usr/bin/env bash

###############################################################################
#
# Apple RNA-seq Analysis Pipeline
#
# Module: 07 - Salmon Quantification
#
# Description:
# Quantifies transcript abundance from paired-end trimmed FASTQ files
# using Salmon and the Ensembl Plants Malus domestica transcriptome.
#
# Author: Jean Resende
#
###############################################################################

set -euo pipefail

########################################
# Input / Output
########################################

TRIM_DIR="../data/trimmed"

INDEX_DIR="../reference/Ensembl_Malus_domestica/salmon_index"

OUTPUT_DIR="../results/salmon_quant"

mkdir -p "$OUTPUT_DIR"

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
    THREADS=$((TOTAL_THREADS - 4))
fi

########################################
# Check index
########################################

if [ ! -f "${INDEX_DIR}/info.json" ]; then
    echo "ERROR: Salmon index not found:"
    echo "  ${INDEX_DIR}"
    exit 1
fi

########################################
# Detect samples
########################################

R1_FILES=( "${TRIM_DIR}"/*_R1.trimmed.fastq.gz )

if [ ! -e "${R1_FILES[0]}" ]; then
    echo "ERROR: No trimmed R1 FASTQ files found in ${TRIM_DIR}"
    exit 1
fi

N_SAMPLES=${#R1_FILES[@]}

########################################
# Header
########################################

echo "============================================================"
echo "           Apple RNA-seq Analysis Pipeline"
echo "             Module 07 - Salmon Quant"
echo "============================================================"
echo
echo "Input FASTQ  : ${TRIM_DIR}"
echo "Salmon index : ${INDEX_DIR}"
echo "Output dir   : ${OUTPUT_DIR}"
echo
echo "Samples found: ${N_SAMPLES}"
echo "Threads used : ${THREADS}"
echo
echo "Start time   : $(date)"
echo
echo "============================================================"

########################################
# Quantification loop
########################################

COUNT=0

for R1 in "${R1_FILES[@]}"
do

    SAMPLE=$(basename "$R1" _R1.trimmed.fastq.gz)

    R2="${TRIM_DIR}/${SAMPLE}_R2.trimmed.fastq.gz"

    if [ ! -f "$R2" ]; then
        echo "WARNING: Missing pair for sample ${SAMPLE}. Skipping."
        continue
    fi

    COUNT=$((COUNT + 1))

    SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE}"

    echo
    echo "------------------------------------------------------------"
    echo "[$COUNT/${N_SAMPLES}] Processing sample: ${SAMPLE}"
    echo "------------------------------------------------------------"
    echo

    salmon quant \
        -i "$INDEX_DIR" \
        -l A \
        -1 "$R1" \
        -2 "$R2" \
        -p "$THREADS" \
        --validateMappings \
        --gcBias \
        --seqBias \
        -o "$SAMPLE_OUT"

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
echo " Salmon quantification completed successfully!"
echo "============================================================"
echo
echo "Samples processed : ${COUNT}"
echo "Threads used      : ${THREADS}"
printf "Execution time    : %02d:%02d:%02d\n" \
       "$HOURS" "$MINUTES" "$SECONDS"
echo
echo "Results directory:"
echo "  ${OUTPUT_DIR}"
echo
echo "Finished at: $(date)"
echo "============================================================"