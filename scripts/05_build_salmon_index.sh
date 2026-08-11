#!/usr/bin/env bash

###############################################################################
#
# Apple RNA-seq Analysis Pipeline
#
# Module: 05 - Build Salmon Index
#
# Description:
# Builds a Salmon transcriptome index using the Ensembl Plants
# Malus domestica transcriptome reference.
#
# Author: Jean Resende
#
###############################################################################

set -euo pipefail
set -x

########################################
# Reference files
########################################

TRANSCRIPTOME="../reference/Ensembl_Malus_domestica/transcriptome/Malus_domestica_golden.ASM211411v1.cdna.all.fa"

INDEX_DIR="../reference/Ensembl_Malus_domestica/salmon_index"

########################################
# Start timer
########################################

START_TIME=$(date +%s)

########################################
# Detect CPUs
########################################

TOTAL_THREADS=$(nproc)

# Use no máximo 8 threads para indexação
if [ "$TOTAL_THREADS" -ge 8 ]; then
    THREADS=8
else
    THREADS=$TOTAL_THREADS
fi

########################################
# Check transcriptome
########################################

if [ ! -f "$TRANSCRIPTOME" ]; then
    echo "ERROR: Transcriptome file not found:"
    echo "  $TRANSCRIPTOME"
    exit 1
fi

########################################
# Create output directory
########################################

mkdir -p "$INDEX_DIR"

########################################
# Header
########################################

echo "============================================================"
echo "           Apple RNA-seq Analysis Pipeline"
echo "             Module 05 - Salmon Index"
echo "============================================================"
echo
echo "Transcriptome : $TRANSCRIPTOME"
echo "Index output  : $INDEX_DIR"
echo
echo "CPUs detected : $TOTAL_THREADS"
echo "Threads in use: $THREADS"
echo
echo "Start time    : $(date)"
echo
echo "============================================================"

########################################
# Build Salmon index
########################################

echo
echo "Building Salmon index..."
echo

salmon index \
    -t "$TRANSCRIPTOME" \
    -i "$INDEX_DIR" \
    -p "$THREADS"

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
echo " Salmon index completed successfully!"
echo "============================================================"
echo
printf "Execution time : %02d:%02d:%02d\n" \
       "$HOURS" "$MINUTES" "$SECONDS"
echo
echo "Index location:"
echo "  $INDEX_DIR"
echo
echo "Finished at: $(date)"
echo "============================================================"