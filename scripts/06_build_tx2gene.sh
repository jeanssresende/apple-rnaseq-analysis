#!/usr/bin/env bash

###############################################################################
#
# Apple RNA-seq Analysis Pipeline
#
# Module: 06 - Build tx2gene Table
#
# Description:
# Extracts transcript-to-gene relationships from the Ensembl Plants
# GFF3 annotation and generates a tx2gene.tsv file for tximport.
#
# Author: Jean Resende
#
###############################################################################

set -euo pipefail

########################################
# Input annotation
########################################

GFF3="../reference/Ensembl_Malus_domestica/annotation/Malus_domestica_golden.ASM211411v1.63.gff3.gz"

OUTPUT_DIR="../reference/Ensembl_Malus_domestica/tx2gene"

OUTPUT_FILE="${OUTPUT_DIR}/tx2gene.tsv"

########################################
# Create output directory
########################################

mkdir -p "$OUTPUT_DIR"

########################################
# Check annotation file
########################################

if [ ! -f "$GFF3" ]; then
    echo "ERROR: Annotation file not found:"
    echo "  $GFF3"
    exit 1
fi

########################################
# Header
########################################

echo "============================================================"
echo "           Apple RNA-seq Analysis Pipeline"
echo "             Module 06 - Build tx2gene"
echo "============================================================"
echo
echo "Input GFF3 : $GFF3"
echo "Output     : $OUTPUT_FILE"
echo
echo "============================================================"

########################################
# Extract transcript-gene mapping
########################################

echo "Extracting transcript-to-gene relationships..."
echo

zgrep -P "\tmRNA\t" "$GFF3" | \
awk -F '\t' '
{
    split($9, attrs, ";")

    tx=""
    gene=""

    for(i in attrs) {

        if(attrs[i] ~ /^transcript_id=/) {
            tx=attrs[i]
            sub(/^transcript_id=/, "", tx)
        }

        if(attrs[i] ~ /^Parent=gene:/) {
            gene=attrs[i]
            sub(/^Parent=gene:/, "", gene)
        }
    }

    if(tx != "" && gene != "") {
        print tx "\t" gene
    }
}' | sort -u > "$OUTPUT_FILE"

########################################
# Count entries
########################################

N_TX=$(wc -l < "$OUTPUT_FILE")

########################################
# Summary
########################################

echo
echo "============================================================"
echo " tx2gene table created successfully!"
echo "============================================================"
echo
echo "Transcript-gene pairs : $N_TX"
echo
echo "File:"
echo "  $OUTPUT_FILE"
echo
echo "Preview:"
head -5 "$OUTPUT_FILE"
echo
echo "============================================================"