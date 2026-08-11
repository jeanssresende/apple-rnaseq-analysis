#!/usr/bin/env bash

###############################################################################
#
# Apple RNA-seq Analysis Pipeline
#
# Module: 04 - Download Reference Files
#
# Description:
# Downloads the Ensembl Plants transcriptome FASTA and GFF3 annotation
# for Malus domestica (ASM211411v1, release 63).
#
# Author: Jean Resende
#
###############################################################################

set -euo pipefail

########################################
# Directories
########################################

BASE_DIR="../reference/Ensembl_Malus_domestica"
TRANSCRIPTOME_DIR="${BASE_DIR}/transcriptome"
ANNOTATION_DIR="${BASE_DIR}/annotation"

mkdir -p "$TRANSCRIPTOME_DIR"
mkdir -p "$ANNOTATION_DIR"

########################################
# URLs
########################################

TRANSCRIPTOME_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/fasta/malus_domestica_golden/cdna/Malus_domestica_golden.ASM211411v1.cdna.all.fa.gz"

ANNOTATION_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-63/gff3/malus_domestica_golden/Malus_domestica_golden.ASM211411v1.63.gff3.gz"

########################################
# Output files
########################################

TRANSCRIPTOME_FILE="${TRANSCRIPTOME_DIR}/Malus_domestica_golden.ASM211411v1.cdna.all.fa.gz"

ANNOTATION_FILE="${ANNOTATION_DIR}/Malus_domestica_golden.ASM211411v1.63.gff3.gz"

########################################
# Download transcriptome
########################################

echo "Downloading transcriptome..."

curl -L --fail --retry 3 --continue-at - \
     -o "$TRANSCRIPTOME_FILE" \
     "$TRANSCRIPTOME_URL"

########################################
# Download annotation
########################################

echo
echo "Downloading annotation..."

curl -L --fail --retry 3 --continue-at - \
     -o "$ANNOTATION_FILE" \
     "$ANNOTATION_URL"

########################################
# Summary
########################################

echo
echo "============================================================"
echo " Reference download completed successfully!"
echo "============================================================"
echo
echo "Transcriptome:"
echo "  $TRANSCRIPTOME_FILE"
echo
echo "Annotation:"
echo "  $ANNOTATION_FILE"
echo
echo "============================================================"