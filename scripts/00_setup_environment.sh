#!/usr/bin/env bash

###############################################################################
# Apple RNA-seq Analysis Pipeline
#
# Script: 00_setup_environment.sh
#
# Description:
# Creates the Conda environment and installs all required software.
#
# Author: Jean Resende
###############################################################################

set -e
set -o pipefail

ENV_NAME="apple-rnaseq"

echo "Creating Conda environment..."

conda create \
    -n ${ENV_NAME} \
    -c conda-forge \
    -c bioconda \
    fastqc \
    multiqc \
    fastp \
    salmon \
    python=3.12 \
    -y

echo
echo "Environment created successfully!"
echo
echo "Activate using:"
echo
echo "conda activate ${ENV_NAME}"