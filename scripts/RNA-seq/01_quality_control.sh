#!/bin/bash
# ============================================================================
# Script: 01_quality_control.sh
# Purpose: Run FastQC on raw FASTQ files
# Input:   All .gz files in current directory
# Output:  FastQC HTML reports and ZIP archives
# ============================================================================

# Activate conda environment (adjust environment name as needed)
source activate rnaseq

# Run FastQC on all compressed FASTQ files in current directory
ls ./*.gz | xargs fastqc -o ./

# Deactivate environment
conda deactivate

echo "FastQC analysis completed. Check the HTML reports in the current directory."
