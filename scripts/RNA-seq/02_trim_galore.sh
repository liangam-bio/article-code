#!/bin/bash
# ============================================================================
# Script: 02_trim_galore.sh
# Purpose: Trim adapters and low-quality bases using Trim Galore
# Input:   config file with paired-end FASTQ paths
# Output:  Trimmed FASTQ files (.clean.fq.gz)
# ============================================================================

# Activate conda environment
source activate rnaseq

# Output directory
dir='./'

# ============================================================================
# STEP 1: Create config file (uncomment and run once, or prepare manually)
# ============================================================================
# ls ./*1.fq.gz > 1
# ls ./*2.fq.gz > 2
# paste 1 2 > config

# ============================================================================
# STEP 2: Run Trim Galore on all sample pairs
# ============================================================================
cat config | while read id
do
    # Split the line into forward and reverse read files
    arr=(${id})
    fq1=${arr[0]}
    fq2=${arr[1]} 
    
    # Run trimming
    # -q 30: quality threshold
    # --phred33: PHRED quality encoding
    # --length 36: minimum read length after trimming
    # -e 0.1: maximum error rate
    # --stringency 3: adapter stringency
    # --paired: paired-end mode
    nohup trim_galore -q 30 --phred33 --length 36 -e 0.1 --stringency 3 --paired -o $dir $fq1 $fq2
done

conda deactivate

echo "Trim Galore completed. Trimmed files are in: $dir"
