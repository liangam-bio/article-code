#!/bin/bash
# ============================================================================
# Script: 03_hisat2_align.sh
# Purpose: Align trimmed reads to Drosophila genome using HISAT2
# Input:   Trimmed FASTQ files: ${sample}_1.clean.fq.gz, ${sample}_2.clean.fq.gz
# Output:  SAM files: ${sample}.hisat.sam
# ============================================================================

# Activate conda environment
source activate rnaseq

# Reference genome index path (adjust to your system)
INDEX="/data1/amliang/reference/index/hisat2/fly_dmel-all-chromosome-r6.31/dmel-all-chromosome-r6.31"

# Sample IDs (adjust as needed)
# S2 cells: S2-S1, S2-S2, S2-S3 (expression) and S2-V1, S2-V2, S2-V3 (control)
# Kc cells: modify sample list accordingly
# sxl mutant flies: modify sample list accordingly
for id in {S2-S1,S2-S2,S2-S3,S2-V1,S2-V2,S2-V3};
do 
    echo "Processing sample ${id}"
    
    # Input paths (adjust based on your directory structure)
    INPUT_DIR="/data1/amliang/projects/sxl-circ/20220418_Sxl.vel_RNA_Seq/rawdata/cleandata/"
    
    # Run HISAT2 alignment
    # -p 18: number of threads
    # -x: reference genome index
    # -1: forward reads
    # -2: reverse reads
    # -S: output SAM file
    hisat2 -p 18 \
           -x $INDEX \
           -1 ${INPUT_DIR}/${id}_1.clean.fq.gz \
           -2 ${INPUT_DIR}/${id}_2.clean.fq.gz \
           -S ./${id}.hisat.sam
done

conda deactivate

echo "HISAT2 alignment completed. SAM files are in the current directory."
