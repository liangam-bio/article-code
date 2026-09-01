#!/bin/bash
# ============================================================================
# Script: 05_samtools_index.sh
# Purpose: Index sorted BAM files for downstream analysis
# Input:   Sorted BAM files (.bam)
# Output:  BAM index files (.bam.bai)
# ============================================================================

# Activate conda environment
source activate chipseq

# Index all BAM files in current directory
ls *.bam | while read id
do 
    samtools index $id
done

conda deactivate

echo "BAM indexing completed. .bam.bai files generated."
