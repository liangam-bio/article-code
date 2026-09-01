#!/bin/bash
# ============================================================================
# Script: 04_samtools_sort.sh
# Purpose: Convert SAM to BAM and sort BAM files
# Input:   SAM files (.sam)
# Output:  Sorted BAM files (.bam)
# ============================================================================

# Activate conda environment (chipseq environment has samtools)
source activate chipseq

# Convert all SAM files to sorted BAM files
ls *.sam | while read id
do 
    samtools sort -O bam -T sorted -o $(basename ${id} ".sam").bam ${id}
done

conda deactivate

echo "SAM to BAM conversion and sorting completed."
