#!/bin/bash
# ============================================================================
# Script: 06_featurecounts.sh
# Purpose: Quantify gene expression using featureCounts
# Input:   Sorted BAM files
# Output:  Gene count matrix (sxl.vel_count.txt)
# ============================================================================

# Activate conda environment
source activate rnaseq

# GTF annotation file (adjust to your system)
GTF="/data1/amliang/annotation/fly/flyase/gtf/last.v.gtf/dmel-all-r6.44.gtf"

# BAM files directory (adjust to your system)
BAM_DIR="/data1/amliang/projects/sxl-circ/20220418_Sxl.vel_RNA_Seq/align-hisat/"

# Run featureCounts
# -T 16: number of threads
# -p: paired-end
# -t exon: feature type
# -g gene_id: attribute type
# -a: annotation file
# -o: output file
featureCounts -T 16 \
              -p \
              -t exon \
              -g gene_id \
              -a $GTF \
              -o sxl.vel_count.txt \
              ${BAM_DIR}/*.bam

conda deactivate

echo "featureCounts completed. Gene count matrix: sxl.vel_count.txt"
