#!/bin/bash
# ============================================================================
# Script: 07_featurecounts.sh
# Purpose: Quantify Ribo-seq reads mapping to CDS regions using featureCounts
#          S2 cells only: 3 control + 3 treatment samples
# Input:   Sorted BAM files (*.sort.bam)
# Output:  CDS count matrix (Ribo.CDS_count.txt)
# ============================================================================

# ============================================================================
# 1. Activate conda environment
# ============================================================================
source activate rnaseq

# ============================================================================
# 2. Set reference paths (MODIFY IF NEEDED)
# ============================================================================
GTF="/data1/amliang/annotation/fly/Drosophila_melanogaster.BDGP6.22.42.gtf"

# Alternative GTF (commented out):
# GTF="/data1/amliang/annotation/fly/flyase/gtf/last.v.gtf/dmel-all-r6.44.gtf"

# ============================================================================
# 3. Run featureCounts on CDS features
# ============================================================================
echo "=== Quantifying Ribo-seq reads on CDS regions ==="

featureCounts -T 16 \
              -t CDS \
              -g gene_id \
              -a $GTF \
              -o Ribo.CDS_count.txt \
              ./*.sort.bam

conda deactivate

echo "=== featureCounts completed ==="
echo "Output: Ribo.CDS_count.txt"
