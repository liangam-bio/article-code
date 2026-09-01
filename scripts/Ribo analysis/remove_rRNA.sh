#!/bin/bash
# ============================================================================
# Script: remove_rRNA.sh
# Purpose: Remove rRNA reads from Ribo-seq data (S2 cells)
# Input:   Trimmed FASTQ files (*_trimmed.fq.gz)
# Output:  rRNA-depleted FASTQ files (*.rm.rRNA.fq.gz)
# ============================================================================

# ============================================================================
# 1. Activate conda environment
# ============================================================================
source activate rnaseq

# ============================================================================
# 2. rRNA reference index (MODIFY IF NEEDED)
# ============================================================================
RRNA_INDEX="/data1/amliang/reference/rRNA/Drosophila.M/rRNA/rRNA"

# ============================================================================
# 3. Remove rRNA reads from all trimmed FASTQ files
# ============================================================================
echo "=== Removing rRNA reads from Ribo-seq data (S2 cells) ==="

ls ./*_trimmed.fq.gz | while read id
do
    # Extract base name (remove _trimmed.fq.gz suffix)
    base=$(basename ${id} "_trimmed.fq.gz")
    
    echo "Processing: ${base}"
    
    # Bowtie2: remove rRNA reads, keep non-rRNA reads
    bowtie2 -p 18 \
            -x $RRNA_INDEX \
            --un-gz ${base}.rm.rRNA.fq.gz \
            -U $id \
            -S ${base}.rRNA.sam
    
    echo "Finished: ${base}"
done

# ============================================================================
# 4. Clean up SAM files (optional - uncomment to save space)
# ============================================================================
# rm *.rRNA.sam

# ============================================================================
# 5. Deactivate conda environment
# ============================================================================
conda deactivate

echo "=== rRNA removal completed ==="
echo "rRNA-depleted files: *.rm.rRNA.fq.gz"
