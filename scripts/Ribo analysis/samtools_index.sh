#!/bin/bash
# ============================================================================
# Script: samtools_index.sh
# Purpose: Index sorted BAM files (S2 cells)
# Input:   Sorted BAM files (*.sort.bam)
# Output:  BAM index files (*.sort.bam.bai)
# ============================================================================

source activate chipseq

ls *.sort.bam | while read id
do
    echo "Indexing: ${id}"
    samtools index $id
done

conda deactivate

echo "=== BAM indexing completed ==="
echo "Index files: *.sort.bam.bai"
