#!/bin/bash
# ============================================================================
# Script: 05_samtools_sort.sh
# Purpose: Sort STAR output BAM files (S2 cells)
# Input:   STAR BAM files (*Aligned.sortedByCoord.out.bam)
# Output:  Sorted BAM files (*.sort.bam)
# ============================================================================

source activate chipseq

ls *.out.bam | while read id
do
    echo "Sorting: ${id}"
    samtools sort -O bam -T sorted -o $(basename ${id} ".out.bam").sort.bam ${id}
done

conda deactivate

echo "=== BAM sorting completed ==="
echo "Sorted files: *.sort.bam"
