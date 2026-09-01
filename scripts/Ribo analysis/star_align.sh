#!/bin/bash
# ============================================================================
# Script: star_align.sh
# Purpose: Align rRNA-depleted Ribo-seq reads to Drosophila genome using STAR
#          S2 cells only: 3 control + 3 treatment samples
# Input:   rRNA-depleted FASTQ files (*.rm.rRNA.fq.gz)
# Output:  BAM files (*.bam), gene counts, and unmapped reads
# ============================================================================

# ============================================================================
# 1. Activate conda environment
# ============================================================================
source activate rnaseq

# ============================================================================
# 2. Set reference paths (MODIFY IF NEEDED)
# ============================================================================
Ref_index="/data1/amliang/reference/index/star/flybase_m6/"
GTF="/data1/amliang/annotation/fly/flyase/gtf/dmel-all-r6.31.gtf"
Read_dir="./"
Run_log="runlog.txt"

# ============================================================================
# 3. S2 cell samples: 3 control + 3 treatment
# ============================================================================
# Modify sample names to match your files:
# Control: CTRL_1, CTRL_2, CTRL_3
# Treatment: SXL_1, SXL_2, SXL_3
for id in {CTRL_1,CTRL_2,CTRL_3,SXL_1,SXL_2,SXL_3};
do
    echo "Processing: ${id}"
    
    # Input: rRNA-depleted FASTQ
    Read=${Read_dir}/${id}.rm.rRNA.fq.gz
    
    # Run STAR alignment
    STAR \
        --runThreadN 10 \
        --genomeDir $Ref_index \
        --readFilesIn $Read \
        --readFilesCommand zcat \
        --sjdbGTFfile $GTF \
        --outFileNamePrefix ${id}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outBAMcompression 6 \
        --outSAMattributes All \
        --outSAMattrRGline ID:1 LB:ribo_seq PL:ILLUMINA SM:${id} \
        --outReadsUnmapped Fastx \
        --quantMode TranscriptomeSAM GeneCounts \
        --outFilterType BySJout \
        --outFilterMismatchNmax 2 \
        --outFilterMultimapNmax 1 \
        --limitBAMsortRAM 85799345920 \
        2>> $Run_log
    
    echo "Finished: ${id}"
done

conda deactivate

echo "=== STAR alignment completed ==="
echo "Output files: ${id}_Aligned.sortedByCoord.out.bam for each sample"
