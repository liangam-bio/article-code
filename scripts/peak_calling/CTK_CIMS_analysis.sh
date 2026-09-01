#!/bin/bash
# ==============================================================================
# Script: 02_CTK_CIMS_analysis.sh
# Purpose: Identify and filter Crosslink-Induced Mutation Sites (CIMS) for SXL
#          PAR-CLIP data using the CLIP Tool Kit (CTK)
# ==============================================================================

# ==============================================================================
# 1. Activate conda environment
# ==============================================================================

source activate chipseq

# ==============================================================================
# 2. Set paths (MODIFY IF NEEDED)
# ==============================================================================

WORK_DIR="/data1/amliang/projects/sxl_DC/20260801_GB_review.revise/4.SXL_CLIP/"
BAM_DIR="${WORK_DIR}/align/"
CIMS_OUT_DIR="${WORK_DIR}/cims_analysis/"
REFERENCE_FA="/data1/amliang/annotation/Drosaphila/omniCLIP/Dro_chrom/dm6.fa"
CHROM_SIZES="/data1/amliang/annotation/Drosaphila/omniCLIP/Dro_chrom/dm6.chrom.sizes"

mkdir -p $CIMS_OUT_DIR

# ==============================================================================
# 3. Define SXL IP BAM files (4 replicates)
# ==============================================================================

SXL_BAMS=(
    "${BAM_DIR}/SXL_rep1-Contain_unmapped_tmDup.bam"
    "${BAM_DIR}/SXL_rep2-Contain_unmapped_tmDup.bam"
    "${BAM_DIR}/SXL_rep3-Contain_unmapped_tmDup.bam"
    "${BAM_DIR}/SXL_rep4-Contain_unmapped_tmDup.bam"
)

# ==============================================================================
# 4. Run CIMS analysis for each replicate
# ==============================================================================

echo "=== Starting CIMS analysis for SXL replicates ==="

for bam_file in "${SXL_BAMS[@]}"; do
    sample_name=$(basename "$bam_file" "-Contain_unmapped_tmDup.bam")
    echo "Processing sample: ${sample_name}"

    ctk-parse-cims \
        --subst \
        --type sub \
        --nuc t \
        --mut c \
        --min-count 2 \
        --min-depth 5 \
        --ref $REFERENCE_FA \
        --bam $bam_file \
        --out-prefix ${CIMS_OUT_DIR}/${sample_name}

    echo "CIMS analysis for ${sample_name} completed."
done

# ==============================================================================
# 5. Merge and filter CIMS sites across replicates
# ==============================================================================

echo "=== Merging and filtering CIMS sites across replicates ==="

cat ${CIMS_OUT_DIR}/*.cims.bed > ${CIMS_OUT_DIR}/all_replicates.cims.raw.bed

bedtools sort -i ${CIMS_OUT_DIR}/all_replicates.cims.raw.bed | \
bedtools merge -c 4 -o count_distinct -i - | \
awk '$4 >= 2 {print $0}' > ${CIMS_OUT_DIR}/SXL_high_confidence_cims.bed

bedtools slop -i ${CIMS_OUT_DIR}/SXL_high_confidence_cims.bed \
              -g $CHROM_SIZES \
              -l 10 -r 10 \
              -s \
              > ${CIMS_OUT_DIR}/SXL_high_confidence_cims_flank10.bed

# ==============================================================================
# 6. Deactivate conda environment
# ==============================================================================

conda deactivate

echo "=== CIMS analysis completed ==="
