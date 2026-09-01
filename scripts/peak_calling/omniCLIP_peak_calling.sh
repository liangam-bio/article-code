#!/bin/bash
# ==============================================================================
# Script: 01_omniCLIP_peak_calling.sh
# Purpose: SXL CLIP-seq peak calling using omniCLIP
#          S2 cells: SXL IP vs IgG control
# Input:   BAM files (SXL IP and IgG control)
# Output:  omniCLIP peak calls
# ==============================================================================


# ==============================================================================
# 1. Activate conda environment
# ==============================================================================

source activate omniCLIP


# ==============================================================================
# 2. Set reference paths (MODIFY IF NEEDED)
# ==============================================================================

DB_FILE="/data1/amliang/annotation/Drosaphila/omniCLIP/omniCLIP_db/Drosophila_melanogaster.BDGP6.54.63.db"
GENOME_DIR="/data1/amliang/annotation/Drosaphila/omniCLIP/Dro_chrom/"


# ==============================================================================
# 3. Set working directory (MODIFY THIS PATH)
# ==============================================================================

WORK_DIR="/data1/amliang/projects/sxl_DC/20260801_GB_review.revise/4.SXL_CLIP/"
BAM_DIR="${WORK_DIR}/align/"
PEAK_DIR="${WORK_DIR}/peak/"


# ==============================================================================
# 4. Create output directory
# ==============================================================================

mkdir -p $PEAK_DIR


# ==============================================================================
# 5. Generate background (IgG control) data
# ==============================================================================

echo "=== Generating background data (IgG control) ==="

omniCLIP parsingBG \
    --db-file $DB_FILE \
    --genome-dir $GENOME_DIR \
    --bg-files ${BAM_DIR}/IgG_rep1-Contain_unmapped_tmDup.bam \
    --bg-files ${BAM_DIR}/IgG_rep2-Contain_unmapped_tmDup.bam \
    --bg-files ${BAM_DIR}/IgG_rep3-Contain_unmapped_tmDup.bam \
    --out-file ${PEAK_DIR}/SXL_IgG_background.dat

echo "Background data: ${PEAK_DIR}/SXL_IgG_background.dat"


# ==============================================================================
# 6. Generate CLIP (SXL IP) data
# ==============================================================================

echo "=== Generating CLIP data (SXL IP) ==="

omniCLIP parsingCLIP \
    --db-file $DB_FILE \
    --genome-dir $GENOME_DIR \
    --clip-files ${BAM_DIR}/SXL_rep1-Contain_unmapped_tmDup.bam \
    --clip-files ${BAM_DIR}/SXL_rep2-Contain_unmapped_tmDup.bam \
    --clip-files ${BAM_DIR}/SXL_rep3-Contain_unmapped_tmDup.bam \
    --clip-files ${BAM_DIR}/SXL_rep4-Contain_unmapped_tmDup.bam \
    --out-file ${PEAK_DIR}/SXL_IP.dat

echo "CLIP data: ${PEAK_DIR}/SXL_IP.dat"


# ==============================================================================
# 7. Run omniCLIP peak calling
# ==============================================================================

echo "=== Running omniCLIP peak calling ==="

omniCLIP run_omniCLIP \
    --db-file $DB_FILE \
    --bg-dat ${PEAK_DIR}/SXL_IgG_background.dat \
    --clip-dat ${PEAK_DIR}/SXL_IP.dat \
    --out-dir ${PEAK_DIR}/SXL_peaks

echo "Peak calling completed."
echo "Results: ${PEAK_DIR}/SXL_peaks/"


# ==============================================================================
# 8. Deactivate conda environment
# ==============================================================================

conda deactivate


# ==============================================================================
# 9. Summary
# ==============================================================================

echo ""
echo "=== omniCLIP SXL peak calling completed ==="
echo "  Input SXL IP BAMs: 4 replicates"
echo "  Input IgG BAMs:    3 replicates"
echo "  Output:            ${PEAK_DIR}/SXL_peaks/"
