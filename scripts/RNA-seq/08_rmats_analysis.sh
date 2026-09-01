#!/bin/bash
# ============================================================================
# Script: 01_rmats_analysis.sh
# Purpose: Detect differential alternative splicing events using rMATS
#          For S2 cells, Kc cells, and sxl mutant flies
# Input:   BAM files (from RNA-seq alignment)
# Output:  rMATS output directory with AS event tables
# ============================================================================

# ============================================================================
# 1. Activate conda environment
# ============================================================================
source activate alter_spl

# ============================================================================
# 2. Set variables (MODIFY THESE PATHS)
# ============================================================================

# GTF annotation file
GTF="/path/to/your/annotation/new_dmel-all-r6.44.gtf"  # <-- MODIFY THIS

# Output directory
OUTPUT_DIR="./rmats_output"
TMP_DIR="./rmats_tmp"

# Read length (match your sequencing data)
READ_LENGTH=150

# Number of threads
THREADS=1

# Create output and tmp directories if they don't exist
mkdir -p $OUTPUT_DIR
mkdir -p $TMP_DIR

# ============================================================================
# 3. S2 cells: SXL expression vs control
# ============================================================================
echo "=== Running rMATS for S2 cells (SXL expression vs control) ==="

# Create BAM lists for S2 cells
# Control samples (VEL_1, VEL_2, VEL_3)
echo "/path/to/S2/VEL_1.bam" > b1_S2.txt
echo "/path/to/S2/VEL_2.bam" >> b1_S2.txt
echo "/path/to/S2/VEL_3.bam" >> b1_S2.txt

# Treatment samples (SXL_1, SXL_2, SXL_3)
echo "/path/to/S2/SXL_1.bam" > b2_S2.txt
echo "/path/to/S2/SXL_2.bam" >> b2_S2.txt
echo "/path/to/S2/SXL_3.bam" >> b2_S2.txt

# Run rMATS for S2 cells
rmats.py \
    --b1 b1_S2.txt \
    --b2 b2_S2.txt \
    --gtf $GTF \
    --od ${OUTPUT_DIR}/S2 \
    --tmp $TMP_DIR \
    -t paired \
    --novelSS \
    --readLength $READ_LENGTH \
    --cstat 0.05 \
    --nthread $THREADS

echo "S2 rMATS analysis completed."

# ============================================================================
# 4. Kc cells: SXL knockdown vs control
# ============================================================================
echo "=== Running rMATS for Kc cells (SXL knockdown vs control) ==="

# Create BAM lists for Kc cells
# Control samples
echo "/path/to/Kc/control_1.bam" > b1_Kc.txt
echo "/path/to/Kc/control_2.bam" >> b1_Kc.txt
echo "/path/to/Kc/control_3.bam" >> b1_Kc.txt

# Knockdown samples
echo "/path/to/Kc/KD_1.bam" > b2_Kc.txt
echo "/path/to/Kc/KD_2.bam" >> b2_Kc.txt
echo "/path/to/Kc/KD_3.bam" >> b2_Kc.txt

# Run rMATS for Kc cells
rmats.py \
    --b1 b1_Kc.txt \
    --b2 b2_Kc.txt \
    --gtf $GTF \
    --od ${OUTPUT_DIR}/Kc \
    --tmp $TMP_DIR \
    -t paired \
    --novelSS \
    --readLength $READ_LENGTH \
    --cstat 0.05 \
    --nthread $THREADS

echo "Kc rMATS analysis completed."

# ============================================================================
# 5. sxl mutant flies: sxl^m/m vs WT
# ============================================================================
echo "=== Running rMATS for sxl mutant flies (sxl^m/m vs WT) ==="

# Create BAM lists for sxl mutant
# WT samples
echo "/path/to/mutant/WT_1.bam" > b1_mutant.txt
echo "/path/to/mutant/WT_2.bam" >> b1_mutant.txt
echo "/path/to/mutant/WT_3.bam" >> b1_mutant.txt

# Mutant samples (sxl^m/m)
echo "/path/to/mutant/mut_1.bam" > b2_mutant.txt
echo "/path/to/mutant/mut_2.bam" >> b2_mutant.txt
echo "/path/to/mutant/mut_3.bam" >> b2_mutant.txt

# Run rMATS for sxl mutant
rmats.py \
    --b1 b1_mutant.txt \
    --b2 b2_mutant.txt \
    --gtf $GTF \
    --od ${OUTPUT_DIR}/mutant \
    --tmp $TMP_DIR \
    -t paired \
    --novelSS \
    --readLength $READ_LENGTH \
    --cstat 0.05 \
    --nthread $THREADS

echo "sxl mutant rMATS analysis completed."

# ============================================================================
# 6. Clean up temporary files (optional)
# ============================================================================
# rm -rf $TMP_DIR
# rm b1_*.txt b2_*.txt

# ============================================================================
# 7. Deactivate conda environment
# ============================================================================
conda deactivate

echo "=== All rMATS analyses completed successfully! ==="
echo "Results are in: $OUTPUT_DIR"
