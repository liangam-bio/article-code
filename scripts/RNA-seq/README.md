# RNA-seq Analysis

This directory contains scripts for processing mRNA-seq data from Drosophila S2 cells, Kc cells, and sxl mutant flies.

## Pipeline Overview

1. Quality control: FastQC
2. Adapter trimming: Trim Galore
3. Alignment: HISAT2 to Drosophila dm6 genome
4. Sorting & indexing: SAMtools
5. Quantification: featureCounts

## Usage

All scripts should be run in numerical order (01 → 06).
