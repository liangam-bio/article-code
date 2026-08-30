# article-code
Code for SXL_porject
This repository contains all custom scripts used for the analysis in SXL-project.

All analyses were performed on a Linux-based computing cluster. The following software and R packages are required:

### Software
- **HISAT2** (v2.1.0) — for mRNA-seq alignment
- **STAR** (v2.6.1) — for alternative splicing and Ribo-seq alignment
- **featureCounts** (v1.6.4) — for gene expression quantification
- **Bowtie2** (v2.3.5.1) — for ChIP-seq alignment
- **MACS2** (v2.2.6) — for ChIP-seq peak calling
- **HOMER** — for de novo motif enrichment analysis
- **BEDTools** (v2.31.0) — for genomic interval operations
- **UMI-tools** (v1.1.5) — for deduplication of CLIP-seq reads
- **CTK** (v1.1.3) — for CLIP-seq crosslink-induced mutation site (CIMS) identification
- **MEME** (v5.0.5) — for de novo motif discovery

### R Packages
- **DESeq2** (v1.42.0) — differential expression analysis
- **rMATS** (v4.1.2) — alternative splicing analysis
- **clusterProfiler** (v4.10.0) — GO enrichment analysis
- **ChIPseeker** (v1.42.0) — ChIP peak annotation
- **ggplot2** — data visualization
- **randomForest** — machine learning for SXL-regulated TA gene prediction

### Python Packages
- **numpy**, **pandas**, **matplotlib** — data processing and visualization (if applicable)

## Data Availability
Sequencing data generated in this study have been submitted to GEO under accession numbers: GSE329406, GSE329409, and GSE329649.

