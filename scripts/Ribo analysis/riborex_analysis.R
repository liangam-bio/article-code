#!/usr/bin/env Rscript
# ==============================================================================
# RiboREX Differential Translation Analysis
# ==============================================================================
# Script   : 08_riborex_analysis.R
# Purpose  : Identify genes with altered translation efficiency (TE) upon SXL
#            overexpression in Drosophila S2 cells
# Method   : RiboREX (DESeq2-based) comparing Ribo-seq vs. cytoplasmic RNA-seq
# ==============================================================================
# KEY CONCEPT
# ==============================================================================
# Translation Efficiency (TE) = Ribo-seq / Cytoplasmic RNA-seq
#
# Cytoplasmic RNA-seq (not total RNA) is used as the input pool because it
# represents the translatable mRNA fraction available to ribosomes, providing
# an accurate baseline for measuring translational regulation.
# ==============================================================================


# ==============================================================================
# 1. Environment
# ==============================================================================

rm(list = ls())

setwd("E:/220625_PC/R workplace/220320_SXL/240115_tranlation/RiboSeq/featureCount_CDS_multi/")


# ==============================================================================
# 2. Load Ribo-seq counts (ribosome-protected fragments)
# ==============================================================================

ribo_raw <- read.table("Ribo.CDS.multi_count.a.txt", header = TRUE)
rownames(ribo_raw) <- ribo_raw$Geneid

colnames(ribo_raw) <- c(
    "Geneid", "Chr", "Start", "End", "Strand", "Length",
    "WT_1", "WT_2", "WT_3", "SXL_1", "SXL_2", "SXL_3"
)

ribo_counts <- ribo_raw[, paste0(rep(c("WT", "SXL"), each = 3), "_", 1:3)]
ribo_counts <- ribo_counts[rowSums(ribo_counts) > 0, ]

message(sprintf("Ribo-seq: %d genes retained", nrow(ribo_counts)))


# ==============================================================================
# 3. Load cytoplasmic RNA-seq counts (translatable mRNA pool)
# ==============================================================================
# Critical: These are cytoplasmic fraction counts, NOT total RNA.
#           This ensures TE reflects true translational activity.
# ==============================================================================

setwd("E:/220625_PC/R workplace/220320_SXL/240115_tranlation/transla_by_Cyt.RNAseq/")

rna_raw <- read.table("RNASeq.CDS_count.txt", header = TRUE)
rownames(rna_raw) <- rna_raw$Geneid

colnames(rna_raw) <- c(
    "Geneid", "Chr", "Start", "End", "Strand", "Length",
    "C.Sxl_1", "C.Sxl_2", "C.Sxl_3",
    "C.WT_1",  "C.WT_2",  "C.WT_3",
    "N.Sxl_1", "N.Sxl_2", "N.Sxl_3",
    "N.WT_1",  "N.WT_2",  "N.WT_3",
    "T.Sxl_1", "T.Sxl_2", "T.Sxl_3",
    "T.WT_1",  "T.WT_2",  "T.WT_3"
)

# Extract cytoplasmic fractions (C.) only
rna_cyto <- rna_raw[, c("C.WT_1", "C.WT_2", "C.WT_3", "C.Sxl_1", "C.Sxl_2", "C.Sxl_3")]
colnames(rna_cyto) <- c("WT_1", "WT_2", "WT_3", "SXL_1", "SXL_2", "SXL_3")

# Keep only genes also present in Ribo-seq
rna_cyto <- rna_cyto[rownames(rna_cyto) %in% rownames(ribo_counts), ]
rna_cyto <- rna_cyto[rowSums(rna_cyto) > 0, ]

message(sprintf("Cytoplasmic RNA-seq: %d genes retained", nrow(rna_cyto)))


# ==============================================================================
# 4. Merge and clean count matrices
# ==============================================================================

common_genes <- intersect(rownames(ribo_counts), rownames(rna_cyto))
message(sprintf("Common genes: %d", length(common_genes)))

ribo_final <- ribo_counts[common_genes, ] %>% na.omit()
rna_final  <- rna_cyto[common_genes, ]  %>% na.omit()


# ==============================================================================
# 5. Run RiboREX
# ==============================================================================

library(fdrtool)
library(riborex)
library(DESeq2)

# Define experimental conditions
conditions <- c(rep("control", 3), rep("treated", 3))

# Run RiboREX (DESeq2 engine)
res <- riborex(
    rna_counts  = rna_final,
    ribo_counts = ribo_final,
    rna_cond    = conditions,
    ribo_cond   = conditions,
    engine      = "DESeq2"
)

# Quality check: p-value distribution
hist(
    res$pvalue,
    main     = "RiboREX: Unadjusted p-value Distribution",
    xlab     = "p-value",
    col      = "#F4C7AB",
    border   = NA,
    cex.main = 1.2
)


# ==============================================================================
# 6. Define differentially translated genes
# ==============================================================================
# Thresholds:
#   |log2FC| > 1  : biologically meaningful change in TE
#   padj < 0.05   : statistically significant (FDR-corrected)
#   baseMean > 10 : sufficient read support
# ==============================================================================

res_df <- as.data.frame(res)

te_sig <- res_df %>%
    filter(abs(log2FoldChange) > 1,
           padj < 0.05,
           baseMean > 10) %>%
    na.omit()

te_up   <- filter(te_sig, log2FoldChange > 0)
te_down <- filter(te_sig, log2FoldChange < 0)

message(sprintf("Differentially translated genes: %d", nrow(te_sig)))
message(sprintf("  ↑ TE (SXL-induced): %d", nrow(te_up)))
message(sprintf("  ↓ TE (SXL-repressed): %d", nrow(te_down)))


# ==============================================================================
# 7. Load genome annotation
# ==============================================================================

setwd("E:/220625_PC/R workplace/220910_annotation/Drosaphila/")

library(rtracklayer)
library(dplyr)

gtf <- import("dmel-all-r6.44.gtf") %>%
    as.data.frame() %>%
    filter(type == "gene") %>%
    select(gene_id, gene_symbol, seqnames, start, end, width, strand) %>%
    setNames(c("GeneID", "gene_symbol", "seqnames", "start", "end", "length", "strand"))


# ==============================================================================
# 8. Annotate full results
# ==============================================================================

res_df$GeneID <- rownames(res_df)

te_annotated <- res_df %>%
    left_join(gtf, by = "GeneID") %>%
    left_join(ribo_final %>% mutate(GeneID = rownames(.)), by = "GeneID") %>%
    left_join(rna_final  %>% mutate(GeneID = rownames(.)), by = "GeneID")


# ==============================================================================
# 9. Secondary filtering (less stringent, for exploratory purposes)
# ==============================================================================
# |log2FC| > 0.6, pvalue < 0.05 (uncorrected)
# ==============================================================================

te_extra <- te_annotated %>%
    filter(abs(log2FoldChange) > 0.6,
           pvalue < 0.05,
           !is.na(GeneID))

te_extra_up   <- filter(te_extra, log2FoldChange > 0)
te_extra_down <- filter(te_extra, log2FoldChange < 0)

message(sprintf("Exploratory filtering (|log2FC| > 0.6, p < 0.05): %d genes", nrow(te_extra)))


# ==============================================================================
# 10. Export results (Excel)
# ==============================================================================

setwd("E:/220625_PC/R workplace/220320_SXL/202404_Fig/240509")

library(openxlsx)

wb <- createWorkbook()
addWorksheet(wb, "All_genes")
addWorksheet(wb, "TE_sig")
addWorksheet(wb, "TE_up")
addWorksheet(wb, "TE_down")
addWorksheet(wb, "Extra_filter")
addWorksheet(wb, "Extra_up")
addWorksheet(wb, "Extra_down")

writeData(wb, "All_genes", te_annotated, rowNames = FALSE)
writeData(wb, "TE_sig",    te_sig,       rowNames = FALSE)
writeData(wb, "TE_up",     te_up,        rowNames = FALSE)
writeData(wb, "TE_down",   te_down,      rowNames = FALSE)
writeData(wb, "Extra_filter", te_extra,  rowNames = FALSE)
writeData(wb, "Extra_up",  te_extra_up,  rowNames = FALSE)
writeData(wb, "Extra_down", te_extra_down, rowNames = FALSE)

saveWorkbook(wb, file = "Translation_S2_240506.xlsx", overwrite = TRUE)

message("Results saved to: Translation_S2_240506.xlsx")


# ==============================================================================
# 11. Session info (for reproducibility)
# ==============================================================================

sessionInfo()

cat("\n=== RiboREX analysis completed successfully ===\n")
cat(sprintf("  Up-regulated TE: %d\n", nrow(te_up)))
cat(sprintf("  Down-regulated TE: %d\n", nrow(te_down)))
cat(sprintf("  Total significant: %d\n", nrow(te_sig)))
