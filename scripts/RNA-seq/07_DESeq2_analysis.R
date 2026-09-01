
#!/usr/bin/env Rscript
# ============================================================================
# Script: 07_DESeq2_analysis.R
# Purpose: Differential expression analysis using DESeq2
#          For S2 cells, Kc cells, and sxl mutant flies
# Input:   featureCounts output (sxl.vel_count.txt)
# Output:  DEG tables (.xlsx) and Volcano plots (.tiff)
# ============================================================================

# ============================================================================
# 1. Load required libraries
# ============================================================================
library(rtracklayer)
library(tidyverse)
library(DESeq2)
library(openxlsx)
library(ggplot2)

# ============================================================================
# 2. Load genome annotation
# ============================================================================
setwd("/path/to/your/working/directory/")  # <-- MODIFY THIS

genome.anno <- import("dmel-all-r6.44.gtf") %>% as.data.frame()
genome.anno <- genome.anno[, c("gene_id", "gene_symbol", "seqnames", "start", "end", "width", "type")]
genome.anno <- genome.anno[genome.anno$type == "gene", ]
colnames(genome.anno) <- c("GeneID", "gene_symbol", "seqnames", "start", "end", "width", "type")

# ============================================================================
# 3. S2 cells: SXL expression vs control
# ============================================================================
cat("\n=== Processing S2 cells (SXL expression vs control) ===\n")

setwd("/path/to/S2/DEGs/directory/")  # <-- MODIFY THIS

# Read featureCounts output
exp_count_S2 <- read.table("sxl.vel_count.txt", header = TRUE)
rownames(exp_count_S2) <- exp_count_S2$Geneid

# Rename columns (SXL_1-3: treat, VEL_1-3: control)
colnames(exp_count_S2) <- c("Geneid", "Chr", "Start", "End", "Strand", "Length",
                            "SXL_1", "SXL_2", "SXL_3", "VEL_1", "VEL_2", "VEL_3")

# Prepare count matrix
exp_S2 <- exp_count_S2[, c(10:12, 7:9)]  # control first, then treat
exp_S2 <- exp_S2[rowSums(exp_S2) > 0, ]

# Create DESeq2 object
condition <- factor(c(rep("control", 3), rep("treat", 3)), levels = c("control", "treat"))
colData <- data.frame(row.names = colnames(exp_S2), condition)
dds_S2 <- DESeqDataSetFromMatrix(exp_S2, colData, design = ~ condition)
dds_S2 <- DESeq(dds_S2)

# Extract results
result_S2 <- results(dds_S2, contrast = c("condition", "treat", "control"))
result_S2 <- result_S2[order(result_S2$padj), ]  # Sort by FDR

# Annotate with gene symbols and expression values
diff_gene_all_S2 <- as.data.frame(result_S2)
diff_gene_all_S2$GeneID <- rownames(diff_gene_all_S2)
new.diff_gene_all_S2 <- merge(diff_gene_all_S2, genome.anno, by = "GeneID", all.x = TRUE)

exp_S2_df <- exp_S2
exp_S2_df$GeneID <- rownames(exp_S2)
new.diff_gene_all_S2 <- merge(new.diff_gene_all_S2, exp_S2_df, by = "GeneID", all.x = TRUE)

# Define DEGs using FDR (padj) instead of pvalue
diff_gene_up_S2 <- new.diff_gene_all_S2[new.diff_gene_all_S2$log2FoldChange > 1 & 
                                         new.diff_gene_all_S2$padj < 0.05, ]  # CHANGED: pvalue -> padj
diff_gene_down_S2 <- new.diff_gene_all_S2[new.diff_gene_all_S2$log2FoldChange < -1 & 
                                           new.diff_gene_all_S2$padj < 0.05, ]  # CHANGED: pvalue -> padj
diff_gene_S2 <- rbind(diff_gene_up_S2, diff_gene_down_S2)

# Assign group labels
S2_DEG <- new.diff_gene_all_S2
S2_DEG$group <- ifelse(S2_DEG$GeneID %in% diff_gene_up_S2$GeneID, "up",
                       ifelse(S2_DEG$GeneID %in% diff_gene_down_S2$GeneID, "down", "non"))

# Export DEG tables
setwd("/path/to/output/directory/")  # <-- MODIFY THIS

wb <- createWorkbook()
addWorksheet(wb, "new.diff_gene_all_S2")
addWorksheet(wb, "diff_gene_S2")
addWorksheet(wb, "diff_gene_up_S2")
addWorksheet(wb, "diff_gene_down_S2")

writeData(wb, "new.diff_gene_all_S2", new.diff_gene_all_S2, startCol = 1, startRow = 1, rowNames = FALSE)
writeData(wb, "diff_gene_S2", diff_gene_S2, startCol = 1, startRow = 1, rowNames = FALSE)
writeData(wb, "diff_gene_up_S2", diff_gene_up_S2, startCol = 1, startRow = 1, rowNames = FALSE)
writeData(wb, "diff_gene_down_S2", diff_gene_down_S2, startCol = 1, startRow = 1, rowNames = FALSE)
saveWorkbook(wb, file = "DEGs_S2_FDR_0.05.xlsx")  # Renamed to indicate FDR

# ============================================================================
# 4. Generate Volcano plot for S2 cells (with FDR-based jitter)
# ============================================================================
cat("\n=== Generating volcano plot for S2 cells ===\n")

# Jitter extremely small FDR values to improve visualization
fdr_high <- subset(S2_DEG, -log10(S2_DEG$padj) > 100)  # CHANGED: pvalue -> padj
if (nrow(fdr_high) > 0) {
    fdr_jitter <- round(rnorm(nrow(fdr_high), mean = 95, sd = 5), 0)
    fdr_jitter <- 1 * 10^-fdr_jitter
    fdr_high$padj <- fdr_jitter  # CHANGED: pvalue -> padj
}
fdr_low <- subset(S2_DEG, -log10(S2_DEG$padj) <= 100)  # CHANGED: pvalue -> padj
new_S2_DEG <- rbind(fdr_low, fdr_high)

# Export jittered data
wb <- createWorkbook()
addWorksheet(wb, "new_S2_DEG")
addWorksheet(wb, "fdr_high")
addWorksheet(wb, "fdr_low")
writeData(wb, "new_S2_DEG", new_S2_DEG, startCol = 1, startRow = 1, rowNames = FALSE)
writeData(wb, "fdr_high", fdr_high, startCol = 1, startRow = 1, rowNames = FALSE)
writeData(wb, "fdr_low", fdr_low, startCol = 1, startRow = 1, rowNames = FALSE)
saveWorkbook(wb, file = "new_S2_DEGs_FDR_241118.xlsx")

# Generate volcano plot (y-axis uses padj)
tiff("S2_DEGs_volcano_FDR.tiff", units = "in", width = 6.2, height = 5.8, res = 600)
ggplot(data = new_S2_DEG, aes(x = log2FoldChange, y = -log10(padj))) +  # CHANGED: pvalue -> padj
    geom_point(alpha = 0.6, size = 2.5, aes(color = group)) +
    scale_color_manual(values = c("#CD6155", "#CCD1D1", "#5499C7")) +
    geom_vline(xintercept = c(-1, 1), lty = 2, col = "black", lwd = 0.3) +
    geom_hline(yintercept = -log10(0.05), lty = 2, col = "black", lwd = 0.3) +  # FDR cutoff = 0.05
    labs(title = "", x = "log2 Fold Change", y = "-log10 FDR") +  # Changed label
    theme_bw() +
    xlim(-5, 10) + ylim(0, 100) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 18, color = "black"),
          axis.text.y = element_text(size = 18, color = "black"),
          axis.title.x = element_text(size = 21, color = "black"),
          axis.title.y = element_text(size = 21, color = "black"),
          legend.text = element_text(size = 12, color = "black"),
          legend.title = element_text(size = 15, color = "black"))
dev.off()

cat("S2 DEG analysis completed.\n")

# ============================================================================
# 5. Kc cells: SXL knockdown vs control
# ============================================================================
cat("\n=== Processing Kc cells (SXL knockdown vs control) ===\n")

setwd("/path/to/Kc/DEGs/directory/")  # <-- MODIFY THIS

# READ: Add your Kc cell data processing here
# exp_count_Kc <- read.table("kc.vel_count.txt", header = TRUE)
# ... (similar to S2 analysis, but using padj)

cat("Kc DEG analysis completed.\n")

# ============================================================================
# 6. sxl mutant flies: sxl^m/m vs WT
# ============================================================================
cat("\n=== Processing sxl mutant flies (sxl^m/m vs WT) ===\n")

setwd("/path/to/sxl_mutant/DEGs/directory/")  # <-- MODIFY THIS

# READ: Add your sxl mutant data processing here
# exp_count_mutant <- read.table("mutant.vel_count.txt", header = TRUE)
# ... (similar to S2 analysis, but using padj)

cat("sxl mutant DEG analysis completed.\n")

# ============================================================================
# 7. Session information
# ============================================================================
cat("\n=== Session Info ===\n")
sessionInfo()

cat("\n=== All analyses completed successfully! ===\n")
