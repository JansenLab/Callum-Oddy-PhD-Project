# ====================================================================
# COMPLETE BULK RNA-SEQ ANALYSIS PIPELINE
# Dual Analysis: Simple Model vs Patient-Adjusted Model
# ====================================================================

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(tibble)
library(EnhancedVolcano)
library(ggrepel)
library(gridExtra)
library(VennDiagram)
library(ggpubr)

setwd('/Users/callu/OneDrive - University College London/Protocol/RNAseq/Claude/New_Results/Attempt_2/')

# Color scheme
color_palette <- c("Gastric" = "#C75A5A", "IM" = "#9370B8", "Duodenal" = "#6B8DB8")
patient_palette <- c("Patient3" = "#E69F00", "Patient4" = "#56B4E9", "Patient6" = "#009E73")

# ====================================================================
# DATA LOADING
# ====================================================================

cat("=== Loading Data ===\n\n")

merged_counts_file <- "/Users/callu/OneDrive - University College London/Protocol/RNAseq/salmon.merged.gene_counts.tsv"
sample_sheet <- "/Users/callu/OneDrive - University College London/Protocol/RNAseq/sample_metadata.csv"

cts <- read.csv(merged_counts_file, sep="\t", header=TRUE, row.names=1)
colData <- read.csv(sample_sheet, row.names=1)

# Prepare count matrix
ctsNew <- cts[, 2:11]
ctsNew <- round(ctsNew)

# Fix sample names if needed (Im13 -> IM13)
if("Im13" %in% rownames(colData)) {
  rownames(colData)[rownames(colData) == "Im13"] <- "IM13"
}
if("Im13" %in% colnames(ctsNew)) {
  colnames(ctsNew)[colnames(ctsNew) == "Im13"] <- "IM13"
}

# Prepare metadata
colData$phenotype <- factor(colData$phenotype, levels = c("Gastric", "IM", "Duodenal"))

# Add patient information
colData$Patient <- factor(case_when(
  grepl("3", rownames(colData)) ~ "Patient3",
  grepl("4", rownames(colData)) ~ "Patient4",
  grepl("6", rownames(colData)) ~ "Patient6"
))

cat("Samples loaded:", ncol(ctsNew), "\n")
cat("Genes before filtering:", nrow(ctsNew), "\n\n")

# ====================================================================
# QUALITY CONTROL (SHARED FOR BOTH MODELS)
# ====================================================================

cat("=== Quality Control Analysis ===\n\n")

# Calculate QC metrics
qc_metrics <- data.frame(
  Sample_ID = colnames(ctsNew),
  Phenotype = colData$phenotype,
  Patient = colData$Patient,
  Total_Reads = colSums(ctsNew),
  Genes_Detected = colSums(ctsNew > 0),
  Percent_Detected = colSums(ctsNew > 0) / nrow(ctsNew) * 100
)

# Summary by phenotype
qc_summary <- qc_metrics %>%
  group_by(Phenotype) %>%
  summarise(
    N = n(),
    Mean_Reads = mean(Total_Reads),
    SD_Reads = sd(Total_Reads),
    Mean_Genes = mean(Genes_Detected),
    .groups = 'drop'
  )

print(qc_summary)
write.csv(qc_summary, "QC_summary_by_phenotype.csv", row.names = FALSE)
write.csv(qc_metrics, "QC_all_samples.csv", row.names = FALSE)

# QC Plots
p1 <- ggplot(qc_metrics, aes(x=reorder(Sample_ID, Total_Reads), y=Total_Reads/1e6, fill=Phenotype)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=Patient), angle=90, hjust=-0.1, size=2.5) +
  scale_fill_manual(values=color_palette) +
  coord_flip() +
  theme_bw(base_size=11) +
  labs(title="Library Size per Sample", y="Millions of Reads", x="")

p2 <- ggplot(qc_metrics, aes(x=Phenotype, y=Total_Reads/1e6, fill=Phenotype)) +
  geom_boxplot(alpha=0.7, outlier.shape=NA) +
  geom_jitter(aes(shape=Patient), width=0.15, size=3) +
  scale_fill_manual(values=color_palette) +
  scale_shape_manual(values=c(15,16,17)) +
  theme_bw(base_size=12) +
  theme(legend.position="right") +
  labs(title="Library Size by Phenotype", y="Millions of Reads")

p3 <- ggplot(qc_metrics, aes(x=Phenotype, y=Genes_Detected, fill=Phenotype)) +
  geom_boxplot(alpha=0.7, outlier.shape=NA) +
  geom_jitter(aes(shape=Patient), width=0.15, size=3) +
  scale_fill_manual(values=color_palette) +
  scale_shape_manual(values=c(15,16,17)) +
  theme_bw(base_size=12) +
  theme(legend.position="right") +
  labs(title="Genes Detected by Phenotype", y="Number of Genes")

pdf("QC_sequencing_metrics.pdf", width=14, height=10)
grid.arrange(p1, arrangeGrob(p2, p3, ncol=2), ncol=1, heights=c(2,1))
dev.off()

cat("QC plots saved.\n\n")

# ====================================================================
# FUNCTION: RUN COMPLETE DESEQ2 ANALYSIS
# ====================================================================

run_complete_analysis <- function(design_formula, model_name, output_prefix) {
  
  cat(paste0("\n========================================\n"))
  cat(paste0("RUNNING: ", model_name, "\n"))
  cat(paste0("Design: ", deparse(design_formula), "\n"))
  cat(paste0("========================================\n\n"))
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = ctsNew,
                                colData = colData,
                                design = design_formula)
  
  # Filter low count genes
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  cat("Genes after filtering:", nrow(dds), "\n\n")
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Get normalized counts
  norm_counts <- counts(dds, normalized=TRUE)
  write.csv(norm_counts, paste0(output_prefix, "_normalized_counts.csv"))
  
  # VST transformation for visualization
  vsd <- vst(dds, blind=TRUE)
  
  # ====================================================================
  # GLOBAL STRUCTURE: PCA
  # ====================================================================
  
  pcaData <- plotPCA(vsd, intgroup="phenotype", returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pcaData$Patient <- colData$Patient
  
  pca_plot <- ggplot(pcaData, aes(PC1, PC2, color=group, shape=Patient)) +
    geom_point(size=6, alpha=0.8) +
    geom_text_repel(aes(label=name), size=3.5, max.overlaps=20) +
    scale_color_manual(values=color_palette, name="Phenotype") +
    scale_shape_manual(values=c(15,16,17)) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_bw(base_size=14) +
    ggtitle(paste("PCA:", model_name)) +
    theme(plot.title = element_text(hjust=0.5, face="bold"))
  
  ggsave(paste0(output_prefix, "_PCA.pdf"), pca_plot, width=11, height=8)
  
  # ====================================================================
  # GLOBAL STRUCTURE: Sample Distance Heatmap
  # ====================================================================
  
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  
  annotation_col <- data.frame(
    Phenotype = colData$phenotype,
    Patient = colData$Patient
  )
  rownames(annotation_col) <- colnames(vsd)
  
  annotation_colors <- list(
    Phenotype = color_palette,
    Patient = patient_palette
  )
  
  pdf(paste0(output_prefix, "_sample_distance_heatmap.pdf"), width=10, height=9)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(255),
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           main = paste("Sample Distance:", model_name))
  dev.off()
  
  # ====================================================================
  # GLOBAL STRUCTURE: Correlation Heatmap
  # ====================================================================
  
  cor_matrix <- cor(assay(vsd), method="pearson")
  
  pdf(paste0(output_prefix, "_correlation_heatmap.pdf"), width=10, height=9)
  pheatmap(cor_matrix,
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           col = colorRampPalette(c("blue", "white", "red"))(100),
           breaks = seq(0.85, 1, length.out=101),
           display_numbers = TRUE,
           number_format = "%.2f",
           fontsize_number = 7,
           main = paste("Sample Correlation:", model_name))
  dev.off()
  
  # ====================================================================
  # DIFFERENTIAL EXPRESSION: All Comparisons
  # ====================================================================
  
  cat("\n=== Running Differential Expression Analysis ===\n")
  
  # IM vs Gastric
  res_IM_Gastric <- results(dds, contrast=c("phenotype", "IM", "Gastric"))
  summary_IM_Gastric <- summary(res_IM_Gastric)
  cat("\nIM vs Gastric:\n")
  print(summary_IM_Gastric)
  
  # Duodenal vs Gastric
  res_Duo_Gastric <- results(dds, contrast=c("phenotype", "Duodenal", "Gastric"))
  summary_Duo_Gastric <- summary(res_Duo_Gastric)
  cat("\nDuodenal vs Gastric:\n")
  print(summary_Duo_Gastric)
  
  # Duodenal vs IM
  res_Duo_IM <- results(dds, contrast=c("phenotype", "Duodenal", "IM"))
  summary_Duo_IM <- summary(res_Duo_IM)
  cat("\nDuodenal vs IM:\n")
  print(summary_Duo_IM)
  
  # Save results
  write.csv(as.data.frame(res_IM_Gastric), 
            paste0(output_prefix, "_DEG_IM_vs_Gastric.csv"))
  write.csv(as.data.frame(res_Duo_Gastric), 
            paste0(output_prefix, "_DEG_Duodenal_vs_Gastric.csv"))
  write.csv(as.data.frame(res_Duo_IM), 
            paste0(output_prefix, "_DEG_Duodenal_vs_IM.csv"))
  
  # ====================================================================
  # DEG SUMMARY STATISTICS
  # ====================================================================
  
  get_deg_counts <- function(res) {
    data.frame(
      Total_DEG = sum(res$padj < 0.1 & abs(res$log2FoldChange) > 1, na.rm=TRUE),
      Upregulated = sum(res$padj < 0.1 & res$log2FoldChange > 1, na.rm=TRUE),
      Downregulated = sum(res$padj < 0.1 & res$log2FoldChange < -1, na.rm=TRUE)
    )
  }
  
  deg_summary <- bind_rows(
    get_deg_counts(res_IM_Gastric) %>% mutate(Comparison = "IM vs Gastric"),
    get_deg_counts(res_Duo_Gastric) %>% mutate(Comparison = "Duodenal vs Gastric"),
    get_deg_counts(res_Duo_IM) %>% mutate(Comparison = "Duodenal vs IM")
  ) %>% select(Comparison, everything())
  
  write.csv(deg_summary, paste0(output_prefix, "_DEG_summary.csv"), row.names=FALSE)
  print(deg_summary)
  
  # ====================================================================
  # MA PLOTS (with padj < 0.1 threshold)
  # ====================================================================
  
  pdf(paste0(output_prefix, "_MA_plots.pdf"), width=15, height=5)
  par(mfrow=c(1,3))
  plotMA(res_IM_Gastric, ylim=c(-5,5), main="IM vs Gastric",
         colNonSig="grey80", colSig=color_palette["IM"], alpha=0.1)
  abline(h=c(-1,1), col="blue", lty=2)
  plotMA(res_Duo_Gastric, ylim=c(-5,5), main="Duodenal vs Gastric",
         colNonSig="grey80", colSig=color_palette["Duodenal"], alpha=0.1)
  abline(h=c(-1,1), col="blue", lty=2)
  plotMA(res_Duo_IM, ylim=c(-5,5), main="Duodenal vs IM",
         colNonSig="grey80", colSig="#E69F00", alpha=0.1)
  abline(h=c(-1,1), col="blue", lty=2)
  dev.off()
  
  # ====================================================================
  # VOLCANO PLOTS
  # ====================================================================
  
  create_volcano <- function(res,
                             title,
                             key_color,
                             top_n = 30,
                             padj_cutoff = 0.1,
                             lfc_cutoff  = 1,
                             extreme_lfc = 15) {
    
    # Convert to dataframe
    res_full <- as.data.frame(res) %>%
      rownames_to_column("gene")
    
    # Top N significant genes by padj
    top_sig <- res_full %>%
      filter(!is.na(padj),
             !is.na(log2FoldChange),
             padj < padj_cutoff,
             abs(log2FoldChange) > lfc_cutoff) %>%
      arrange(padj) %>%
      head(top_n) %>%
      pull(gene)
    
    # Extreme log2FC genes
    extreme_genes <- res_full %>%
      filter(!is.na(log2FoldChange),
             log2FoldChange < -extreme_lfc | log2FoldChange > extreme_lfc) %>%
      pull(gene)
    
    # Combine unique labels
    label_genes <- unique(c(top_sig, extreme_genes))
    
    # Colour vector
    keyvals <- rep("grey70", nrow(res_full))
    names(keyvals) <- rep("NS", nrow(res_full))
    
    sig_up <- which(!is.na(res_full$padj) &
                      !is.na(res_full$log2FoldChange) &
                      res_full$padj < padj_cutoff &
                      res_full$log2FoldChange > lfc_cutoff)
    
    sig_down <- which(!is.na(res_full$padj) &
                        !is.na(res_full$log2FoldChange) &
                        res_full$padj < padj_cutoff &
                        res_full$log2FoldChange < -lfc_cutoff)
    
    keyvals[sig_up] <- key_color
    names(keyvals)[sig_up] <- "Up"
    
    keyvals[sig_down] <- "#4169E1"
    names(keyvals)[sig_down] <- "Down"
    
    EnhancedVolcano(
      res,
      lab = rownames(res),
      x = "log2FoldChange",
      y = "padj",
      selectLab = label_genes,
      title = title,
      subtitle = paste0(
        "Top ", length(top_sig), " DEGs (padj < ", padj_cutoff, 
        ", |LFC| > ", lfc_cutoff, "); plus extreme LFC > ±", extreme_lfc
      ),
      pCutoff = padj_cutoff,
      FCcutoff = lfc_cutoff,
      pointSize = 2.5,
      labSize = 4.2,
      colCustom = keyvals,
      colAlpha = 0.6,
      legendPosition = "right",
      drawConnectors = TRUE,
      widthConnectors = 0.4,
      max.overlaps = 20
    )
  }
  
  
  
  pdf(paste0(output_prefix, "_volcano_IM_vs_Gastric.pdf"), width=12, height=10)
  print(create_volcano(res_IM_Gastric, "IM vs Gastric", color_palette["IM"]))
  dev.off()
  
  pdf(paste0(output_prefix, "_volcano_Duodenal_vs_Gastric.pdf"), width=12, height=10)
  print(create_volcano(res_Duo_Gastric, "Duodenal vs Gastric", color_palette["Duodenal"]))
  dev.off()
  
  pdf(paste0(output_prefix, "_volcano_Duodenal_vs_IM.pdf"), width=12, height=10)
  print(create_volcano(res_Duo_IM, "Duodenal vs IM", "#E69F00"))
  dev.off()
  
  # ====================================================================
  # VENN DIAGRAM
  # ====================================================================
  
  get_sig_genes <- function(res, padj_cutoff=0.1, lfc_cutoff=1) {
    sig <- res[!is.na(res$padj) & res$padj < padj_cutoff & 
                 abs(res$log2FoldChange) > lfc_cutoff, ]
    return(rownames(sig))
  }
  
  deg_IM_Gastric <- get_sig_genes(res_IM_Gastric)
  deg_Duo_Gastric <- get_sig_genes(res_Duo_Gastric)
  deg_Duo_IM <- get_sig_genes(res_Duo_IM)
  
  venn.diagram(
    x = list(
      "IM vs Gastric" = deg_IM_Gastric,
      "Duo vs Gastric" = deg_Duo_Gastric,
      "Duo vs IM" = deg_Duo_IM
    ),
    filename = paste0(output_prefix, "_venn_diagram.png"),
    col = "black",
    fill = c(color_palette["IM"], color_palette["Duodenal"], "#E69F00"),
    alpha = 0.5,
    cex = 1.5,
    cat.cex = 1.3,
    cat.fontface = "bold",
    margin = 0.1,
    imagetype = "png"
  )
  
  # ====================================================================
  # TOP VARIABLE GENES HEATMAP
  # ====================================================================
  
  topVarGenes <- head(order(rowVars(assay(vsd)), decreasing=TRUE), 500)
  mat <- assay(vsd)[topVarGenes, ]
  mat <- t(scale(t(mat)))
  
  pdf(paste0(output_prefix, "_top500_variable_genes_heatmap.pdf"), width=10, height=12)
  pheatmap(mat,
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           show_rownames = FALSE,
           cluster_cols = TRUE,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           breaks = seq(-2, 2, length.out=101),
           main = paste("Top 500 Variable Genes:", model_name))
  dev.off()
  
  # ====================================================================
  # BAR PLOTS OF TOP DEGs
  # ====================================================================
  
  plot_top_degs <- function(res, title, n=20) {
    res_df <- as.data.frame(res) %>%
      rownames_to_column("Gene") %>%
      filter(!is.na(padj), padj < 0.1, abs(log2FoldChange) > 1) %>%
      arrange(log2FoldChange)
    
    if(nrow(res_df) == 0) {
      cat("No significant DEGs for", title, "\n")
      return(NULL)
    }
    
    top_genes <- rbind(head(res_df, min(n/2, nrow(res_df))), 
                       tail(res_df, min(n/2, nrow(res_df))))
    top_genes$Direction <- ifelse(top_genes$log2FoldChange > 0, "Up", "Down")
    
    ggplot(top_genes, aes(x=log2FoldChange, y=reorder(Gene, log2FoldChange), fill=Direction)) +
      geom_bar(stat="identity") +
      scale_fill_manual(values=c("Down"="#4169E1", "Up"="#C75A5A")) +
      theme_bw(base_size=11) +
      labs(title=title, x="log2 Fold Change", y="") +
      theme(plot.title = element_text(hjust=0.5, face="bold"),
            axis.text.y = element_text(size=8))
  }
  
  pdf(paste0(output_prefix, "_top_DEGs_barplots.pdf"), width=10, height=14)
  p1 <- plot_top_degs(res_IM_Gastric, "IM vs Gastric: Top DEGs")
  p2 <- plot_top_degs(res_Duo_Gastric, "Duodenal vs Gastric: Top DEGs")
  p3 <- plot_top_degs(res_Duo_IM, "Duodenal vs IM: Top DEGs")
  plots <- list(p1, p2, p3)
  plots <- plots[!sapply(plots, is.null)]
  if(length(plots) > 0) {
    do.call(grid.arrange, c(plots, ncol=1))
  }
  dev.off()
  
  # Return results for comparison
  return(list(
    dds = dds,
    vsd = vsd,
    res_IM_Gastric = res_IM_Gastric,
    res_Duo_Gastric = res_Duo_Gastric,
    res_Duo_IM = res_Duo_IM,
    deg_summary = deg_summary,
    deg_lists = list(
      IM_Gastric = deg_IM_Gastric,
      Duo_Gastric = deg_Duo_Gastric,
      Duo_IM = deg_Duo_IM
    )
  ))
}

# ====================================================================
# RUN BOTH MODELS
# ====================================================================

# Model 1: Simple design (phenotype only)
results_simple <- run_complete_analysis(
  design_formula = ~ phenotype,
  model_name = "Simple Model (Phenotype Only)",
  output_prefix = "Simple"
)

# Model 2: Patient-adjusted design
results_paired <- run_complete_analysis(
  design_formula = ~ Patient + phenotype,
  model_name = "Patient-Adjusted Model",
  output_prefix = "Paired"
)

# ====================================================================
# COMPARISON BETWEEN MODELS
# ====================================================================

cat("\n========================================\n")
cat("COMPARING SIMPLE VS PATIENT-ADJUSTED MODELS\n")
cat("========================================\n\n")

# Function to compare two result sets
compare_results <- function(res1, res2, comparison_name) {
  
  comp_df <- data.frame(
    gene = rownames(res1),
    lfc_simple = res1$log2FoldChange,
    padj_simple = res1$padj,
    lfc_paired = res2$log2FoldChange,
    padj_paired = res2$padj
  ) %>%
    filter(!is.na(padj_simple), !is.na(padj_paired)) %>%
    mutate(
      sig_simple = padj_simple < 0.1 & abs(lfc_simple) > 1,
      sig_paired = padj_paired < 0.1 & abs(lfc_paired) > 1,
      Status = case_when(
        sig_simple & sig_paired ~ "Both Significant",
        sig_simple & !sig_paired ~ "Lost in Paired",
        !sig_simple & sig_paired ~ "Gained in Paired",
        TRUE ~ "Not Significant"
      )
    )
  
  # Summary
  status_summary <- comp_df %>%
    group_by(Status) %>%
    summarise(Count = n(), .groups='drop')
  
  print(comparison_name)
  print(status_summary)
  write.csv(status_summary, 
            paste0("Comparison_", gsub(" ", "_", comparison_name), "_summary.csv"),
            row.names=FALSE)
  
  # Scatter plot: Log2FC
  p1 <- ggplot(comp_df, aes(x=lfc_simple, y=lfc_paired, color=Status)) +
    geom_point(alpha=0.5, size=1.5) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
    scale_color_manual(values=c("Both Significant"="darkgreen",
                                "Lost in Paired"="orange",
                                "Gained in Paired"="purple",
                                "Not Significant"="grey70")) +
    theme_bw(base_size=14) +
    labs(title=paste("LFC Comparison:", comparison_name),
         x="log2FC (Simple Model)",
         y="log2FC (Patient-Adjusted Model)") +
    theme(plot.title = element_text(hjust=0.5, face="bold"))
  
  # Scatter plot: P-values
  p2 <- ggplot(comp_df, aes(x=-log10(padj_simple), y=-log10(padj_paired), color=Status)) +
    geom_point(alpha=0.5, size=1.5) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
    geom_hline(yintercept=-log10(0.05), linetype="dotted") +
    geom_vline(xintercept=-log10(0.05), linetype="dotted") +
    scale_color_manual(values=c("Both Significant"="darkgreen",
                                "Lost in Paired"="orange",
                                "Gained in Paired"="purple",
                                "Not Significant"="grey70")) +
    theme_bw(base_size=14) +
    labs(title=paste("Significance Comparison:", comparison_name),
         x="-log10(padj) Simple",
         y="-log10(padj) Paired") +
    theme(plot.title = element_text(hjust=0.5, face="bold"))
  
  return(list(p1=p1, p2=p2, data=comp_df, summary=status_summary))
}

# Compare all three contrasts
comp_IM <- compare_results(results_simple$res_IM_Gastric, 
                           results_paired$res_IM_Gastric,
                           "IM vs Gastric")

comp_Duo <- compare_results(results_simple$res_Duo_Gastric,
                            results_paired$res_Duo_Gastric,
                            "Duodenal vs Gastric")

comp_DuoIM <- compare_results(results_simple$res_Duo_IM,
                              results_paired$res_Duo_IM,
                              "Duodenal vs IM")

# Save comparison plots
pdf("Comparison_IM_vs_Gastric.pdf", width=14, height=7)
grid.arrange(comp_IM$p1, comp_IM$p2, ncol=2)
dev.off()

pdf("Comparison_Duodenal_vs_Gastric.pdf", width=14, height=7)
grid.arrange(comp_Duo$p1, comp_Duo$p2, ncol=2)
dev.off()

pdf("Comparison_Duodenal_vs_IM.pdf", width=14, height=7)
grid.arrange(comp_DuoIM$p1, comp_DuoIM$p2, ncol=2)
dev.off()

# ====================================================================
# DEG COUNT COMPARISON BAR PLOT
# ====================================================================

deg_comparison <- bind_rows(
  results_simple$deg_summary %>% mutate(Model = "Simple"),
  results_paired$deg_summary %>% mutate(Model = "Patient-Adjusted")
)

p_deg_comp <- ggplot(deg_comparison, aes(x=Comparison, y=Total_DEG, fill=Model)) +
  geom_bar(stat="identity", position="dodge") +
  geom_text(aes(label=Total_DEG), position=position_dodge(width=0.9), vjust=-0.5) +
  scale_fill_manual(values=c("Simple"="grey60", "Patient-Adjusted"="#4A90E2")) +
  theme_bw(base_size=14) +
  labs(title="DEG Count Comparison: Simple vs Patient-Adjusted Models",
       y="Number of DEGs (padj<0.1, |LFC|>1)") +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        axis.text.x = element_text(angle=45, hjust=1))

ggsave("Comparison_DEG_counts.pdf", p_deg_comp, width=10, height=7)

# ====================================================================
# FINAL SUMMARY
# ====================================================================

cat("\n========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat("SIMPLE MODEL FILES:\n")
cat("  - Simple_*.pdf/csv\n\n")

cat("PATIENT-ADJUSTED MODEL FILES:\n")
cat("  - Paired_*.pdf/csv\n\n")

cat("COMPARISON FILES:\n")
cat("  - Comparison_*.pdf/csv\n\n")

cat("QC FILES (shared):\n")
cat("  - QC_*.pdf/csv\n\n")

# Print final DEG summary
cat("\nFINAL DEG COUNTS:\n")
print(deg_comparison)

cat("\nAll analyses complete. Files saved to working directory.\n")





# ====================================================================
# FUNCTION: EXTRACT VENN SEGMENTS AS TABLES
# ====================================================================

extract_venn_segments <- function(set_A, set_B, set_C,
                                  name_A = "A",
                                  name_B = "B",
                                  name_C = "C") {
  
  all_genes <- unique(c(set_A, set_B, set_C))
  
  tibble(
    gene = all_genes,
    in_A = gene %in% set_A,
    in_B = gene %in% set_B,
    in_C = gene %in% set_C
  ) %>%
    mutate(
      Segment = case_when(
        in_A & !in_B & !in_C ~ paste0(name_A, "_only"),
        !in_A & in_B & !in_C ~ paste0(name_B, "_only"),
        !in_A & !in_B & in_C ~ paste0(name_C, "_only"),
        in_A & in_B & !in_C ~ paste0(name_A, "_and_", name_B),
        in_A & !in_B & in_C ~ paste0(name_A, "_and_", name_C),
        !in_A & in_B & in_C ~ paste0(name_B, "_and_", name_C),
        in_A & in_B & in_C ~ paste0(name_A, "_and_", name_B, "_and_", name_C)
      )
    ) %>%
    select(gene, Segment)
}

# ====================================================================
# VENN SEGMENT TABLES — SIMPLE MODEL
# ====================================================================

venn_simple_tbl <- extract_venn_segments(
  set_A = results_simple$deg_lists$IM_Gastric,
  set_B = results_simple$deg_lists$Duo_Gastric,
  set_C = results_simple$deg_lists$Duo_IM,
  name_A = "IM_vs_Gastric",
  name_B = "Duo_vs_Gastric",
  name_C = "Duo_vs_IM"
)

# Save full table
write.csv(venn_simple_tbl,
          "Simple_Venn_Segments_All_Genes.csv",
          row.names = FALSE)

# Save one file per segment (nice for inspection)
venn_simple_tbl %>%
  group_by(Segment) %>%
  group_walk(~ write.csv(.x,
                         paste0("Simple_Venn_", .y$Segment, ".csv"),
                         row.names = FALSE))



# ====================================================================
# VENN SEGMENT TABLES — PATIENT-ADJUSTED MODEL
# ====================================================================

venn_paired_tbl <- extract_venn_segments(
  set_A = results_paired$deg_lists$IM_Gastric,
  set_B = results_paired$deg_lists$Duo_Gastric,
  set_C = results_paired$deg_lists$Duo_IM,
  name_A = "IM_vs_Gastric",
  name_B = "Duo_vs_Gastric",
  name_C = "Duo_vs_IM"
)

write.csv(venn_paired_tbl,
          "Paired_Venn_Segments_All_Genes.csv",
          row.names = FALSE)

venn_paired_tbl %>%
  group_by(Segment) %>%
  group_walk(~ write.csv(.x,
                         paste0("Paired_Venn_", .y$Segment, ".csv"),
                         row.names = FALSE))


# ====================================================================
# SIMPLE vs PAIRED DEG OVERLAP (IM vs Gastric example)
# ====================================================================

deg_overlap_IM <- tibble(
  gene = union(results_simple$deg_lists$IM_Gastric,
               results_paired$deg_lists$IM_Gastric),
  Simple = gene %in% results_simple$deg_lists$IM_Gastric,
  Paired = gene %in% results_paired$deg_lists$IM_Gastric
) %>%
  mutate(
    Status = case_when(
      Simple & Paired ~ "Shared",
      Simple & !Paired ~ "Simple_only",
      !Simple & Paired ~ "Paired_only"
    )
  )

write.csv(deg_overlap_IM,
          "IM_vs_Gastric_Simple_vs_Paired_DEG_overlap.csv",
          row.names = FALSE)

res_IM_Gas <- read.csv("Paired_DEG_IM_vs_Gastric.csv", row.names=1)
























# ====================================================================
# GENE SET ENRICHMENT ANALYSIS (GSEA)
# GO and KEGG Analysis for Upregulated and Downregulated Genes
# ====================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(tibble)

setwd('/Users/callu/OneDrive - University College London/Protocol/RNAseq/Claude/New_Results/Attempt_2/')

# ====================================================================
# LOAD DEG RESULTS (Choose Simple or Paired model)
# ====================================================================

# Change this to "Simple" or "Paired" depending on which model you want to analyze
MODEL <- "Paired"

cat(paste0("Running enrichment analysis for: ", MODEL, " model\n\n"))

# Load the DEG results
res_IM_Gastric <- read.csv(paste0(MODEL, "_DEG_IM_vs_Gastric.csv"), row.names=1)
res_Duo_Gastric <- read.csv(paste0(MODEL, "_DEG_Duodenal_vs_Gastric.csv"), row.names=1)
res_Duo_IM <- read.csv(paste0(MODEL, "_DEG_Duodenal_vs_IM.csv"), row.names=1)

# ====================================================================
# FUNCTION: EXTRACT GENE LISTS
# ====================================================================

extract_gene_lists <- function(res, padj_cutoff = 0.1, lfc_cutoff = 1) {
  
  # Filter significant genes
  sig_genes <- res %>%
    filter(!is.na(padj), padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff)
  
  # Upregulated genes
  up_genes <- sig_genes %>%
    filter(log2FoldChange > lfc_cutoff) %>%
    rownames()
  
  # Downregulated genes
  down_genes <- sig_genes %>%
    filter(log2FoldChange < -lfc_cutoff) %>%
    rownames()
  
  return(list(
    upregulated = up_genes,
    downregulated = down_genes,
    all_sig = rownames(sig_genes)
  ))
}

# ====================================================================
# EXTRACT GENE LISTS FOR ALL COMPARISONS
# ====================================================================

genes_IM_Gastric <- extract_gene_lists(res_IM_Gastric)
genes_Duo_Gastric <- extract_gene_lists(res_Duo_Gastric)
genes_Duo_IM <- extract_gene_lists(res_Duo_IM)

cat("Gene counts extracted:\n")
cat("IM vs Gastric - Up:", length(genes_IM_Gastric$upregulated), 
    "Down:", length(genes_IM_Gastric$downregulated), "\n")
cat("Duodenal vs Gastric - Up:", length(genes_Duo_Gastric$upregulated), 
    "Down:", length(genes_Duo_Gastric$downregulated), "\n")
cat("Duodenal vs IM - Up:", length(genes_Duo_IM$upregulated), 
    "Down:", length(genes_Duo_IM$downregulated), "\n\n")

# ====================================================================
# FUNCTION: CONVERT GENE SYMBOLS TO ENTREZ IDs
# ====================================================================

convert_to_entrez <- function(gene_symbols) {
  entrez_ids <- mapIds(org.Hs.eg.db,
                       keys = gene_symbols,
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first")
  
  # Remove NAs
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  
  cat("Converted", length(gene_symbols), "symbols to", length(entrez_ids), "Entrez IDs\n")
  
  return(entrez_ids)
}

# ====================================================================
# FUNCTION: RUN GO ENRICHMENT
# ====================================================================

run_go_enrichment <- function(genes, comparison_name, direction, ont = "BP") {
  
  cat(paste0("\nRunning GO (", ont, ") enrichment: ", comparison_name, " - ", direction, "\n"))
  
  if(length(genes) == 0) {
    cat("No genes provided. Skipping.\n")
    return(NULL)
  }
  
  entrez_ids <- convert_to_entrez(genes)
  
  if(length(entrez_ids) < 5) {
    cat("Too few genes after conversion. Skipping.\n")
    return(NULL)
  }
  
  ego <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  ont = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)
  
  if(is.null(ego) || nrow(ego@result) == 0) {
    cat("No significant enrichment found.\n")
    return(NULL)
  }
  
  cat("Found", nrow(ego@result), "significant GO terms\n")
  
  return(ego)
}

# ====================================================================
# FUNCTION: RUN KEGG ENRICHMENT
# ====================================================================

run_kegg_enrichment <- function(genes, comparison_name, direction) {
  
  cat(paste0("\nRunning KEGG enrichment: ", comparison_name, " - ", direction, "\n"))
  
  if(length(genes) == 0) {
    cat("No genes provided. Skipping.\n")
    return(NULL)
  }
  
  entrez_ids <- convert_to_entrez(genes)
  
  if(length(entrez_ids) < 5) {
    cat("Too few genes after conversion. Skipping.\n")
    return(NULL)
  }
  
  ekegg <- enrichKEGG(gene = entrez_ids,
                      organism = 'hsa',
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)
  
  if(is.null(ekegg) || nrow(ekegg@result) == 0) {
    cat("No significant enrichment found.\n")
    return(NULL)
  }
  
  cat("Found", nrow(ekegg@result), "significant KEGG pathways\n")
  
  return(ekegg)
}

# ====================================================================
# FUNCTION: CREATE ENRICHMENT PLOT
# ====================================================================

create_enrichment_plot <- function(enrich_result, title, max_terms = 20, color = "#C75A5A") {
  
  if(is.null(enrich_result)) {
    return(NULL)
  }
  
  # Get top terms and order by gene ratio for nice cascading effect
  plot_data <- enrich_result@result %>%
    arrange(p.adjust) %>%
    head(max_terms) %>%
    mutate(
      GeneRatio_numeric = sapply(GeneRatio, function(x) {
        nums <- as.numeric(unlist(strsplit(as.character(x), "/")))
        nums[1] / nums[2]
      })
    ) %>%
    arrange(GeneRatio_numeric) %>%  # Order by gene ratio (lowest to highest)
    mutate(Description = factor(Description, levels = Description))
  
  p <- ggplot(plot_data, aes(x = GeneRatio_numeric, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "#C75A5A", high = "#6B8DB8", name = "Adjusted\np-value") +
    scale_size_continuous(name = "Gene\nCount", range = c(3, 10)) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right"
    ) +
    labs(
      title = title,
      x = "Gene Ratio",
      y = ""
    )
  
  return(p)
}

# ====================================================================
# RUN ALL GO ENRICHMENT ANALYSES (Biological Process)
# ====================================================================

cat("\n========================================\n")
cat("GENE ONTOLOGY (BP) ENRICHMENT ANALYSIS\n")
cat("========================================\n")

# IM vs Gastric
go_IM_Gas_up <- run_go_enrichment(genes_IM_Gastric$upregulated, 
                                  "IM vs Gastric", "Upregulated", "BP")
go_IM_Gas_down <- run_go_enrichment(genes_IM_Gastric$downregulated, 
                                    "IM vs Gastric", "Downregulated", "BP")

# Duodenal vs Gastric
go_Duo_Gas_up <- run_go_enrichment(genes_Duo_Gastric$upregulated, 
                                   "Duodenal vs Gastric", "Upregulated", "BP")
go_Duo_Gas_down <- run_go_enrichment(genes_Duo_Gastric$downregulated, 
                                     "Duodenal vs Gastric", "Downregulated", "BP")

# Duodenal vs IM
go_Duo_IM_up <- run_go_enrichment(genes_Duo_IM$upregulated, 
                                  "Duodenal vs IM", "Upregulated", "BP")
go_Duo_IM_down <- run_go_enrichment(genes_Duo_IM$downregulated, 
                                    "Duodenal vs IM", "Downregulated", "BP")

# ====================================================================
# RUN ALL KEGG ENRICHMENT ANALYSES
# ====================================================================

cat("\n========================================\n")
cat("KEGG PATHWAY ENRICHMENT ANALYSIS\n")
cat("========================================\n")

# IM vs Gastric
kegg_IM_Gas_up <- run_kegg_enrichment(genes_IM_Gastric$upregulated, 
                                      "IM vs Gastric", "Upregulated")
kegg_IM_Gas_down <- run_kegg_enrichment(genes_IM_Gastric$downregulated, 
                                        "IM vs Gastric", "Downregulated")

# Duodenal vs Gastric
kegg_Duo_Gas_up <- run_kegg_enrichment(genes_Duo_Gastric$upregulated, 
                                       "Duodenal vs Gastric", "Upregulated")
kegg_Duo_Gas_down <- run_kegg_enrichment(genes_Duo_Gastric$downregulated, 
                                         "Duodenal vs Gastric", "Downregulated")

# Duodenal vs IM
kegg_Duo_IM_up <- run_kegg_enrichment(genes_Duo_IM$upregulated, 
                                      "Duodenal vs IM", "Upregulated")
kegg_Duo_IM_down <- run_kegg_enrichment(genes_Duo_IM$downregulated, 
                                        "Duodenal vs IM", "Downregulated")

# ====================================================================
# SAVE ENRICHMENT RESULTS AS CSV
# ====================================================================

save_enrichment_csv <- function(enrich_result, filename) {
  if(!is.null(enrich_result) && nrow(enrich_result@result) > 0) {
    write.csv(enrich_result@result, filename, row.names = FALSE)
    cat("Saved:", filename, "\n")
  }
}

cat("\nSaving enrichment results as CSV...\n")

# GO results
save_enrichment_csv(go_IM_Gas_up, paste0(MODEL, "_GO_IM_vs_Gastric_UP.csv"))
save_enrichment_csv(go_IM_Gas_down, paste0(MODEL, "_GO_IM_vs_Gastric_DOWN.csv"))
save_enrichment_csv(go_Duo_Gas_up, paste0(MODEL, "_GO_Duodenal_vs_Gastric_UP.csv"))
save_enrichment_csv(go_Duo_Gas_down, paste0(MODEL, "_GO_Duodenal_vs_Gastric_DOWN.csv"))
save_enrichment_csv(go_Duo_IM_up, paste0(MODEL, "_GO_Duodenal_vs_IM_UP.csv"))
save_enrichment_csv(go_Duo_IM_down, paste0(MODEL, "_GO_Duodenal_vs_IM_DOWN.csv"))

# KEGG results
save_enrichment_csv(kegg_IM_Gas_up, paste0(MODEL, "_KEGG_IM_vs_Gastric_UP.csv"))
save_enrichment_csv(kegg_IM_Gas_down, paste0(MODEL, "_KEGG_IM_vs_Gastric_DOWN.csv"))
save_enrichment_csv(kegg_Duo_Gas_up, paste0(MODEL, "_KEGG_Duodenal_vs_Gastric_UP.csv"))
save_enrichment_csv(kegg_Duo_Gas_down, paste0(MODEL, "_KEGG_Duodenal_vs_Gastric_DOWN.csv"))
save_enrichment_csv(kegg_Duo_IM_up, paste0(MODEL, "_KEGG_Duodenal_vs_IM_UP.csv"))
save_enrichment_csv(kegg_Duo_IM_down, paste0(MODEL, "_KEGG_Duodenal_vs_IM_DOWN.csv"))

# ====================================================================
# CREATE ALL 12 PLOTS
# ====================================================================

cat("\n========================================\n")
cat("CREATING ENRICHMENT PLOTS\n")
cat("========================================\n")

# Define colors
up_color <- "#C75A5A"    # Red for upregulated
down_color <- "#4169E1"  # Blue for downregulated

# GO Plots
cat("\nCreating GO plots...\n")

if(!is.null(go_IM_Gas_up)) {
  p <- create_enrichment_plot(go_IM_Gas_up, 
                              "GO: IM vs Gastric (Upregulated)", 
                              color = up_color)
  ggsave(paste0(MODEL, "_GO_IM_vs_Gastric_UP.pdf"), p, width = 12, height = 8)
  cat("Saved: GO IM vs Gastric UP\n")
}

if(!is.null(go_IM_Gas_down)) {
  p <- create_enrichment_plot(go_IM_Gas_down, 
                              "GO: IM vs Gastric (Downregulated)", 
                              color = down_color)
  ggsave(paste0(MODEL, "_GO_IM_vs_Gastric_DOWN.pdf"), p, width = 12, height = 8)
  cat("Saved: GO IM vs Gastric DOWN\n")
}

if(!is.null(go_Duo_Gas_up)) {
  p <- create_enrichment_plot(go_Duo_Gas_up, 
                              "GO: Duodenal vs Gastric (Upregulated)", 
                              color = up_color)
  ggsave(paste0(MODEL, "_GO_Duodenal_vs_Gastric_UP.pdf"), p, width = 12, height = 8)
  cat("Saved: GO Duodenal vs Gastric UP\n")
}

if(!is.null(go_Duo_Gas_down)) {
  p <- create_enrichment_plot(go_Duo_Gas_down, 
                              "GO: Duodenal vs Gastric (Downregulated)", 
                              color = down_color)
  ggsave(paste0(MODEL, "_GO_Duodenal_vs_Gastric_DOWN.pdf"), p, width = 12, height = 8)
  cat("Saved: GO Duodenal vs Gastric DOWN\n")
}

if(!is.null(go_Duo_IM_up)) {
  p <- create_enrichment_plot(go_Duo_IM_up, 
                              "GO: Duodenal vs IM (Upregulated)", 
                              color = up_color)
  ggsave(paste0(MODEL, "_GO_Duodenal_vs_IM_UP.pdf"), p, width = 12, height = 8)
  cat("Saved: GO Duodenal vs IM UP\n")
}

if(!is.null(go_Duo_IM_down)) {
  p <- create_enrichment_plot(go_Duo_IM_down, 
                              "GO: Duodenal vs IM (Downregulated)", 
                              color = down_color)
  ggsave(paste0(MODEL, "_GO_Duodenal_vs_IM_DOWN.pdf"), p, width = 12, height = 8)
  cat("Saved: GO Duodenal vs IM DOWN\n")
}

# KEGG Plots
cat("\nCreating KEGG plots...\n")

if(!is.null(kegg_IM_Gas_up)) {
  p <- create_enrichment_plot(kegg_IM_Gas_up, 
                              "KEGG: IM vs Gastric (Upregulated)", 
                              color = up_color)
  ggsave(paste0(MODEL, "_KEGG_IM_vs_Gastric_UP.pdf"), p, width = 12, height = 8)
  cat("Saved: KEGG IM vs Gastric UP\n")
}

if(!is.null(kegg_IM_Gas_down)) {
  p <- create_enrichment_plot(kegg_IM_Gas_down, 
                              "KEGG: IM vs Gastric (Downregulated)", 
                              color = down_color)
  ggsave(paste0(MODEL, "_KEGG_IM_vs_Gastric_DOWN.pdf"), p, width = 12, height = 8)
  cat("Saved: KEGG IM vs Gastric DOWN\n")
}

if(!is.null(kegg_Duo_Gas_up)) {
  p <- create_enrichment_plot(kegg_Duo_Gas_up, 
                              "KEGG: Duodenal vs Gastric (Upregulated)", 
                              color = up_color)
  ggsave(paste0(MODEL, "_KEGG_Duodenal_vs_Gastric_UP.pdf"), p, width = 12, height = 8)
  cat("Saved: KEGG Duodenal vs Gastric UP\n")
}

if(!is.null(kegg_Duo_Gas_down)) {
  p <- create_enrichment_plot(kegg_Duo_Gas_down, 
                              "KEGG: Duodenal vs Gastric (Downregulated)", 
                              color = down_color)
  ggsave(paste0(MODEL, "_KEGG_Duodenal_vs_Gastric_DOWN.pdf"), p, width = 12, height = 8)
  cat("Saved: KEGG Duodenal vs Gastric DOWN\n")
}

if(!is.null(kegg_Duo_IM_up)) {
  p <- create_enrichment_plot(kegg_Duo_IM_up, 
                              "KEGG: Duodenal vs IM (Upregulated)", 
                              color = up_color)
  ggsave(paste0(MODEL, "_KEGG_Duodenal_vs_IM_UP.pdf"), p, width = 12, height = 8)
  cat("Saved: KEGG Duodenal vs IM UP\n")
}

if(!is.null(kegg_Duo_IM_down)) {
  p <- create_enrichment_plot(kegg_Duo_IM_down, 
                              "KEGG: Duodenal vs IM (Downregulated)", 
                              color = down_color)
  ggsave(paste0(MODEL, "_KEGG_Duodenal_vs_IM_DOWN.pdf"), p, width = 12, height = 8)
  cat("Saved: KEGG Duodenal vs IM DOWN\n")
}

# ====================================================================
# SUMMARY
# ====================================================================

cat("\n========================================\n")
cat("ENRICHMENT ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat("Generated files:\n")
cat("  - 12 PDF plots (GO and KEGG for up/down in all comparisons)\n")
cat("  - 12 CSV files with complete enrichment results\n\n")

cat("Model analyzed:", MODEL, "\n\n")

cat("To analyze the other model, change MODEL <- \"Paired\" to MODEL <- \"Simple\" at the top of the script.\n")
