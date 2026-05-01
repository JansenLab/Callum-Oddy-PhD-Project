workdir= "setdirectory"
setwd(workdir)

library(ChAMP)
library(knitr)
library(dplyr)
library(tidyr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(minfi)
library(shinyMethyl)
library(ggplot2)
library(DNAmArray)
library(Gviz)
library(DMRcate)
library(stringr)
library(IlluminaHumanMethylationEPICv2manifest)
library("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
library("ggsignif")
library("gplots")
library("DNAmArray")
library("sva")
library("DNAcopy")
library("impute")
library("wateRmelon")
library("ggfortify")
library("irlba")
library("devtools")
library("fastcluster")
library("ggpubr")
library("qqman")

####ChAMP Load in of data, this section will provide QC plots for both pre and post normalization
Sample_Data <- read.csv("setdirectory/IDATs/All Data/Sample_Data.csv")
myImport <- champ.import("setdirectory/IDATs/All Data", arraytype = "EPICv2")
myLoad <-  champ.filter(beta = myImport$beta,
                        pd=myImport$pd,
                        detP=myImport$detP,
                        beadcount=myImport$beadcount,
                        ProbeCutoff=0.1,
                        arraytype = "EPICv2")
myNorm <- champ.norm(myLoad$beta, method = "BMIQ", arraytype = "EPICv2")


############################################
# CORRECTED: Differential Methylation Analysis
############################################

library(limma)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

############################################
# STEP 1: Prepare sample data and design
############################################

# Use your Sample_Data (not myLoad$pd)
pd <- Sample_Data

# Ensure factors are set correctly
pd$Sample_Type <- factor(
  pd$Sample_Type,
  levels = c("Gastric", "IM", "Duodenal")
)
pd$Patient <- factor(pd$Patient)

# Verify
print(head(pd))
print(table(pd$Sample_Type))
print(table(pd$Patient))

############################################
# STEP 2: Convert beta to M-values
############################################

# Clamp beta values to avoid log(0) or log(Inf)
beta <- pmin(pmax(myNorm, 1e-6), 1 - 1e-6)

# Convert to M-values
Mvals <- log2(beta / (1 - beta))

# Check for any NA or Inf values
sum(is.na(Mvals))
sum(is.infinite(Mvals))

# If you have any NA/Inf, remove those probes
if (sum(is.na(Mvals)) > 0 | sum(is.infinite(Mvals)) > 0) {
  bad_probes <- apply(Mvals, 1, function(x) any(is.na(x) | is.infinite(x)))
  cat("Removing", sum(bad_probes), "probes with NA/Inf values\n")
  Mvals <- Mvals[!bad_probes, ]
}

cat("M-values matrix dimensions:", dim(Mvals), "\n")

############################################
# STEP 3: Get and align annotation
############################################

annEPICv2 <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

# Match probes
match_idx <- match(rownames(Mvals), annEPICv2$Name)
cat("Unmatched probes:", sum(is.na(match_idx)), "\n")

# Subset to matched probes BEFORE modeling
Mvals_sub <- Mvals[!is.na(match_idx), ]
annEPICv2Sub <- annEPICv2[match_idx[!is.na(match_idx)], ]

# Verify alignment
stopifnot(all(rownames(Mvals_sub) == annEPICv2Sub$Name))
cat("Final M-values matrix for modeling:", dim(Mvals_sub), "\n")

############################################
# STEP 4: Create design matrix
############################################

design <- model.matrix(~ Sample_Type, data = pd)
colnames(design)
cat("\nDesign matrix:\n")
print(design)

############################################
# STEP 5: Estimate within-patient correlation
############################################

cat("\nEstimating within-patient correlation...\n")
corfit <- duplicateCorrelation(
  Mvals_sub,
  design,
  block = pd$Patient
)

cat("Consensus correlation:", corfit$consensus, "\n")

# Check if correlation is reasonable
if (corfit$consensus < 0) {
  warning("Negative correlation detected. Patient blocking may not be necessary.")
} else if (corfit$consensus > 0.9) {
  warning("Very high correlation (>0.9). Check for potential issues.")
} else {
  cat("Correlation looks reasonable (0.3-0.8 is typical)\n")
}

############################################
# STEP 6: Fit linear model
############################################

cat("\nFitting linear model...\n")
fit <- lmFit(
  Mvals_sub,
  design,
  block = pd$Patient,
  correlation = corfit$consensus
)

fit <- eBayes(fit)

############################################
# STEP 7: Extract annotated DMPs
############################################

cat("\nExtracting DMPs for all comparisons...\n")

# 1. IM vs Gastric
dmp_IM_vs_Gas <- topTable(
  fit,
  coef = "Sample_TypeIM",
  number = Inf,
  adjust.method = "BH",
  genelist = annEPICv2Sub
)

# 2. Duodenal vs Gastric
dmp_Duo_vs_Gas <- topTable(
  fit,
  coef = "Sample_TypeDuodenal",
  number = Inf,
  adjust.method = "BH",
  genelist = annEPICv2Sub
)

# 3. IM vs Duodenal (via contrast)
contrast.matrix <- makeContrasts(
  IM_vs_Duo = Sample_TypeIM - Sample_TypeDuodenal,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

dmp_IM_vs_Duo <- topTable(
  fit2,
  coef = "IM_vs_Duo",
  number = Inf,
  adjust.method = "BH",
  genelist = annEPICv2Sub
)

############################################
# STEP 8: Verify results
############################################

cat("\n=== VERIFICATION ===\n")

# Check dimensions
cat("\nDimensions:\n")
cat("IM vs Gastric:", dim(dmp_IM_vs_Gas), "\n")
cat("Duodenal vs Gastric:", dim(dmp_Duo_vs_Gas), "\n")
cat("IM vs Duodenal:", dim(dmp_IM_vs_Duo), "\n")

# Check significant DMPs (FDR < 0.05)
cat("\nSignificant DMPs (adj.P < 0.05):\n")
cat("IM vs Gastric:", sum(dmp_IM_vs_Gas$adj.P.Val < 0.05), "\n")
cat("Duodenal vs Gastric:", sum(dmp_Duo_vs_Gas$adj.P.Val < 0.05), "\n")
cat("IM vs Duodenal:", sum(dmp_IM_vs_Duo$adj.P.Val < 0.05), "\n")

# Check significant DMPs with effect size (adj.P < 0.05 & |logFC| > 0.2)
cat("\nSignificant DMPs (adj.P < 0.05 & |logFC| > 0.2):\n")
cat("IM vs Gastric:", sum(dmp_IM_vs_Gas$adj.P.Val < 0.05 & abs(dmp_IM_vs_Gas$logFC) > 0.2), "\n")
cat("Duodenal vs Gastric:", sum(dmp_Duo_vs_Gas$adj.P.Val < 0.05 & abs(dmp_Duo_vs_Gas$logFC) > 0.2), "\n")
cat("IM vs Duodenal:", sum(dmp_IM_vs_Duo$adj.P.Val < 0.05 & abs(dmp_IM_vs_Duo$logFC) > 0.2), "\n")

# Check annotation columns
cat("\nAnnotation columns in dmp_IM_vs_Gas:\n")
print(colnames(dmp_IM_vs_Gas))

# View top results
cat("\nTop 10 DMPs (IM vs Gastric):\n")
print(head(dmp_IM_vs_Gas, 10))

############################################
# STEP 9: Save results (without Delta_Beta)
############################################

cat("\nSaving initial results...\n")
write.csv(dmp_IM_vs_Gas, "DMPs_IM_vs_Gastric_initial.csv", row.names = FALSE)
write.csv(dmp_Duo_vs_Gas, "DMPs_Duodenal_vs_Gastric_initial.csv", row.names = FALSE)
write.csv(dmp_IM_vs_Duo, "DMPs_IM_vs_Duodenal_initial.csv", row.names = FALSE)


############################################
# STEP 7b: Calculate Beta Differences (Delta Beta)
############################################

# 1. Get the original Beta values for the matched probes
# (Use the subset we created in Step 3 so rows match the fit object)
Beta_sub <- myNorm[rownames(Mvals_sub), ]

# 2. Calculate Group Means (using Beta values)
# Note: simple rowMeans is usually sufficient for interpretation
mean_Gastric <- rowMeans(Beta_sub[, pd$Sample_Type == "Gastric"], na.rm=TRUE)
mean_IM      <- rowMeans(Beta_sub[, pd$Sample_Type == "IM"], na.rm=TRUE)
mean_Duo     <- rowMeans(Beta_sub[, pd$Sample_Type == "Duodenal"], na.rm=TRUE)

# 3. Add these to your topTables

# --- Process IM vs Gastric ---
dmp_IM_vs_Gas$Mean_Gastric <- mean_Gastric[rownames(dmp_IM_vs_Gas)]
dmp_IM_vs_Gas$Mean_IM      <- mean_IM[rownames(dmp_IM_vs_Gas)]
dmp_IM_vs_Gas$Delta_Beta   <- dmp_IM_vs_Gas$Mean_IM - dmp_IM_vs_Gas$Mean_Gastric

# --- Process Duodenal vs Gastric ---
dmp_Duo_vs_Gas$Mean_Gastric <- mean_Gastric[rownames(dmp_Duo_vs_Gas)]
dmp_Duo_vs_Gas$Mean_Duodenal<- mean_Duo[rownames(dmp_Duo_vs_Gas)]
dmp_Duo_vs_Gas$Delta_Beta   <- dmp_Duo_vs_Gas$Mean_Duodenal - dmp_Duo_vs_Gas$Mean_Gastric

# --- Process IM vs Duodenal ---
dmp_IM_vs_Duo$Mean_IM       <- mean_IM[rownames(dmp_IM_vs_Duo)]
dmp_IM_vs_Duo$Mean_Duodenal <- mean_Duo[rownames(dmp_IM_vs_Duo)]
dmp_IM_vs_Duo$Delta_Beta    <- dmp_IM_vs_Duo$Mean_IM - dmp_IM_vs_Duo$Mean_Duodenal

# Save updated files WITH Delta_Beta
cat("\nSaving final results WITH Delta_Beta...\n")
write.csv(dmp_IM_vs_Gas, "DMPs_IM_vs_Gastric.csv", row.names = FALSE)
write.csv(dmp_Duo_vs_Gas, "DMPs_Duodenal_vs_Gastric.csv", row.names = FALSE)
write.csv(dmp_IM_vs_Duo, "DMPs_IM_vs_Duodenal.csv", row.names = FALSE)

cat("\n✅ COMPLETE! All DMP tables saved with Delta_Beta.\n")


############################################
# Load RNA-seq data for integration
############################################

library(ggplot2)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(grid)
library(tidyr)
library(stringr)

# Load the RNA-seq data
DEG_IM_vs_Gastric <- read.csv("setdirectory/Simple_DEG_IM_vs_Gastric.csv", stringsAsFactors = FALSE)
DEG_Duodenal_vs_IM <- read.csv("setdirectory/Simple_DEG_Duodenal_vs_IM.csv", stringsAsFactors = FALSE)
DEG_Duodenal_vs_Gastric <- read.csv("setdirectory/Simple_DEG_Duodenal_vs_Gastric.csv", stringsAsFactors = FALSE)

# NOTE: We already have dmp_IM_vs_Gas, dmp_Duo_vs_Gas, dmp_IM_vs_Duo in memory with Delta_Beta calculated

# Convert to data.frames (if needed)
dmp <- as.data.frame(dmp_IM_vs_Gas)
dmp2 <- as.data.frame(dmp_Duo_vs_Gas)
dmp3 <- as.data.frame(dmp_IM_vs_Duo)

# Define thresholds for significant changes
meth_threshold <- 0.1  
fc_threshold <- 1      

# Function to create integrative plot with flexible labeling
create_integrative_plot <- function(deg_data, dmp_data, comparison_name, 
                                    meth_threshold = 0.1, fc_threshold = 1,
                                    label_strategy = "top_per_quadrant", top_n = 15,
                                    add_labels = TRUE) {
  
  cat(paste0("\n========== Processing ", comparison_name, " ==========\n"))
  
  ############################################################
  ## Sanity checks
  ############################################################
  stopifnot("UCSC_RefGene_Name"  %in% colnames(dmp_data))
  stopifnot("UCSC_RefGene_Group" %in% colnames(dmp_data))
  stopifnot("Delta_Beta" %in% colnames(dmp_data))
  
  ############################################################
  ## 1. Keep CpGs with gene annotation
  ############################################################
  dmp_clean <- dmp_data[
    !is.na(dmp_data$UCSC_RefGene_Name) &
      dmp_data$UCSC_RefGene_Name != "",
  ]
  
  ############################################################
  ## 2. Expand CpGs mapping to multiple genes
  ############################################################
  dmp_clean$UCSC_RefGene_Name <- gsub(",", ";", dmp_clean$UCSC_RefGene_Name)
  dmp_gene_long <- separate_rows(
    dmp_clean,
    UCSC_RefGene_Name,
    sep = ";"
  )
  colnames(dmp_gene_long)[
    colnames(dmp_gene_long) == "UCSC_RefGene_Name"
  ] <- "gene_symbol"
  
  ############################################################
  ## 3. Restrict to promoter CpGs
  ############################################################
  dmp_gene_long <- dmp_gene_long[
    grepl("TSS200|TSS1500", dmp_gene_long$UCSC_RefGene_Group),
  ]
  
  ############################################################
  ## 4. Aggregate CpGs → gene level using Delta_Beta
  ############################################################
  gene_meth <- dmp_gene_long %>%
    group_by(gene_symbol) %>%
    summarise(
      meth_delta = mean(Delta_Beta, na.rm = TRUE),
      n_cpg      = n(),
      .groups    = "drop"
    )
  
  ############################################################
  ## 5. Prepare DEG data
  ############################################################
  deg <- deg_data
  if("X" %in% colnames(deg)) {
    colnames(deg)[colnames(deg) == "X"] <- "gene_symbol"
  }
  
  ############################################################
  ## 6. Merge methylation and expression data
  ############################################################
  deg_meth <- merge(
    deg,
    gene_meth,
    by = "gene_symbol"
  )
  
  # Print data range
  cat("\nData ranges:\n")
  cat(paste0("Methylation (Delta_Beta): ", round(min(deg_meth$meth_delta, na.rm=TRUE), 3), 
             " to ", round(max(deg_meth$meth_delta, na.rm=TRUE), 3), "\n"))
  cat(paste0("Expression (log2FC): ", round(min(deg_meth$log2FoldChange, na.rm=TRUE), 3), 
             " to ", round(max(deg_meth$log2FoldChange, na.rm=TRUE), 3), "\n"))
  
  ############################################################
  ## 7. Classify genes - genes in extreme quadrants only
  ############################################################
  deg_meth <- deg_meth %>%
    mutate(
      quadrant = case_when(
        # Top right: Hypomethylated & Upregulated (blue)
        meth_delta < -meth_threshold & log2FoldChange > fc_threshold ~ "Hypomethylated-Upregulated",
        
        # Top left: Hypermethylated & Upregulated (green)
        meth_delta > meth_threshold & log2FoldChange > fc_threshold ~ "Hypermethylated-Upregulated",
        
        # Bottom left: Hypermethylated & Downregulated (red)
        meth_delta > meth_threshold & log2FoldChange < -fc_threshold ~ "Hypermethylated-Downregulated",
        
        # Bottom right: Hypomethylated & Downregulated (yellow)
        meth_delta < -meth_threshold & log2FoldChange < -fc_threshold ~ "Hypomethylated-Downregulated",
        
        # Everything else
        TRUE ~ "Not significant"
      ),
      # Calculate combined effect size (distance from origin for genes in extreme quadrants)
      effect_size = sqrt(meth_delta^2 + log2FoldChange^2)
    )
  
  ############################################################
  ## 8. Select genes to label based on strategy
  ############################################################
  if(label_strategy == "all_colored") {
    # Label ALL genes in colored quadrants
    genes_to_label <- deg_meth %>%
      filter(quadrant != "Not significant") %>%
      pull(gene_symbol)
    
    cat(paste0("\nLabeling strategy: ALL colored genes (", length(genes_to_label), " genes)\n"))
    
  } else if(label_strategy == "top_per_quadrant") {
    # Label top N genes from EACH colored quadrant
    genes_to_label <- deg_meth %>%
      filter(quadrant != "Not significant") %>%
      group_by(quadrant) %>%
      arrange(desc(effect_size)) %>%
      slice_head(n = top_n) %>%
      ungroup() %>%
      pull(gene_symbol)
    
    cat(paste0("\nLabeling strategy: Top ", top_n, " genes per quadrant (", 
               length(genes_to_label), " genes total)\n"))
    
    # Print breakdown by quadrant
    genes_by_quadrant <- deg_meth %>%
      filter(gene_symbol %in% genes_to_label) %>%
      group_by(quadrant) %>%
      summarise(n = n(), .groups = "drop")
    print(genes_by_quadrant)
  }
  
  deg_meth <- deg_meth %>%
    mutate(to_label = gene_symbol %in% genes_to_label)
  
  # Set factor levels
  deg_meth$quadrant <- factor(deg_meth$quadrant,
                              levels = c("Hypomethylated-Downregulated",
                                         "Hypomethylated-Upregulated",
                                         "Hypermethylated-Downregulated",
                                         "Hypermethylated-Upregulated",
                                         "Not significant"))
  
  # Define colors
  color_map <- c(
    "Hypomethylated-Downregulated" = "#FFD700",
    "Hypomethylated-Upregulated" = "#0000FF",
    "Hypermethylated-Downregulated" = "#FF0000",
    "Hypermethylated-Upregulated" = "#228B22",
    "Not significant" = "#808080"
  )
  
  # Print summary
  cat("\nQuadrant summary:\n")
  print(table(deg_meth$quadrant))
  cat(paste0("Total genes in extreme quadrants: ", sum(deg_meth$quadrant != "Not significant"), "\n"))
  cat(paste0("Genes to be labeled: ", sum(deg_meth$to_label), "\n\n"))
  
  ############################################################
  ## 9. Create plot
  ############################################################
  p <- ggplot(deg_meth, aes(x = meth_delta, y = log2FoldChange, color = quadrant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = color_map) +
    # Threshold lines showing the cutoffs
    geom_vline(xintercept = c(-meth_threshold, meth_threshold), 
               linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_hline(yintercept = c(-fc_threshold, fc_threshold), 
               linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    # Center lines
    geom_vline(xintercept = 0, linetype = "solid", colour = "black", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "solid", colour = "black", linewidth = 0.3) +
    theme_bw() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      plot.title = element_text(size = 11, face = "bold")
    ) +
    labs(
      x = "Mean promoter methylation difference (Delta-Beta)",
      y = "log2 expression fold change",
      title = comparison_name,
      color = "Classification"
    )
  
  # Add labels only if requested
  if(add_labels && sum(deg_meth$to_label) > 0) {
    p <- p + geom_text_repel(
      data = deg_meth %>% filter(to_label),
      aes(label = gene_symbol),
      size = 2.5,
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "grey50",
      segment.size = 0.3,
      min.segment.length = 0,
      color = "black",
      fontface = "bold",
      force = 2
    )
  }
  
  ############################################################
  ## 10. Save plots
  ############################################################
  ggsave(paste0("Integrative_plot_", gsub(" ", "_", comparison_name), ".png"),
         plot = p, width = 11, height = 8, units = "in", dpi = 300)
  
  ggsave(paste0("Integrative_plot_", gsub(" ", "_", comparison_name), ".pdf"),
         plot = p, width = 11, height = 8, units = "in")
  
  ############################################################
  ## 11. Export gene lists - only extreme quadrants
  ############################################################
  genes_extreme <- deg_meth %>% 
    filter(quadrant != "Not significant") %>%
    arrange(desc(effect_size)) %>%
    dplyr::select(gene_symbol, meth_delta, log2FoldChange, quadrant, effect_size, n_cpg, to_label, everything())
  
  write.csv(genes_extreme, 
            paste0("Genes_extreme_quadrants_", gsub(" ", "_", comparison_name), ".csv"), 
            row.names = FALSE)
  
  # Print top labeled genes
  cat("\nTop labeled genes by effect size:\n")
  labeled_genes <- deg_meth %>%
    filter(to_label) %>%
    arrange(desc(effect_size)) %>%
    dplyr::select(gene_symbol, meth_delta, log2FoldChange, quadrant, effect_size)
  print(head(labeled_genes, 20))
  cat("\n")
  
  return(list(plot = p, data = deg_meth))
}

############################################################
## Create plots for all three comparisons
############################################################
cat("\n############################################")
cat("\n### CREATING INTEGRATIVE ANALYSIS PLOTS ###")
cat("\n############################################\n")

# 1. IM vs Gastric - Label ALL colored genes
result1 <- create_integrative_plot(
  deg_data = DEG_IM_vs_Gastric,
  dmp_data = dmp,
  comparison_name = "IM vs Gastric",
  meth_threshold = meth_threshold,
  fc_threshold = fc_threshold,
  label_strategy = "all_colored",  # Label ALL colored points
  add_labels = TRUE
)

# 2. Duodenal vs Gastric - Label top 15 per quadrant
result2 <- create_integrative_plot(
  deg_data = DEG_Duodenal_vs_Gastric,
  dmp_data = dmp2,
  comparison_name = "Duodenal vs Gastric",
  meth_threshold = meth_threshold,
  fc_threshold = fc_threshold,
  label_strategy = "top_per_quadrant",
  top_n = 15,
  add_labels = TRUE
)

# 3. Duodenal vs IM - Label top 15 per quadrant
# Reverse Delta_Beta for correct orientation
dmp3_reversed <- dmp3
dmp3_reversed$Delta_Beta <- -dmp3_reversed$Delta_Beta

result3 <- create_integrative_plot(
  deg_data = DEG_Duodenal_vs_IM,
  dmp_data = dmp3_reversed,
  comparison_name = "Duodenal vs IM",
  meth_threshold = meth_threshold,
  fc_threshold = fc_threshold,
  label_strategy = "top_per_quadrant",  
  top_n = 15,
  add_labels = TRUE
)

cat("\n\n### All plots completed! ###\n")

############################################################
## Create combined figures
############################################################

# Combined plot without legend
combined_plot <- grid.arrange(
  result1$plot + theme(legend.position = "none"),
  result2$plot + theme(legend.position = "none"),
  result3$plot + theme(legend.position = "none"),
  ncol = 2,
  nrow = 2,
  top = textGrob("Integrative Analysis of Gene Expression and DNA Methylation", 
                 gp = gpar(fontsize = 14, fontface = "bold"))
)

ggsave("Integrative_plot_ALL_COMPARISONS.png",
       plot = combined_plot, width = 18, height = 14, units = "in", dpi = 300)

ggsave("Integrative_plot_ALL_COMPARISONS.pdf",
       plot = combined_plot, width = 18, height = 14, units = "in")

# Extract legend
get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(result1$plot)

# Combined plot WITH legend
combined_plot_with_legend <- grid.arrange(
  arrangeGrob(
    result1$plot + theme(legend.position = "none"),
    result2$plot + theme(legend.position = "none"),
    result3$plot + theme(legend.position = "none"),
    ncol = 2,
    nrow = 2,
    top = textGrob("Integrative Analysis of Gene Expression and DNA Methylation", 
                   gp = gpar(fontsize = 14, fontface = "bold"))
  ),
  legend,
  ncol = 2,
  widths = c(10, 1.5)
)

ggsave("Integrative_plot_ALL_COMPARISONS_with_legend.png",
       plot = combined_plot_with_legend, width = 20, height = 14, units = "in", dpi = 300)

ggsave("Integrative_plot_ALL_COMPARISONS_with_legend.pdf",
       plot = combined_plot_with_legend, width = 20, height = 14, units = "in")

# Script to create publication-ready tables for IM vs Gastric comparison
# Tables organized by methylation direction (Hypermethylated vs Hypomethylated)

library(dplyr)
library(tidyr)
library(openxlsx)

# Set working directory
workdir <- "setdirectory"
setwd(workdir)

cat("\n=== LOADING DATA ===\n")

############################################
# Load data files
############################################

# Load the integrative analysis results
genes_extreme <- read.csv("Genes_extreme_quadrants_IM_vs_Gastric.csv", stringsAsFactors = FALSE)

# Load the original DMP file to get methylation FDR
dmp_original <- read.csv("DMPs_IM_vs_Gastric.csv", stringsAsFactors = FALSE)

cat("Genes extreme loaded:", nrow(genes_extreme), "rows\n")
cat("DMP original loaded:", nrow(dmp_original), "rows\n")

############################################
# Create methylation FDR lookup
############################################

cat("\n=== CREATING METHYLATION FDR LOOKUP ===\n")

# For each gene, get the minimum (most significant) methylation FDR across all promoter CpGs
meth_fdr_lookup <- dmp_original %>%
  filter(!is.na(UCSC_RefGene_Name) & UCSC_RefGene_Name != "") %>%
  mutate(UCSC_RefGene_Name = gsub(",", ";", UCSC_RefGene_Name)) %>%
  separate_rows(UCSC_RefGene_Name, sep = ";") %>%
  filter(grepl("TSS200|TSS1500", UCSC_RefGene_Group)) %>%
  group_by(UCSC_RefGene_Name) %>%
  summarise(
    Methylation_FDR = min(adj.P.Val, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(gene_symbol = UCSC_RefGene_Name)

cat("Methylation FDR lookup created for", nrow(meth_fdr_lookup), "genes\n")

# Merge with genes_extreme
genes_extreme <- genes_extreme %>%
  left_join(meth_fdr_lookup, by = "gene_symbol")

cat("Methylation FDR merged\n")

############################################
# Table 1: Hypermethylated Genes
############################################

cat("\n=== CREATING TABLE 1: HYPERMETHYLATED GENES ===\n")

table1 <- genes_extreme %>%
  filter(grepl("Hypermethylated", quadrant)) %>%
  arrange(desc(effect_size)) %>%
  mutate(
    Expression_Direction = case_when(
      log2FoldChange > 0 ~ "Upregulated",
      log2FoldChange < 0 ~ "Downregulated",
      TRUE ~ "No change"
    )
  ) %>%
  dplyr::select(
    Gene = gene_symbol,
    `Promoter CpGs` = n_cpg,
    `Delta-Beta` = meth_delta,
    `log2FC` = log2FoldChange,
    `Expression Direction` = Expression_Direction,
    `Expression FDR` = padj,
    `Methylation FDR` = Methylation_FDR,
    `Effect Size` = effect_size
  ) %>%
  mutate(
    `Delta-Beta` = round(`Delta-Beta`, 4),
    log2FC = round(log2FC, 4),
    `Expression FDR` = round(`Expression FDR`, 4),
    `Methylation FDR` = round(`Methylation FDR`, 4),
    `Effect Size` = round(`Effect Size`, 4)
  )

cat("Table 1 created with", nrow(table1), "genes\n")
cat("  Downregulated:", sum(table1$`Expression Direction` == "Downregulated"), "\n")
cat("  Upregulated:", sum(table1$`Expression Direction` == "Upregulated"), "\n")

############################################
# Table 2: Hypomethylated Genes
############################################

cat("\n=== CREATING TABLE 2: HYPOMETHYLATED GENES ===\n")

table2 <- genes_extreme %>%
  filter(grepl("Hypomethylated", quadrant)) %>%
  arrange(desc(effect_size)) %>%
  mutate(
    Expression_Direction = case_when(
      log2FoldChange > 0 ~ "Upregulated",
      log2FoldChange < 0 ~ "Downregulated",
      TRUE ~ "No change"
    )
  ) %>%
  dplyr::select(
    Gene = gene_symbol,
    `Promoter CpGs` = n_cpg,
    `Delta-Beta` = meth_delta,
    `log2FC` = log2FoldChange,
    `Expression Direction` = Expression_Direction,
    `Expression FDR` = padj,
    `Methylation FDR` = Methylation_FDR,
    `Effect Size` = effect_size
  ) %>%
  mutate(
    `Delta-Beta` = round(`Delta-Beta`, 4),
    log2FC = round(log2FC, 4),
    `Expression FDR` = round(`Expression FDR`, 4),
    `Methylation FDR` = round(`Methylation FDR`, 4),
    `Effect Size` = round(`Effect Size`, 4)
  )

cat("Table 2 created with", nrow(table2), "genes\n")
cat("  Upregulated:", sum(table2$`Expression Direction` == "Upregulated"), "\n")
cat("  Downregulated:", sum(table2$`Expression Direction` == "Downregulated"), "\n")

############################################
# Save as CSV files
############################################

cat("\n=== SAVING CSV FILES ===\n")

write.csv(table1, "Table1_Hypermethylated_Genes_IM_vs_Gastric.csv", row.names = FALSE)
cat("Saved: Table1_Hypermethylated_Genes_IM_vs_Gastric.csv\n")

write.csv(table2, "Table2_Hypomethylated_Genes_IM_vs_Gastric.csv", row.names = FALSE)
cat("Saved: Table2_Hypomethylated_Genes_IM_vs_Gastric.csv\n")

############################################
# Save as Excel file with both sheets
############################################

cat("\n=== CREATING EXCEL FILE ===\n")

wb <- createWorkbook()

# Add Table 1
addWorksheet(wb, "Hypermethylated Genes")
writeData(wb, sheet = 1, table1)

# Add Table 2
addWorksheet(wb, "Hypomethylated Genes")
writeData(wb, sheet = 2, table2)

# Create header style
headerStyle <- createStyle(
  fontSize = 11,
  fontColour = "#FFFFFF",
  halign = "center",
  valign = "center",
  textDecoration = "bold",
  fgFill = "#4F81BD",
  border = "TopBottomLeftRight",
  borderColour = "#4F81BD"
)

# Apply header style
addStyle(wb, sheet = 1, headerStyle, rows = 1, cols = 1:ncol(table1), gridExpand = TRUE)
addStyle(wb, sheet = 2, headerStyle, rows = 1, cols = 1:ncol(table2), gridExpand = TRUE)

# Set column widths
setColWidths(wb, sheet = 1, cols = 1:ncol(table1), widths = "auto")
setColWidths(wb, sheet = 2, cols = 1:ncol(table2), widths = "auto")

# Save workbook
saveWorkbook(wb, "IM_vs_Gastric_Publication_Tables.xlsx", overwrite = TRUE)

cat("Saved: IM_vs_Gastric_Publication_Tables.xlsx\n")

############################################
# Print summary statistics
############################################

cat("\n=== SUMMARY STATISTICS ===\n\n")

cat("Table 1 - Hypermethylated Genes:\n")
cat("  Total genes:", nrow(table1), "\n")
cat("  Downregulated (concordant):", sum(table1$`Expression Direction` == "Downregulated"), 
    paste0("(", round(100*sum(table1$`Expression Direction` == "Downregulated")/nrow(table1), 1), "%)"), "\n")
cat("  Upregulated (discordant):", sum(table1$`Expression Direction` == "Upregulated"),
    paste0("(", round(100*sum(table1$`Expression Direction` == "Upregulated")/nrow(table1), 1), "%)"), "\n")
cat("  Mean Delta-Beta:", round(mean(table1$`Delta-Beta`), 3), "\n")
cat("  Mean log2FC:", round(mean(table1$log2FC), 3), "\n\n")

cat("Table 2 - Hypomethylated Genes:\n")
cat("  Total genes:", nrow(table2), "\n")
cat("  Upregulated (concordant):", sum(table2$`Expression Direction` == "Upregulated"),
    paste0("(", round(100*sum(table2$`Expression Direction` == "Upregulated")/nrow(table2), 1), "%)"), "\n")
cat("  Downregulated (discordant):", sum(table2$`Expression Direction` == "Downregulated"),
    paste0("(", round(100*sum(table2$`Expression Direction` == "Downregulated")/nrow(table2), 1), "%)"), "\n")
cat("  Mean Delta-Beta:", round(mean(table2$`Delta-Beta`), 3), "\n")
cat("  Mean log2FC:", round(mean(table2$log2FC), 3), "\n\n")

############################################
# Print top 10 genes from each table
############################################

cat("=== TOP 10 HYPERMETHYLATED GENES (by effect size) ===\n")
print(head(table1 %>% dplyr::select(Gene, `Delta-Beta`, log2FC, `Expression Direction`, `Effect Size`), 10))

cat("\n=== TOP 10 HYPOMETHYLATED GENES (by effect size) ===\n")
print(head(table2 %>% dplyr::select(Gene, `Delta-Beta`, log2FC, `Expression Direction`, `Effect Size`), 10))

cat("\n✅ ALL TABLES GENERATED SUCCESSFULLY! ✅\n")
