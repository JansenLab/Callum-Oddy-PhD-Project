if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Unknownmanifest")

# ============================================================================
# Comprehensive DNA Methylation Analysis Pipeline for EPICv2
# For Gastric, IM, and Duodenal Tissue Comparison
# ============================================================================

# Load required packages
library(ENmix)
library(sesame)
library(minfi)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(IlluminaHumanMethylationEPICv2manifest)
library(DMRcate)
library(limma)
library(sva)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(missMethyl)
library(Unknownmanifest)

# Set working directory (for outputs)
setwd("setworkingdirectory")

# Define path to raw data
data_path <- "setworkingdirectory/IDATs/All Data"

# ============================================================================
# 1. DATA IMPORT AND INITIAL QC
# ============================================================================

# Read in IDAT files using ENmix (superior to ChAMP for EPICv2), directly read the manifest file
rgSet <- readidat(path = data_path, 
                  recursive = TRUE,
                  manifestfile = NULL)  # Will auto-detect EPICv2

# Read sample sheet
targets <- read.csv(file.path(data_path, "Sample_Data_noblood.csv"), 
                    stringsAsFactors = FALSE)

# Create proper rownames from Basename column
rownames(targets) <- targets$Basename

# Ensure targets match rgSet samples
# rgSet sample names should match the Basename column
sample_names <- colnames(rgSet)
if(!all(sample_names %in% rownames(targets))) {
  cat("Warning: Not all samples in rgSet found in targets file\n")
  cat("Attempting to match samples...\n")
  # Try to match by extracting sentrix info from colnames
  targets <- targets[targets$Basename %in% sample_names, ]
  targets <- targets[match(sample_names, targets$Basename), ]
  rownames(targets) <- targets$Basename
}

# Add full path to Basename for minfi compatibility
targets$Basename <- file.path(data_path, targets$Basename)

colData(rgSet) <- DataFrame(targets)

# Initial QC plots
pdf("01_QC_plots.pdf", width = 12, height = 8)

# Control probe plots
plotCtrl(rgSet, IDorder = colnames(rgSet))

# Quality metrics using ENmix
qc <- QCinfo(rgSet)
print(summary(qc))

# Plot sample quality
plotQC(qc, IDorder = colnames(rgSet))

# Detection p-values
detP <- pval(rgSet)
barplot(colMeans(detP > 0.01), las = 2, 
        main = "Mean detection P-values > 0.01 per sample",
        ylab = "Proportion", col = "steelblue")
abline(h = 0.01, col = "red", lty = 2)

dev.off()

# ============================================================================
# 2. PREPROCESSING WITH ENmix
# ============================================================================

# ENmix preprocessing pipeline (superior performance vs ChAMP)
# This includes:
# - Background correction (exponential-normal mixture model)
# - Dye bias correction
# - Probe type bias correction (RCP method)

mset <- preprocessENmix(rgSet, 
                        bgParaEst = "oob",  # Use out-of-band for background
                        dyeCorr = "RELIC",   # Dye bias correction
                        QCinfo = qc,
                        nCores = 4)

# Quality filtering
# Remove poor quality probes and samples
mset_filtered <- qcfilter(mset,
                          qcscore = NULL,
                          badSampleCutoff = 0.05,  # Remove samples with >5% poor probes
                          badCpgCutoff = 0.05,     # Remove probes failing in >5% samples
                          detPcut = 0.01,
                          detSamplecut = 0.1)

# ============================================================================
# 3. PROBE FILTERING
# ============================================================================

# Get beta values for filtering
beta <- getB(mset_filtered)

# Filter probes based on:
# 1. Cross-reactive probes
# 2. SNP-associated probes
# 3. Sex chromosomes (if not relevant to your study)
# 4. Non-CpG probes
# 5. Multi-mapping probes

# Load probe annotations
ann <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

# Remove sex chromosome probes (adjust if sex is relevant)
keep_probes <- !(ann$chr %in% c("chrX", "chrY"))

# Remove probes with SNPs at CpG site
keep_probes <- keep_probes & is.na(ann$Probe_rs)

# Remove cross-reactive probes (from Peters et al. 2024)
# You may need to load this separately
# cross_reactive <- read.csv("EPICv2_cross_reactive_probes.csv")
# keep_probes <- keep_probes & !(rownames(beta) %in% cross_reactive$probe)

# Remove ch probes (non-CpG)
keep_probes <- keep_probes & !grepl("^ch", rownames(beta))

# Apply filtering
beta_filtered <- beta[keep_probes, ]
mset_filtered <- mset_filtered[keep_probes, ]

cat(sprintf("Retained %d probes after filtering\n", nrow(beta_filtered)))

# ============================================================================
# 4. NORMALIZATION AND BATCH CORRECTION
# ============================================================================

# Get M-values for statistical analysis (better for linear modeling)
mval <- getM(mset_filtered)

# Check for batch effects
# Visualize with PCA before correction
pca <- prcomp(t(mval), scale = TRUE)

pdf("02_PCA_before_correction.pdf", width = 10, height = 8)
pca_data <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
                       Tissue = targets$Sample_Type,
                       Slide = as.factor(targets$Sentrix_ID),
                       Patient = as.factor(targets$Patient))

ggplot(pca_data, aes(x = PC1, y = PC2, color = Tissue, shape = Slide)) +
  geom_point(size = 4) +
  theme_bw() +
  labs(title = "PCA before batch correction",
       x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)"))

ggplot(pca_data, aes(x = PC1, y = PC2, color = Patient, shape = Tissue)) +
  geom_point(size = 4) +
  theme_bw() +
  labs(title = "PCA - by Patient")
dev.off()

# Batch correction using ComBat (if needed based on PCA)
# Only correct for technical batches, not biological variables
if(length(unique(targets$Sentrix_ID)) > 1) {
  mod <- model.matrix(~ Sample_Type, data = targets)
  mval_combat <- ComBat(dat = mval, 
                        batch = targets$Sentrix_ID, 
                        mod = mod)
  beta_combat <- ilogit2(mval_combat)
} else {
  mval_combat <- mval
  beta_combat <- beta_filtered
}

# ============================================================================
# 5. EXPLORATORY ANALYSIS
# ============================================================================

pdf("03_exploratory_analysis.pdf", width = 12, height = 10)

# Density plots
par(mfrow = c(2,2))
densityPlot(beta_filtered, sampGroups = targets$Sample_Type,
            main = "Beta values by tissue type", legend = TRUE)
densityPlot(mval, sampGroups = targets$Sample_Type,
            main = "M-values by tissue type", legend = TRUE)

# Hierarchical clustering
dist_mat <- dist(t(mval_combat))
hc <- hclust(dist_mat)
colors <- brewer.pal(3, "Set1")[as.factor(targets$Sample_Type)]
plot(hc, labels = targets$Sample_Name, main = "Hierarchical Clustering")

# Heatmap of most variable probes
var_probes <- head(order(apply(mval_combat, 1, var), decreasing = TRUE), 1000)
pheatmap(mval_combat[var_probes, ],
         annotation_col = data.frame(
           Tissue = targets$Sample_Type,
           Patient = targets$Patient,
           row.names = colnames(mval_combat)
         ),
         show_rownames = FALSE,
         scale = "row",
         main = "Top 1000 most variable CpGs")

# PCA after correction
pca_corrected <- prcomp(t(mval_combat), scale = TRUE)
pca_corr_data <- data.frame(
  PC1 = pca_corrected$x[,1], 
  PC2 = pca_corrected$x[,2],
  PC3 = pca_corrected$x[,3],
  Tissue = targets$Sample_Type,
  Patient = as.factor(targets$Patient)
)

print(ggplot(pca_corr_data, aes(x = PC1, y = PC2, color = Tissue, shape = Patient)) +
        geom_point(size = 4) +
        theme_bw() +
        labs(title = "PCA after batch correction"))

print(ggplot(pca_corr_data, aes(x = PC2, y = PC3, color = Tissue, shape = Patient)) +
        geom_point(size = 4) +
        theme_bw())

dev.off()

# ============================================================================
# 6. DIFFERENTIAL METHYLATION ANALYSIS - PAIRED DESIGN
# ============================================================================

# Since you have patient-matched samples, use paired design
# Compare: Gastric vs IM, IM vs Duodenal, Gastric vs Duodenal

# Create design matrix for paired analysis
# Based on your data: Patient column, Sample_Type for tissue
Patient <- factor(targets$Patient)
Tissue <- factor(targets$Sample_Type, 
                 levels = c("Gastric", "IM", "Duodenal"))

design <- model.matrix(~ Patient + Tissue)

# Fit linear model
fit <- lmFit(mval_combat, design)
fit <- eBayes(fit, robust = TRUE)

# Extract results for each comparison
# Adjust contrast based on your specific tissue types

# IM vs Gastric
contrast_IM_Gastric <- makeContrasts(TissueIM - TissueGastric, levels = design)
fit_IM_Gastric <- contrasts.fit(fit, contrast_IM_Gastric)
fit_IM_Gastric <- eBayes(fit_IM_Gastric, robust = TRUE)

# Duodenal vs IM
contrast_Duo_IM <- makeContrasts(TissueDuodenal - TissueIM, levels = design)
fit_Duo_IM <- contrasts.fit(fit, contrast_Duo_IM)
fit_Duo_IM <- eBayes(fit_Duo_IM, robust = TRUE)

# Duodenal vs Gastric
contrast_Duo_Gastric <- makeContrasts(TissueDuodenal - TissueGastric, levels = design)
fit_Duo_Gastric <- contrasts.fit(fit, contrast_Duo_Gastric)
fit_Duo_Gastric <- eBayes(fit_Duo_Gastric, robust = TRUE)

# Get DMPs with appropriate significance threshold
# For EPIC arrays: p < 9e-8 (Bonferroni corrected)
sig_threshold <- 9e-8
deltabeta_threshold <- 0.1  # 10% methylation difference

# IM vs Gastric DMPs
dmp_IM_Gastric <- topTable(fit_IM_Gastric, num = Inf, coef = 1)
dmp_IM_Gastric$deltaBeta <- rowMeans(beta_combat[rownames(dmp_IM_Gastric), 
                                                 targets$Sample_Type == "IM"]) -
  rowMeans(beta_combat[rownames(dmp_IM_Gastric), 
                       targets$Sample_Type == "Gastric"])
dmp_IM_Gastric_sig <- dmp_IM_Gastric[dmp_IM_Gastric$adj.P.Val < 0.05 & 
                                       abs(dmp_IM_Gastric$deltaBeta) > deltabeta_threshold, ]

# Duodenal vs IM DMPs
dmp_Duo_IM <- topTable(fit_Duo_IM, num = Inf, coef = 1)
dmp_Duo_IM$deltaBeta <- rowMeans(beta_combat[rownames(dmp_Duo_IM), 
                                             targets$Sample_Type == "Duodenal"]) -
  rowMeans(beta_combat[rownames(dmp_Duo_IM), 
                       targets$Sample_Type == "IM"])
dmp_Duo_IM_sig <- dmp_Duo_IM[dmp_Duo_IM$adj.P.Val < 0.05 & 
                               abs(dmp_Duo_IM$deltaBeta) > deltabeta_threshold, ]

# Duodenal vs Gastric DMPs
dmp_Duo_Gastric <- topTable(fit_Duo_Gastric, num = Inf, coef = 1)
dmp_Duo_Gastric$deltaBeta <- rowMeans(beta_combat[rownames(dmp_Duo_Gastric), 
                                                  targets$Sample_Type == "Duodenal"]) -
  rowMeans(beta_combat[rownames(dmp_Duo_Gastric), 
                       targets$Sample_Type == "Gastric"])
dmp_Duo_Gastric_sig <- dmp_Duo_Gastric[dmp_Duo_Gastric$adj.P.Val < 0.05 & 
                                         abs(dmp_Duo_Gastric$deltaBeta) > deltabeta_threshold, ]

# Save results
write.csv(dmp_IM_Gastric_sig, "DMP_IM_vs_Gastric_significant.csv")
write.csv(dmp_Duo_IM_sig, "DMP_Duodenal_vs_IM_significant.csv")
write.csv(dmp_Duo_Gastric_sig, "DMP_Duodenal_vs_Gastric_significant.csv")

cat(sprintf("Found %d significant DMPs: IM vs Gastric\n", nrow(dmp_IM_Gastric_sig)))
cat(sprintf("Found %d significant DMPs: Duodenal vs IM\n", nrow(dmp_Duo_IM_sig)))
cat(sprintf("Found %d significant DMPs: Duodenal vs Gastric\n", nrow(dmp_Duo_Gastric_sig)))

# ============================================================================
# 7. VISUALIZATIONS OF DMPs
# ============================================================================

pdf("04_DMP_visualizations.pdf", width = 12, height = 10)

# Volcano plots
make_volcano <- function(results, title, fc_col = "deltaBeta") {
  results$sig <- ifelse(results$adj.P.Val < 0.05 & abs(results[[fc_col]]) > 0.1,
                        "Significant", "Not Significant")
  
  ggplot(results, aes(x = get(fc_col), y = -log10(P.Value), color = sig)) +
    geom_point(alpha = 0.4, size = 1) +
    scale_color_manual(values = c("grey", "red")) +
    geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_bw() +
    labs(title = title, x = "Delta Beta", y = "-log10(P-value)") +
    theme(legend.position = "bottom")
}

print(make_volcano(dmp_IM_Gastric, "IM vs Gastric"))
print(make_volcano(dmp_Duo_IM, "Duodenal vs IM"))
print(make_volcano(dmp_Duo_Gastric, "Duodenal vs Gastric"))

# Heatmap of significant DMPs
if(nrow(dmp_IM_Gastric_sig) > 2) {
  top_dmps <- rownames(head(dmp_IM_Gastric_sig[order(dmp_IM_Gastric_sig$adj.P.Val), ], 50))
  
  pheatmap(beta_combat[top_dmps, ],
           annotation_col = data.frame(
             Tissue = targets$Sample_Type,
             Patient = targets$Patient,
             row.names = colnames(beta_combat)
           ),
           scale = "row",
           main = "Top 50 DMPs: IM vs Gastric",
           show_rownames = FALSE)
}

dev.off()

# ============================================================================
# 8. DIFFERENTIALLY METHYLATED REGIONS (DMR) ANALYSIS
# ============================================================================

# DMRcate for region-level analysis
# Note: For EPICv2, ensure you have the latest DMRcate version

# Annotate CpGs with genomic context
myAnnotation <- cpg.annotate(object = mval_combat, 
                             datatype = "array", 
                             what = "M",
                             arraytype = "EPICv2",
                             analysis.type = "differential",
                             design = design,
                             coef = "TissueIM")  # Adjust coefficient name

# Find DMRs
dmrs_IM_Gastric <- dmrcate(myAnnotation, lambda = 1000, C = 2)
dmr_results_IM_Gastric <- extractRanges(dmrs_IM_Gastric)

# Save DMR results
write.csv(as.data.frame(dmr_results_IM_Gastric), "DMR_IM_vs_Gastric.csv")

# Visualize top DMRs
pdf("05_DMR_plots.pdf", width = 12, height = 8)
if(length(dmr_results_IM_Gastric) > 0) {
  # Plot top DMRs
  for(i in 1:min(5, length(dmr_results_IM_Gastric))) {
    DMR.plot(ranges = dmr_results_IM_Gastric, 
             dmr = i,
             CpGs = beta_combat,
             what = "Beta",
             arraytype = "EPICv2",
             phen.col = as.numeric(as.factor(targets$Tissue_Type)),
             genome = "hg38")
  }
}
dev.off()

# ============================================================================
# 9. GENE SET ENRICHMENT ANALYSIS
# ============================================================================

# Using missMethyl for proper gene set testing, accounts for probe number bias

# Get significant CpGs
sig_cpgs_IM_Gastric <- rownames(dmp_IM_Gastric_sig)
all_cpgs <- rownames(dmp_IM_Gastric)

# GO enrichment
gst_go <- gometh(sig.cpg = sig_cpgs_IM_Gastric, 
                 all.cpg = all_cpgs,
                 collection = "GO",
                 array.type = "EPICv2",
                 plot.bias = TRUE)

# KEGG pathway enrichment
gst_kegg <- gometh(sig.cpg = sig_cpgs_IM_Gastric,
                   all.cpg = all_cpgs,
                   collection = "KEGG",
                   array.type = "EPICv2")

# Save enrichment results
write.csv(gst_go[gst_go$FDR < 0.05, ], "GO_enrichment_IM_vs_Gastric.csv")
write.csv(gst_kegg[gst_kegg$FDR < 0.05, ], "KEGG_enrichment_IM_vs_Gastric.csv")

# Plot top enriched terms
pdf("06_pathway_enrichment.pdf", width = 10, height = 8)

top_go <- head(gst_go[order(gst_go$FDR), ], 20)
top_go$Term <- factor(top_go$TERM, levels = rev(top_go$TERM))

ggplot(top_go, aes(x = -log10(FDR), y = Term)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_bw() +
  labs(title = "Top GO Terms Enriched in IM vs Gastric",
       x = "-log10(FDR)", y = "")

dev.off()

# ============================================================================
# 10. SUMMARY REPORT
# ============================================================================

sink("analysis_summary.txt")
cat("=== DNA Methylation Analysis Summary ===\n\n")
cat(sprintf("Total samples analyzed: %d\n", ncol(beta_combat)))
cat(sprintf("Total probes after filtering: %d\n", nrow(beta_combat)))
cat("\n--- Tissue Distribution ---\n")
print(table(targets$Sample_Type))
cat("\n--- Patient Distribution ---\n")
print(table(targets$Patient))
cat("\n--- Significant DMPs (FDR < 0.05, |deltaBeta| > 0.1) ---\n")
cat(sprintf("IM vs Gastric: %d\n", nrow(dmp_IM_Gastric_sig)))
cat(sprintf("Duodenal vs IM: %d\n", nrow(dmp_Duo_IM_sig)))
cat(sprintf("Duodenal vs Gastric: %d\n", nrow(dmp_Duo_Gastric_sig)))
cat("\n--- DMR Results ---\n")
cat(sprintf("DMRs identified (IM vs Gastric): %d\n", length(dmr_results_IM_Gastric)))
cat("\n--- GO Enrichment ---\n")
cat(sprintf("Significant GO terms (FDR < 0.05): %d\n", sum(gst_go$FDR < 0.05)))
cat("\n--- KEGG Enrichment ---\n")
cat(sprintf("Significant pathways (FDR < 0.05): %d\n", sum(gst_kegg$FDR < 0.05)))
sink()

cat("\n=== Analysis Complete! ===\n")
cat("Check the output PDF files and CSV tables for detailed results.\n")
