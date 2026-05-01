###############################################################################
# RNA-seq DESeq2 Analysis + Correct Gene Plotting (Wald + VST)
###############################################################################
setwd('setworkingdirection')


suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(ggplot2)
})

# ============================================================================
# COLOR SCHEMES
# ============================================================================
color_palette <- c("Gastric" = "#C75A5A",
                   "IM"      = "#9370B8",
                   "Duodenal" = "#6B8DB8")

patient_palette <- c("Patient3" = "#E69F00",
                     "Patient4" = "#56B4E9",
                     "Patient6" = "#009E73")

# ============================================================================
# DATA LOADING
# ============================================================================
cat("=== Loading Data ===\n\n")

merged_counts_file <- "setworkingdirectory/salmon.merged.gene_counts.tsv"
sample_sheet <- "setworkingdirectory/sample_metadata.csv"

cts <- read.csv(merged_counts_file, sep="\t", header=TRUE, row.names=1)
colData <- read.csv(sample_sheet, row.names=1)

# Keep sample columns only
ctsNew <- cts[, 2:11]
ctsNew <- round(ctsNew)

# Fix sample names if needed
if("Im13" %in% rownames(colData)) {
  rownames(colData)[rownames(colData) == "Im13"] <- "IM13"
}
if("Im13" %in% colnames(ctsNew)) {
  colnames(ctsNew)[colnames(ctsNew) == "Im13"] <- "IM13"
}

# Phenotype factor
colData$phenotype <- factor(colData$phenotype,
                            levels=c("Gastric","IM","Duodenal"))

# Patient grouping
colData$Patient <- factor(case_when(
  grepl("3", rownames(colData)) ~ "Patient3",
  grepl("4", rownames(colData)) ~ "Patient4",
  grepl("6", rownames(colData)) ~ "Patient6"
))

cat("Samples loaded:", ncol(ctsNew), "\n")
cat("Genes before filtering:", nrow(ctsNew), "\n\n")

# ============================================================================
# QUALITY CONTROL
# ============================================================================
cat("=== Quality Control Analysis ===\n\n")

qc_metrics <- data.frame(
  Sample_ID       = colnames(ctsNew),
  Phenotype       = colData$phenotype,
  Patient         = colData$Patient,
  Total_Reads     = colSums(ctsNew),
  Genes_Detected  = colSums(ctsNew > 0),
  Percent_Detected = colSums(ctsNew > 0) / nrow(ctsNew) * 100
)

qc_summary <- qc_metrics %>%
  group_by(Phenotype) %>%
  summarise(
    N = n(),
    Mean_Reads = mean(Total_Reads),
    SD_Reads   = sd(Total_Reads),
    Mean_Genes = mean(Genes_Detected),
    .groups = "drop"
  )

print(qc_summary)

# ============================================================================
# DESeq2 MODEL (WITH PATIENT BLOCKING)
# ============================================================================
cat("\n=== Running DESeq2 (Design: ~ Patient + phenotype) ===\n\n")

dds <- DESeqDataSetFromMatrix(countData = ctsNew,
                              colData = colData,
                              design = ~ Patient + phenotype)

dds <- DESeq(dds)

# ============================================================================
# CONTRASTS
# ============================================================================
cat("=== Running contrasts ===\n\n")

res_GvIM  <- results(dds, contrast=c("phenotype","IM","Gastric"))
res_GvDuo <- results(dds, contrast=c("phenotype","Duodenal","Gastric"))
res_IMvD  <- results(dds, contrast=c("phenotype","Duodenal","IM"))

# ============================================================================
# VST TRANSFORMATION
# ============================================================================
cat("=== Performing VST ===\n\n")

vsd <- vst(dds, blind=FALSE)
vst_mat <- assay(vsd)

###############################################################################
# IMPROVED GENE PLOTTING WITH BETTER TEXT POSITIONING
###############################################################################

library(ggplot2)
library(gridExtra)
library(grid)


# ============================================================================
# IMPROVED PLOT FUNCTION - Option 2: Text in Legend Area
# ============================================================================
plot_gene_text_side <- function(gene) {
  
  if(!gene %in% rownames(vst_mat)) {
    stop(paste("Gene", gene, "not found in count matrix."))
  }
  
  df <- data.frame(
    Expression = vst_mat[gene, ],
    Sample     = colnames(vst_mat),
    Phenotype  = colData$phenotype,
    Patient    = colData$Patient
  )
  
  # Get DESeq2 significance
  padj_GvIM  <- res_GvIM[gene, "padj"]
  padj_GvDuo <- res_GvDuo[gene, "padj"]
  padj_IMvD  <- res_IMvD[gene, "padj"]
  
  # Convert to stars
  p_to_star <- function(p) {
    if (is.na(p)) return("ns")
    if (p < 0.001) return("***")
    if (p < 0.01)  return("**")
    if (p < 0.05)  return("*")
    return("ns")
  }
  
  star_GvIM  <- p_to_star(padj_GvIM)
  star_GvDuo <- p_to_star(padj_GvDuo)
  star_IMvD  <- p_to_star(padj_IMvD)
  
  ymax <- max(df$Expression)
  ymin <- min(df$Expression)
  y_range <- ymax - ymin
  
  # Create p-value text for caption
  caption_text <- paste0(
    "IM vs Gastric: padj = ", format.pval(padj_GvIM, digits=3), "\n",
    "Duodenal vs Gastric: padj = ", format.pval(padj_GvDuo, digits=3), "\n",
    "Duodenal vs IM: padj = ", format.pval(padj_IMvD, digits=3)
  )
  
  p <- ggplot(df, aes(x=Phenotype, y=Expression, fill=Phenotype)) +
    geom_boxplot(outlier.shape=NA, alpha=0.7) +
    geom_jitter(aes(color=Patient), width=0.18, size=2.2) +
    scale_fill_manual(values=color_palette) +
    scale_color_manual(values=patient_palette) +
    theme_classic(base_size = 14) +
    labs(title=paste0("Gene: ", gene), 
         y="VST Expression",
         caption=caption_text) +
    
    # Significance brackets
    geom_segment(aes(x=1, xend=2, y=ymax+0.15*y_range, yend=ymax+0.15*y_range), 
                 linewidth=0.5) +
    annotate("text", x=1.5, y=ymax+0.18*y_range, label=star_GvIM, size=5) +
    
    geom_segment(aes(x=1, xend=3, y=ymax+0.28*y_range, yend=ymax+0.28*y_range), 
                 linewidth=0.5) +
    annotate("text", x=2, y=ymax+0.31*y_range, label=star_GvDuo, size=5) +
    
    geom_segment(aes(x=2, xend=3, y=ymax+0.41*y_range, yend=ymax+0.41*y_range), 
                 linewidth=0.5) +
    annotate("text", x=2.5, y=ymax+0.44*y_range, label=star_IMvD, size=5) +
    
    coord_cartesian(ylim=c(ymin - 0.05*y_range, ymax + 0.5*y_range), clip="off") +
    theme(
      plot.margin = margin(t=10, r=10, b=10, l=10),
      plot.caption = element_text(hjust=0, size=9, lineheight=1.2)
    )
  
  return(p)
}

# ============================================================================
# BATCH PLOTTING FUNCTIONS
# ============================================================================

# Plot and save with text as caption
batch_plot_text_side <- function(genes, folder) {
  dir.create(folder, showWarnings=FALSE)
  for (g in genes) {
    pdf(file = paste0(folder, "/", g, ".pdf"), width=6, height=5.5)
    print(plot_gene_text_side(g))
    dev.off()
  }
}


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

cat("\n=== Three improved plot functions available: ===\n")
cat("2. plot_gene_text_side(gene)    - P-values in caption area\n")

cat("Batch plotting functions:\n")
cat("- batch_plot_text_side(genes, folder)\n")

# ============================================================================
# EXAMPLE: Batch plot top 20 genes with your preferred method
# ============================================================================
# Helper: get top N DE genes from a DESeq2 results object
get_top_genes <- function(res, n=20) {
  res_df <- as.data.frame(res)
  res_df <- res_df[!is.na(res_df$padj), ]
  res_df <- res_df[order(res_df$padj), ]
  head(rownames(res_df), n)
}

# Get top genes for each comparison
top40_GvIM  <- get_top_genes(res_GvIM, 40)
top20_GvDuo <- get_top_genes(res_GvDuo, 20)
top20_IMvD  <- get_top_genes(res_IMvD, 20)

cat("Top 20 genes (IM vs Gastric):\n"); print(top40_GvIM)
cat("\nTop 20 genes (Duodenal vs Gastric):\n"); print(top20_GvDuo)
cat("\nTop 20 genes (Duodenal vs IM):\n"); print(top20_IMvD)

# Method 2: Text as caption (RECOMMENDED - most compact)
batch_plot_text_side(top40_GvIM, "plots_improved_IM_vs_Gastric")
batch_plot_text_side(top20_GvDuo, "plots_improved_Gastric_vs_Duo")
batch_plot_text_side(top20_IMvD, "plots_improved_IM_vs_Duo")



# ============================================================================
# Plot some markers
# ============================================================================
###############################################################################
# LINEAGE MARKER GENE SETS (Literature-curated)
###############################################################################

# ---------------------------
# INTESTINAL MARKERS
# ---------------------------
intestinal_markers <- c(
  # Master intestinal transcription factors
  "CDX1", "CDX2",
  
  # Enterocyte / absorptive lineage
  "VIL1", "ALPI", "FABP1", "FABP2",
  "APOA1", "APOA4",
  "SI", "SLC26A3", "SLC5A1",
  
  # Goblet cell / secretory
  "MUC2", "SPDEF", "AGR2",
  
  # Paneth cell
  "LYZ", "DEFA5", "DEFA6",
  
  # Intestinal stem / Wnt
  "LGR5", "ASCL2", "OLFM4",
  
  # Intestinal epithelial adhesion
  "CDH17"
)

# ---------------------------
# GASTRIC MARKERS
# ---------------------------
gastric_markers <- c(
  # Gastric lineage transcription factors
  "SOX2", "GATA4", "GATA6",
  
  # Gastric pit / foveolar
  "MUC5AC", "TFF1", "GKN1", "GKN2",
  
  # Gastric gland / neck
  "MUC6", "TFF2", "BHLHA15",
  
  # Chief cell
  "PGC", "GIF", "CPA1",
  
  # Parietal cell
  "ATP4A", "ATP4B", "SLC26A7",
  
  # Gastric epithelial identity
  "CLDN18", "KRT18"
)

# ---------------------------
# INTESTINAL METAPLASIA (IM) MARKERS
# ---------------------------
IM_markers <- c(
  # Core IM transcriptional drivers
  "CDX2", "SOX9", "HNF4A",
  
  # Hybrid gastric–intestinal mucins (REQUESTED)
  "MUC2", "MUC1", "MUC5AC", "MUC6",
  
  # Intestinal stem / progenitor features (REQUESTED)
  "LGR5", "OLFM4",
  
  # Intestinal epithelial identity (REQUESTED)
  "CDH17", "KRT20",
  
  # IM-associated secretory / trefoil factors (REQUESTED)
  "TFF3",
  
  # Hep antigen (IM / Barrett’s / intestinal differentiation)
  "ANXA10",
  
  # Metaplasia-associated differentiation
  "AGR2", "SPDEF", "REG4",
  
  # Wnt / regeneration
  "AXIN2", "RNF43",
  
  # Tight junction / IM-associated claudins
  "CLDN7", "CLDN10"
)

batch_plot_text_side(intestinal_markers, "plots_improved_intestinal_markers")
batch_plot_text_side(gastric_markers, "plots_improved_gastric_markers")
batch_plot_text_side(IM_markers, "plots_improved_IM_markers")


###############################################################################
# 1. SUBSET TO MARKER GENES
###############################################################################

marker_genes <- unique(c(gastric_markers, IM_markers, intestinal_markers))

heatmap_mat <- vst_mat[rownames(vst_mat) %in% marker_genes, ]

###############################################################################
# 2. Z-SCORE PER GENE
###############################################################################

heatmap_mat_z <- t(scale(t(heatmap_mat)))

###############################################################################
# 3. REMOVE NA / ZERO-VARIANCE GENES
###############################################################################

keep_rows <- apply(heatmap_mat_z, 1, function(x) all(is.finite(x)))

heatmap_mat_z <- heatmap_mat_z[keep_rows, , drop = FALSE]

###############################################################################
# 4. ASSIGN MARKER COLLECTION
###############################################################################

marker_set <- ifelse(
  rownames(heatmap_mat_z) %in% gastric_markers, "Gastric",
  ifelse(
    rownames(heatmap_mat_z) %in% IM_markers, "IM",
    ifelse(
      rownames(heatmap_mat_z) %in% intestinal_markers, "Intestinal",
      NA
    )
  )
)

marker_set <- factor(
  marker_set,
  levels = c("Gastric", "IM", "Intestinal")
)

###############################################################################
# 5. ORDER GENES BY COLLECTION
###############################################################################

order_rows <- order(marker_set)
heatmap_mat_z <- heatmap_mat_z[order_rows, , drop = FALSE]
marker_set    <- marker_set[order_rows]

###############################################################################
# 6. DEFINE GAPS BETWEEN COLLECTIONS
###############################################################################

gap_positions <- cumsum(table(marker_set))

###############################################################################
# 7. COLUMN ANNOTATION (METADATA)
###############################################################################

annotation_col <- data.frame(
  Phenotype = colData$phenotype,
  Patient   = colData$Patient
)

rownames(annotation_col) <- rownames(colData)

# Ensure column order matches matrix
annotation_col <- annotation_col[colnames(heatmap_mat_z), , drop = FALSE]

###############################################################################
# 8. ANNOTATION COLOURS
###############################################################################

annotation_colors <- list(
  Phenotype = c(
    "Gastric"   = "#C75A5A",
    "IM"        = "#9370B8",
    "Duodenal"  = "#6B8DB8"
  ),
  Patient = c(
    "Patient3" = "#E69F00",
    "Patient4" = "#56B4E9",
    "Patient6" = "#009E73"
  )
)

###############################################################################
# 9. PLOT HEATMAP
###############################################################################
pheatmap(
  heatmap_mat_z,
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100),
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  scale = "none",
  show_colnames = FALSE,
  fontsize_row = 9,
  border_color = NA,
  gaps_row = gap_positions,
  main = "Lineage Marker Expression (VST, Z-scored per gene)",
  filename = "Lineage_marker_expression_heatmap.pdf",
  width = 8,
  height = 10
)


















###############################################################################
# DOING STOMACH AND INTESTINAL SIGNATURES
###############################################################################
###############################################################################
# STOMACH AND INTESTINAL SIGNATURE SCORE ANALYSIS
# Based on Human Protein Atlas organ-specific gene signatures
###############################################################################

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Color schemes
color_palette <- c("Gastric" = "#C75A5A", "IM" = "#9370B8", "Duodenal" = "#6B8DB8")
patient_palette <- c("Patient3" = "#E69F00", "Patient4" = "#56B4E9", "Patient6" = "#009E73")

cat("=== Calculating Organ Signature Scores ===\n\n")

# ============================================================================
# ORGAN-SPECIFIC GENE SIGNATURES FROM HUMAN PROTEIN ATLAS
# ============================================================================

stomach_signature <- c(
  "PGA4", "GIF", "GKN1", "PGA3", "PGA5", "ATP4A", "GAST", "ATP4B", "LIPF", "CLDN18",
  "GAGE12G", "MUC5AC", "GAGE2E", "DAZ4", "GKN2", "TFF1", "TFF2", "DAZ2", "KCNE2",
  "NKX6-2", "MUC6", "C6orf58", "PGC", "TAAR1", "CCKBR", "DPCR1", "FUT9", "CHIA",
  "CYP2C18", "CCKAR", "VSIG1", "COL2A1", "CTSE", "NTS", "A4GNT", "ANXA10", "AQP5", "CAPN8",
  "FAM159B", "GHRL", "PSCA", "ALDH3A1", "BARX1", "EPS8L3", "ERN2", "ESRRB", "PAX6",
  "PSAPL1", "C2orf70", "HNF4A", "IGFALS", "IHH", "SDR16C5", "SOX21", "TPH1", "TRIM31",
  "ACER2", "AKR1C2", "AKR7A3", "FER1L6", "FOXA3", "ISL1", "LRRC31", "MYEOV", "NPAS1",
  "RFLNA", "SLC26A9", "SLC9A4", "SULT1B1", "ABC13-47488600E17.1", "ACSM6",
  "ADH1C", "AGR2", "AKR1B10", "AMTN", "ANKRD22", "ANO7", "ARL14", "B3GNT6",
  "B4GALNT3", "BCAS1", "BPIFB1", "CA2", "CA9", "CAPN13", "CAPN9", "CFC1", "CLDN23",
  "CLIC6", "CRYBA2", "CXCL17", "CYP2S1", "DRD5", "EPS8L1", "FA2H", "FAM177B",
  "FAM83E", "FEV", "FOXA2", "FOXQ1", "FRMD1", "GALNT5", "GALNT6", "GIPR", "GPR25",
  "GUCA1C", "HAP1", "HTR1B", "INSM1", "KCNK16", "KPNA7", "LA16c-312E8.5",
  "LGALS9B", "LGALS9C", "LIME1", "LINC00675", "LIPH", "MBOAT4", "MFSD4A", "MIA",
  "MUC1", "MYRF", "NEUROD1", "NKX2-2", "NKX6-3", "NMUR2", "NPW", "NQO1", "OASL",
  "ONECUT3", "OVOL2", "PIK3C2G", "PLA2G10", "PRSS22", "PTGDR2", "RAB27B",
  "RASSF6", "REP15", "RFX6", "RNF223", "RP11-599B13.6", "S100P", "SCGN", "SHH",
  "SLC5A5", "SLC9A2", "SMIM6", "SOSTDC1", "SPTSSB", "SST", "SSTR1", "SULT1C2",
  "SYTL2", "SYTL5", "TESC", "TMED6", "TMEM211", "TMEM238", "TRIM50", "TRIM74",
  "TRNP1", "UNC5CL", "UPK1B", "VILL", "VSIG2", "ZSCAN4"
)

intestinal_signature <- c(
  "INSL5", "AQP8", "MEP1A", "PRAC1", "ISX", "CA1", "TBX10", "BTNL3", "MUC12", "TMIGD1",
  "KRTAP13-2", "MS4A12", "NOX1", "CD177", "GUCA2A", "PYY", "CHST5", "GUCY2C", "MYO1A",
  "NAT2", "NXPE1", "SDHD", "LYPD8", "MOGAT2", "SLC39A5", "AC011513.3", "C10orf99",
  "PIGY", "CHP2", "GLRA2", "HHLA2", "NR1I2", "ATOH1", "CDHR5", "CDX1", "LRRC19", "NXPE4",
  "REG4", "TMEM236", "CDH17", "CLCA1", "CLRN3", "DHRS11", "GAL3ST2", "GPA33", "HOXD13",
  "CEACAM7", "EPS8L3", "ERN2", "FABP1", "FSIP1", "GALNT8", "LRRC26", "PHGR1", "PPP1R14D",
  "SH2D7", "VIL1", "AMN", "BTNL8", "CASP5", "CDX2", "CEACAM6", "CLCA4", "EFNA2", "HAVCR1",
  "HNF4A", "IHH", "LINC01207", "MISP", "MOGAT3", "MUC13", "NEU4", "SLC17A4", "SPINK4",
  "TPH1", "TRIM31", "URAD", "ZG16", "B3GALT5", "FOXA3", "IL22RA1", "LRRC31", "NOS2", "SATB2",
  "SULT1B1", "TRPM5", "VIP", "AC009133.22", "ADAMDEC1", "AIFM3", "ATP10B", "B3GALT1",
  "B3GNT6", "B4GALNT2", "BCL2L15", "BEST2", "BEST4", "C15orf48", "C2orf72", "CA2", "CA4",
  "CA7", "CDC42EP5", "CEACAM1", "CEACAM5", "CES3", "CH17-360D5.1", "CKMT1B", "CLDN23",
  "CLDN3", "CTC-273B12.7", "DHRS9", "ENTPD8", "FAM3D", "FCGBP", "FFAR4", "FOXD2", "FRMD1",
  "FUT3", "FXYD3", "GCNT3", "GPR15", "GUCA2B", "HEPACAM2", "HOXD12", "HSD11B2", "ITLN1",
  "KLK15", "KRT20", "LEFTY1", "LGALS4", "LINC00675", "MAB21L2", "MUC4", "NOXO1", "NPY4R",
  "NRARP", "NXPE2", "OTOP2", "PADI2", "PIGR", "PIGZ", "PKIB", "PLA2G10", "PLA2G2A", "REP15",
  "RETNLB", "RP11-599B13.6", "RXFP4", "SLC22A18AS", "SLC26A2", "SLC26A3", "SLC9A2",
  "SLC9A3", "ST6GALNAC1", "TFF3", "TMEM171", "TPSG1", "TRABD2A", "TRIM15", "TRPM6",
  "TSPAN1", "TSPAN8", "UGT1A10", "UGT1A8", "UGT2B1"
)

cat("Stomach signature genes:", length(stomach_signature), "\n")
cat("Intestinal signature genes:", length(intestinal_signature), "\n\n")

# ============================================================================
# LOAD VST-TRANSFORMED DATA (log2-like scale similar to log2 TPM)
# ============================================================================

# Load VST data from the existing analysis
# VST provides log2-like transformation similar to log2 TPM used in the paper
vst_mat <- assay(vsd)  # This should be from your existing script

# Load metadata
merged_counts_file <- "setworkingdirectory/salmon.merged.gene_counts.tsv"
sample_sheet <- "setworkingdirectory/sample_metadata.csv"

colData <- read.csv(sample_sheet, row.names=1)

# Fix sample name
if("Im13" %in% rownames(colData)) {
  rownames(colData)[rownames(colData) == "Im13"] <- "IM13"
}

# Prepare metadata
colData$phenotype <- factor(colData$phenotype, levels = c("Gastric", "IM", "Duodenal"))
colData$Patient <- factor(case_when(
  grepl("3", rownames(colData)) ~ "Patient3",
  grepl("4", rownames(colData)) ~ "Patient4",
  grepl("6", rownames(colData)) ~ "Patient6"
))

# ============================================================================
# CHECK GENE AVAILABILITY
# ============================================================================

stomach_genes_present <- stomach_signature[stomach_signature %in% rownames(vst_mat)]
intestinal_genes_present <- intestinal_signature[intestinal_signature %in% rownames(vst_mat)]

cat("Stomach signature genes found in data:", length(stomach_genes_present), 
    "/", length(stomach_signature), 
    "(", round(100*length(stomach_genes_present)/length(stomach_signature), 1), "%)\n")

cat("Intestinal signature genes found in data:", length(intestinal_genes_present), 
    "/", length(intestinal_signature),
    "(", round(100*length(intestinal_genes_present)/length(intestinal_signature), 1), "%)\n\n")

# Report missing genes
stomach_missing <- setdiff(stomach_signature, stomach_genes_present)
intestinal_missing <- setdiff(intestinal_signature, intestinal_genes_present)

if(length(stomach_missing) > 0) {
  cat("Missing stomach genes:", paste(head(stomach_missing, 10), collapse=", "), "...\n")
}
if(length(intestinal_missing) > 0) {
  cat("Missing intestinal genes:", paste(head(intestinal_missing, 10), collapse=", "), "...\n\n")
}

# ============================================================================
# CALCULATE SIGNATURE SCORES
# Following the paper's method exactly
# ============================================================================

cat("=== Calculating Signature Scores ===\n\n")

# Extract expression for signature genes
stomach_expr <- vst_mat[stomach_genes_present, , drop=FALSE]
intestinal_expr <- vst_mat[intestinal_genes_present, , drop=FALSE]

# Calculate scores as mean of log2-like values (VST)
stomach_score <- colMeans(stomach_expr, na.rm=TRUE)
intestinal_score <- colMeans(intestinal_expr, na.rm=TRUE)

# Calculate aggregate score
# Stomach genes multiplied by -1, then averaged with intestinal genes
stomach_expr_negative <- -1 * stomach_expr
combined_expr <- rbind(stomach_expr_negative, intestinal_expr)
aggregate_score <- colMeans(combined_expr, na.rm=TRUE)

# Create results dataframe
signature_scores <- data.frame(
  Sample = names(stomach_score),
  Phenotype = colData$phenotype,
  Patient = colData$Patient,
  Stomach_Score = stomach_score,
  Intestinal_Score = intestinal_score,
  Aggregate_Score = aggregate_score
)

# Save results
write.csv(signature_scores, "Signature_Scores.csv", row.names=FALSE)

cat("Signature scores calculated and saved.\n\n")
print(signature_scores)

# ============================================================================
# VISUALIZATIONS
# ============================================================================

cat("\n=== Creating Visualizations ===\n\n")

# 1. Scatter plot: Stomach vs Intestinal Score
p1 <- ggplot(signature_scores, aes(x=Stomach_Score, y=Intestinal_Score, 
                                   color=Phenotype, shape=Patient)) +
  geom_point(size=5, alpha=0.8) +
  scale_color_manual(values=color_palette) +
  scale_shape_manual(values=c(15,16,17)) +
  geom_vline(xintercept=median(stomach_score), linetype="dashed", alpha=0.5) +
  geom_hline(yintercept=median(intestinal_score), linetype="dashed", alpha=0.5) +
  theme_bw(base_size=14) +
  labs(title="Stomach vs Intestinal Signature Scores",
       x="Stomach Signature Score",
       y="Intestinal Signature Score") +
  theme(plot.title = element_text(hjust=0.5, face="bold"))

ggsave("Signature_Scores_Scatter.pdf", p1, width=10, height=8)

# 2. Boxplots for each score by phenotype
scores_long <- signature_scores %>%
  pivot_longer(cols=c(Stomach_Score, Intestinal_Score, Aggregate_Score),
               names_to="Score_Type", values_to="Score") %>%
  mutate(Score_Type = factor(Score_Type, 
                             levels=c("Stomach_Score", "Intestinal_Score", "Aggregate_Score"),
                             labels=c("Stomach", "Intestinal", "Aggregate")))

p2 <- ggplot(scores_long, aes(x=Phenotype, y=Score, fill=Phenotype)) +
  geom_boxplot(alpha=0.7, outlier.shape=NA) +
  geom_jitter(aes(shape=Patient), width=0.15, size=3) +
  facet_wrap(~Score_Type, scales="free_y", nrow=1) +
  scale_fill_manual(values=color_palette) +
  scale_shape_manual(values=c(15,16,17)) +
  theme_bw(base_size=14) +
  labs(title="Signature Scores Across Phenotypes",
       y="Score (Mean VST Expression)") +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        axis.text.x = element_text(angle=45, hjust=1))

ggsave("Signature_Scores_Boxplots.pdf", p2, width=14, height=6)

# 3. Heatmap of signature gene expression
# Combine both signatures
all_sig_genes <- c(stomach_genes_present, intestinal_genes_present)
sig_expr <- vst_mat[all_sig_genes, ]

# Create annotation
gene_annotation <- data.frame(
  Signature = c(rep("Stomach", length(stomach_genes_present)),
                rep("Intestinal", length(intestinal_genes_present)))
)
rownames(gene_annotation) <- all_sig_genes

sample_annotation <- data.frame(
  Phenotype = colData$phenotype,
  Patient = colData$Patient
)
rownames(sample_annotation) <- colnames(sig_expr)

annotation_colors <- list(
  Phenotype = color_palette,
  Patient = patient_palette,
  Signature = c("Stomach"="#C75A5A", "Intestinal"="#6B8DB8")
)

# Z-score for visualization
sig_expr_scaled <- t(scale(t(sig_expr)))

pdf("Signature_Genes_Heatmap.pdf", width=10, height=16)
pheatmap(sig_expr_scaled,
         annotation_col = sample_annotation,
         annotation_row = gene_annotation,
         annotation_colors = annotation_colors,
         show_rownames = FALSE,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-2, 2, length.out=101),
         main = "Organ Signature Gene Expression")
dev.off()

# 4. Aggregate score progression
p3 <- ggplot(signature_scores, aes(x=Phenotype, y=Aggregate_Score, fill=Phenotype)) +
  geom_boxplot(alpha=0.7, outlier.shape=NA, width=0.6) +
  geom_jitter(aes(shape=Patient), width=0.15, size=4) +
  geom_hline(yintercept=0, linetype="dashed", color="black") +
  scale_fill_manual(values=color_palette) +
  scale_shape_manual(values=c(15,16,17)) +
  theme_bw(base_size=16) +
  labs(title="Aggregate Score: Stomach ← → Intestinal",
       subtitle="Negative = Stomach-like, Positive = Intestinal-like",
       y="Aggregate Score") +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5),
        legend.position="right")

ggsave("Aggregate_Score_Boxplot.pdf", p3, width=10, height=8)

# 5. Individual sample trajectories
p4 <- ggplot(signature_scores, aes(x=Phenotype, y=Aggregate_Score, group=Patient, color=Patient)) +
  geom_line(linewidth=1, alpha=0.7) +
  geom_point(size=4, alpha=0.9) +
  scale_color_manual(values=patient_palette) +
  geom_hline(yintercept=0, linetype="dashed", color="black") +
  theme_bw(base_size=14) +
  labs(title="Patient-Specific Trajectories: Aggregate Score",
       subtitle="Stomach-like → Intestinal-like progression",
       y="Aggregate Score") +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        plot.subtitle = element_text(hjust=0.5))

ggsave("Aggregate_Score_Trajectories.pdf", p4, width=10, height=7)

# ============================================================================
# STATISTICAL ANALYSIS
# ============================================================================

cat("\n=== Statistical Analysis ===\n\n")

# ANOVA for each score type
anova_stomach <- aov(Stomach_Score ~ Phenotype, data=signature_scores)
anova_intestinal <- aov(Intestinal_Score ~ Phenotype, data=signature_scores)
anova_aggregate <- aov(Aggregate_Score ~ Phenotype, data=signature_scores)

cat("ANOVA Results:\n\n")
cat("Stomach Score:\n")
print(summary(anova_stomach))
cat("\nIntestinal Score:\n")
print(summary(anova_intestinal))
cat("\nAggregate Score:\n")
print(summary(anova_aggregate))

# Pairwise comparisons
tukey_stomach <- TukeyHSD(anova_stomach)
tukey_intestinal <- TukeyHSD(anova_intestinal)
tukey_aggregate <- TukeyHSD(anova_aggregate)

# Create summary table
stat_results <- data.frame(
  Score = rep(c("Stomach", "Intestinal", "Aggregate"), each=3),
  Comparison = rep(c("IM-Gastric", "Duodenal-Gastric", "Duodenal-IM"), 3),
  Diff = c(tukey_stomach$Phenotype[,"diff"],
           tukey_intestinal$Phenotype[,"diff"],
           tukey_aggregate$Phenotype[,"diff"]),
  padj = c(tukey_stomach$Phenotype[,"p adj"],
           tukey_intestinal$Phenotype[,"p adj"],
           tukey_aggregate$Phenotype[,"p adj"])
)

write.csv(stat_results, "Signature_Scores_Statistics.csv", row.names=FALSE)

cat("\n")
print(stat_results)

# ============================================================================
# SUMMARY REPORT
# ============================================================================

cat("\n========================================\n")
cat("SIGNATURE SCORE ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat("Files created:\n")
cat("  - Signature_Scores.csv (all scores for each sample)\n")
cat("  - Signature_Scores_Scatter.pdf\n")
cat("  - Signature_Scores_Boxplots.pdf\n")
cat("  - Signature_Genes_Heatmap.pdf\n")
cat("  - Aggregate_Score_Boxplot.pdf\n")
cat("  - Aggregate_Score_Trajectories.pdf\n")
cat("  - Signature_Scores_Statistics.csv\n\n")

cat("Key Findings:\n")
cat("  - Stomach signature genes used:", length(stomach_genes_present), "\n")
cat("  - Intestinal signature genes used:", length(intestinal_genes_present), "\n")
cat("  - Aggregate score represents progression from stomach to intestinal phenotype\n\n")

cat("Analysis complete!\n")


