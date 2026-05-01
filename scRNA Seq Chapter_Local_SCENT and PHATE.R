# Load packages
library(Seurat)
library(SeuratObject)
library(harmony)
library(ggplot2)
library(dplyr)
library(EnhancedVolcano)
library(pheatmap)

# Set working directory
setwd("setdirectory")

set.seed(1234)

# Load matrices
mtx.6Gas <- Read10X_h5(filename = "./CO_06_GAS_GEX/filtered_feature_bc_matrix.h5")
mtx.6IM1 <- Read10X_h5(filename = "./CO_06_IM1_GEX/filtered_feature_bc_matrix.h5")
mtx.6IM2 <- Read10X_h5(filename = "./CO_06_IM2_GEX/filtered_feature_bc_matrix.h5")
mtx.6Duo <- Read10X_h5(filename = "./CO_06_DUO_GEX/filtered_feature_bc_matrix.h5")
mtx.7Gas <- Read10X_h5(filename = "./CO_07_GAS_GEX/filtered_feature_bc_matrix.h5")
mtx.7IM1 <- Read10X_h5(filename = "./CO_07_IM1_GEX/filtered_feature_bc_matrix.h5")
mtx.7IM2 <- Read10X_h5(filename = "./CO_07_IM2_GEX/filtered_feature_bc_matrix.h5")
mtx.7Duo <- Read10X_h5(filename = "./CO_07_DUO_GEX/filtered_feature_bc_matrix.h5")

# Create Seurat Objects
seu.6Gas <- CreateSeuratObject(mtx.6Gas, project = "CO_06_GAS")
seu.6IM1 <- CreateSeuratObject(mtx.6IM1, project = "CO_06_IM1")
seu.6IM2 <- CreateSeuratObject(mtx.6IM2, project = "CO_06_IM2")
seu.6Duo <- CreateSeuratObject(mtx.6Duo, project = "CO_06_DUO")
seu.7Gas <- CreateSeuratObject(mtx.7Gas, project = "CO_07_GAS")
seu.7IM1 <- CreateSeuratObject(mtx.7IM1, project = "CO_07_IM1")
seu.7IM2 <- CreateSeuratObject(mtx.7IM2, project = "CO_07_IM2")
seu.7Duo <- CreateSeuratObject(mtx.7Duo, project = "CO_07_DUO")

# Calculate mitochondrial %
seu.6Gas <- PercentageFeatureSet(seu.6Gas, pattern = "^MT-", col.name = "percent.mito")
seu.6IM1 <- PercentageFeatureSet(seu.6IM1, pattern = "^MT-", col.name = "percent.mito")
seu.6IM2 <- PercentageFeatureSet(seu.6IM2, pattern = "^MT-", col.name = "percent.mito")
seu.6Duo <- PercentageFeatureSet(seu.6Duo, pattern = "^MT-", col.name = "percent.mito")
seu.7Gas <- PercentageFeatureSet(seu.7Gas, pattern = "^MT-", col.name = "percent.mito")
seu.7IM1 <- PercentageFeatureSet(seu.7IM1, pattern = "^MT-", col.name = "percent.mito")
seu.7IM2 <- PercentageFeatureSet(seu.7IM2, pattern = "^MT-", col.name = "percent.mito")
seu.7Duo <- PercentageFeatureSet(seu.7Duo, pattern = "^MT-", col.name = "percent.mito")

# Filter cells
seu.6Gas <- subset(seu.6Gas, nFeature_RNA > 400 & percent.mito < 10)
seu.6IM1 <- subset(seu.6IM1, nFeature_RNA > 400 & percent.mito < 10)
seu.6IM2 <- subset(seu.6IM2, nFeature_RNA > 400 & percent.mito < 10)
seu.6Duo <- subset(seu.6Duo, nFeature_RNA > 400 & percent.mito < 10)
seu.7Gas <- subset(seu.7Gas, nFeature_RNA > 400 & percent.mito < 10)
seu.7IM1 <- subset(seu.7IM1, nFeature_RNA > 400 & percent.mito < 10)
seu.7IM2 <- subset(seu.7IM2, nFeature_RNA > 400 & percent.mito < 10)
seu.7Duo <- subset(seu.7Duo, nFeature_RNA > 400 & percent.mito < 10)

# Merge samples
seu.int <- merge(seu.6Gas, y = c(seu.6IM1, seu.6IM2, seu.6Duo, seu.7Gas, seu.7IM1, seu.7IM2, seu.7Duo))

# Normalize and find variable features BEFORE Harmony
seu.int <- NormalizeData(seu.int)
seu.int <- FindVariableFeatures(seu.int)

# Scale data WITHOUT regression (we regress later after Harmony)
seu.int <- ScaleData(seu.int, features = VariableFeatures(seu.int))

# Run PCA
seu.int <- RunPCA(seu.int, features = VariableFeatures(seu.int))

# Add metadata for patient and tissue info
seu.int$patient <- sapply(strsplit(seu.int$orig.ident, "_"), function(x) paste(x[1], x[2], sep = "_"))
seu.int$tissue_type <- sapply(strsplit(seu.int$orig.ident, "_"), function(x) x[3])
seu.int$tissue_group <- ifelse(seu.int$tissue_type %in% c("IM1", "IM2"), "IM", seu.int$tissue_type)

# Run Harmony for batch correction by patient
seu.int <- RunHarmony(seu.int, group.by.vars = "patient")

# Cell cycle scoring BEFORE regression (using default Seurat gene lists)
seu.int <- CellCycleScoring(
  object = seu.int,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident = TRUE
)

# Regress out cell cycle scores during scaling using Harmony embeddings
seu.int <- ScaleData(
  seu.int,
  vars.to.regress = c("S.Score", "G2M.Score"),
  features = VariableFeatures(seu.int)
)

# Run UMAP and clustering on Harmony embeddings
seu.int <- RunUMAP(seu.int, reduction = "harmony", dims = 1:30)
seu.int <- FindNeighbors(seu.int, reduction = "harmony", dims = 1:30)
seu.int <- FindClusters(seu.int, resolution = 2)

# Visualize clusters and cell cycle phase
DimPlot(seu.int, reduction = "umap", label = TRUE, pt.size = 0.5) + ggtitle("Clusters (Harmony Batch Corrected)")
DimPlot(seu.int, reduction = "umap", group.by = "Phase", pt.size = 0.5) + ggtitle("Cell Cycle Phase")
DimPlot(seu.int, reduction = "umap", group.by = "patient", pt.size = 0.5) + ggtitle("Patient")
DimPlot(seu.int, reduction = "umap", group.by = "tissue_group", pt.size = 0.5) + ggtitle("Organoid Identity")



# CCAT Analysis Script for Seurat Objects
# This script calculates CCAT scores for differentiation potency analysis
# and creates visualization plots

# Load required libraries
library(Seurat)
library(SCENT)
library(org.Hs.eg.db)

# ============================================================================
# STEP 1: PREPARE EXPRESSION DATA FROM SEURAT OBJECT
# ============================================================================

# Join all layers in the RNA assay (required for Seurat v5)
seu.int[["RNA"]] <- JoinLayers(seu.int[["RNA"]])

# Extract expression data (keep as sparse matrix to save memory)
expression_matrix <- LayerData(seu.int, assay = "RNA", layer = "data")

# Check dimensions
cat("Expression matrix dimensions:\n")
cat("Genes:", nrow(expression_matrix), "\n")
cat("Cells:", ncol(expression_matrix), "\n\n")

# ============================================================================
# STEP 2: CONVERT GENE IDENTIFIERS TO ENTREZ IDs
# ============================================================================

# Get current gene names
current_genes <- rownames(expression_matrix)

# Separate gene symbols from Ensembl IDs
gene_symbols <- current_genes[!grepl("^ENSG", current_genes)]
ensembl_genes <- current_genes[grepl("^ENSG", current_genes)]

cat("Gene identifier breakdown:\n")
cat("Gene symbols:", length(gene_symbols), "\n")
cat("Ensembl IDs:", length(ensembl_genes), "\n\n")

# Convert gene symbols to Entrez IDs using org.Hs.eg.db
symbol_to_entrez <- mapIds(org.Hs.eg.db, 
                           keys = gene_symbols,
                           column = "ENTREZID", 
                           keytype = "SYMBOL",
                           multiVals = "first")

# Remove NAs
symbol_to_entrez <- symbol_to_entrez[!is.na(symbol_to_entrez)]

cat("Conversion results:\n")
cat("Gene symbols successfully converted:", length(symbol_to_entrez), "\n")

# Create filtered expression matrix with converted gene IDs
genes_to_keep <- intersect(gene_symbols, names(symbol_to_entrez))
expression_matrix_filtered <- expression_matrix[genes_to_keep, ]

# Convert rownames to Entrez IDs
rownames(expression_matrix_filtered) <- symbol_to_entrez[genes_to_keep]

# ============================================================================
# STEP 3: CHECK OVERLAP WITH PPI NETWORK AND CALCULATE CCAT
# ============================================================================

# Check overlap with PPI network
overlap_genes <- intersect(rownames(expression_matrix_filtered), rownames(net13Jun12.m))
cat("Overlap with PPI network:", length(overlap_genes), "\n")

# Verify we have enough genes (SCENT recommends >5000)
if (length(overlap_genes) < 5000) {
  stop("Not enough overlapping genes for reliable CCAT calculation. Need >5000 genes.")
} else {
  cat("✓ Sufficient gene overlap for CCAT calculation\n\n")
}

# Calculate CCAT scores
cat("Calculating CCAT scores...\n")
ccat_scores <- CompCCAT(exp = expression_matrix_filtered, ppiA = net13Jun12.m)

# Check results
cat("CCAT calculation complete!\n")
cat("Number of scores:", length(ccat_scores), "\n")
cat("Score range:", round(range(ccat_scores, na.rm = TRUE), 4), "\n\n")

# ============================================================================
# STEP 4: ADD CCAT SCORES TO SEURAT METADATA
# ============================================================================

# Add CCAT scores to Seurat object metadata
seu.int$CCAT <- ccat_scores

# Verify assignment
cat("CCAT scores added to metadata:\n")
cat("Summary statistics:\n")
print(summary(seu.int$CCAT))
cat("\n")

####Convert into a filetype that can be used in Python
#if (!requireNamespace("remotes", quietly = TRUE)) {
#  install.packages("remotes")
#}
#remotes::install_github("mojaveazure/seurat-disk")

#library(SeuratDisk)
#SaveH5Seurat(seu.int, filename = "seu_int.h5Seurat")
#Convert("seu_int.h5Seurat", dest = "h5ad")

# Extract the neighbor graph (usually shared nearest neighbor)
# Extract the SNN graph
library(Matrix)

# Extract sparse graph
nn_graph <- seu.int@graphs$RNA_snn

# Write to disk in Matrix Market format
Matrix::writeMM(nn_graph, file = "seurat_RNA_snn_graph.mtx")

# Also save cell names
writeLines(colnames(nn_graph), "seurat_RNA_snn_cells.txt")


# Load the PHATE coordinates
phate_coords <- read.csv("PHATE_coordinates.csv", row.names = 1)

# Add to your Seurat object
seu.int[["phate"]] <- CreateDimReducObject(
  embeddings = as.matrix(phate_coords),
  key = "PHATE_",
  assay = "RNA"
)

# Visualize!
DimPlot(seu.int, reduction = "phate")
DimPlot(seu.int, reduction = "phate", group.by = "tissue_group")
DimPlot(seu.int, reduction = "phate", group.by = "patient")

library(ggplot2)

# Extract data
phate_df <- as.data.frame(Embeddings(seu.int, "phate"))
phate_df$CCAT <- seu.int$CCAT

# Method 1: Using ggplot (best control)
ggplot(phate_df, aes(x = PHATE_1, y = PHATE_2, color = CCAT)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_viridis_c(option = "plasma", name = "CCAT") +
  theme_classic() +
  ggtitle("CCAT Score on PHATE") +
  coord_fixed()  




# ============================================================================
# STEP 5: CREATE VISUALIZATION PLOTS
# ============================================================================

cat("Creating visualization plots...\n")

# 1. Boxplot by Original Identity (sample)
boxplot(seu.int$CCAT ~ seu.int$orig.ident,
        main = "CCAT Potency Estimates by Sample",
        xlab = "Sample",
        ylab = "CCAT Score",
        col = rainbow(length(unique(seu.int$orig.ident))),
        las = 2,
        cex.axis = 0.8)

cat("Sample sizes:\n")
print(table(seu.int$orig.ident))
cat("\n")

# 2. Boxplot by Patient
boxplot(seu.int$CCAT ~ seu.int$patient,
        main = "CCAT Potency Estimates by Patient",
        xlab = "Patient",
        ylab = "CCAT Score",
        col = c("lightblue", "lightcoral"))

cat("Patient sizes:\n")
print(table(seu.int$patient))
cat("\n")

# 3. Boxplot by Tissue Type
boxplot(seu.int$CCAT ~ seu.int$tissue_group,
        main = "CCAT Potency Estimates by Tissue Type",
        xlab = "Tissue Type",
        ylab = "CCAT Score",
        col = c("lightgreen", "orange", "purple", "yellow"),
        las = 2)

cat("Tissue type sizes:\n")
print(table(seu.int$tissue_type))
cat("\n")

# 4. Combined multi-panel plot
par(mfrow = c(1, 3), mar = c(8, 4, 4, 2))

# By sample
boxplot(seu.int$CCAT ~ seu.int$orig.ident,
        main = "By Sample",
        xlab = "",
        ylab = "CCAT Score",
        col = rainbow(8),
        las = 2,
        cex.axis = 0.7)

# By patient
boxplot(seu.int$CCAT ~ seu.int$patient,
        main = "By Patient",
        xlab = "",
        ylab = "CCAT Score",
        col = c("lightblue", "lightcoral"),
        las = 2)

# By tissue type
boxplot(seu.int$CCAT ~ seu.int$tissue_type,
        main = "By Tissue Type",
        xlab = "",
        ylab = "CCAT Score",
        col = c("lightgreen", "orange", "purple", "yellow"),
        las = 2,
        cex.axis = 0.8)

# Reset plotting parameters
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

cat("✓ CCAT analysis complete!\n")
cat("CCAT scores have been added to seu.int$CCAT\n")
cat("Higher CCAT values indicate higher differentiation potency\n")

###More plots with CCAT
# ============================================================================
# 1. BOXPLOT: CCAT by Cluster
# ============================================================================

# Basic boxplot
p1 <- ggplot(seu.int@meta.data, aes(x = seurat_clusters, y = CCAT, fill = seurat_clusters)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  labs(title = "CCAT Scores by Cluster",
       x = "Cluster",
       y = "CCAT Score (Differentiation Potency)") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis(discrete = TRUE, option = "turbo")

print(p1)

# ============================================================================
# 2. VIOLIN PLOT: CCAT by Cluster
# ============================================================================

p2 <- ggplot(seu.int@meta.data, aes(x = seurat_clusters, y = CCAT, fill = seurat_clusters)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3) +
  theme_classic() +
  labs(title = "CCAT Score Distribution by Cluster",
       x = "Cluster",
       y = "CCAT Score") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis(discrete = TRUE, option = "turbo")

print(p2)

# ============================================================================
# 3. UMAP COLOURED BY CCAT
# ============================================================================

p3 <- FeaturePlot(seu.int, features = "CCAT", pt.size = 0.5) +
  scale_color_viridis(option = "magma") +
  labs(title = "CCAT Scores on UMAP") +
  theme(plot.title = element_text(hjust = 0.5))

print(p3)

# ============================================================================
# 4. UMAP WITH CLUSTER LABELS + CCAT OVERLAY
# ============================================================================

# Side by side: clusters and CCAT
p4a <- DimPlot(seu.int, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("Clusters") +
  theme(plot.title = element_text(hjust = 0.5))

p4b <- FeaturePlot(seu.int, features = "CCAT", pt.size = 0.5) +
  scale_color_viridis(option = "magma") +
  ggtitle("CCAT Scores") +
  theme(plot.title = element_text(hjust = 0.5))

library(patchwork)
p4 <- p4a | p4b
print(p4)

# ============================================================================
# 5. BAR PLOT: MEAN CCAT PER CLUSTER
# ============================================================================

# Calculate mean CCAT per cluster
mean_ccat <- seu.int@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(mean_CCAT = mean(CCAT, na.rm = TRUE),
            se_CCAT = sd(CCAT, na.rm = TRUE) / sqrt(n()))

p5 <- ggplot(mean_ccat, aes(x = seurat_clusters, y = mean_CCAT, fill = seurat_clusters)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean_CCAT - se_CCAT, ymax = mean_CCAT + se_CCAT), 
                width = 0.2) +
  theme_classic() +
  labs(title = "Mean CCAT Score by Cluster",
       x = "Cluster",
       y = "Mean CCAT Score (± SE)") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis(discrete = TRUE, option = "turbo")

print(p5)

# ============================================================================
# 6. STATISTICAL SUMMARY TABLE
# ============================================================================

# Create summary statistics table
ccat_summary <- seu.int@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_cells = n(),
    mean_CCAT = mean(CCAT, na.rm = TRUE),
    median_CCAT = median(CCAT, na.rm = TRUE),
    sd_CCAT = sd(CCAT, na.rm = TRUE),
    min_CCAT = min(CCAT, na.rm = TRUE),
    max_CCAT = max(CCAT, na.rm = TRUE)
  ) %>%
  arrange(desc(mean_CCAT))

print("CCAT Summary by Cluster (sorted by mean CCAT):")
print(ccat_summary)

# ============================================================================
# 7. HEATMAP: CCAT VALUES
# ============================================================================

# Create a matrix of CCAT values by cluster
library(pheatmap)

# Sample cells for visualization (if you have too many)
set.seed(123)
if(ncol(seu.int) > 5000) {
  cells_to_plot <- sample(colnames(seu.int), 5000)
} else {
  cells_to_plot <- colnames(seu.int)
}

# Prepare data for heatmap
heatmap_data <- seu.int@meta.data[cells_to_plot, c("seurat_clusters", "CCAT")]
heatmap_data <- heatmap_data[order(heatmap_data$seurat_clusters, heatmap_data$CCAT), ]

# Create annotation
annotation_row <- data.frame(
  Cluster = heatmap_data$seurat_clusters,
  row.names = rownames(heatmap_data)
)

# Create matrix
ccat_matrix <- matrix(heatmap_data$CCAT, ncol = 1)
rownames(ccat_matrix) <- rownames(heatmap_data)
colnames(ccat_matrix) <- "CCAT"

pheatmap(ccat_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         annotation_row = annotation_row,
         color = viridis(100, option = "magma"),
         main = "CCAT Scores by Cell (Grouped by Cluster)")

# ============================================================================
# 8. DENSITY PLOT: CCAT DISTRIBUTION BY CLUSTER
# ============================================================================

p8 <- ggplot(seu.int@meta.data, aes(x = CCAT, fill = seurat_clusters)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  labs(title = "CCAT Score Density Distribution by Cluster",
       x = "CCAT Score",
       y = "Density",
       fill = "Cluster") +
  scale_fill_viridis(discrete = TRUE, option = "turbo")

print(p8)

# ============================================================================
# 9. IDENTIFY HIGH AND LOW POTENCY CLUSTERS
# ============================================================================

# Find clusters with highest and lowest mean CCAT
high_potency_clusters <- ccat_summary %>%
  slice_max(mean_CCAT, n = 3)

low_potency_clusters <- ccat_summary %>%
  slice_min(mean_CCAT, n = 3)

cat("\n========================================\n")
cat("TOP 3 HIGH POTENCY CLUSTERS:\n")
cat("========================================\n")
print(high_potency_clusters)

cat("\n========================================\n")
cat("TOP 3 LOW POTENCY CLUSTERS:\n")
cat("========================================\n")
print(low_potency_clusters)

###PLOTS WITH CCAT ON PHATE
library(ggplot2)

# Extract data
phate_df <- as.data.frame(Embeddings(seu.int, "phate"))
phate_df$CCAT <- seu.int$CCAT

# Method 1: Using ggplot (best control)
ggplot(phate_df, aes(x = PHATE_1, y = PHATE_2, color = CCAT)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_viridis_c(option = "plasma", name = "CCAT") +
  theme_classic() +
  ggtitle("CCAT Score on PHATE") +
  coord_fixed() 

# Method 2: Force FeaturePlot to match DimPlot dimensions
FeaturePlot(seu.int, 
            features = "CCAT", 
            reduction = "phate",
            pt.size = 0.5,
            order = TRUE) +  # This plots high values on top
  scale_color_viridis_c(option = "plasma") +
  coord_fixed() +
  theme(aspect.ratio = 1)

# Method 3: Create multiple plots with same dimensions
p1 <- DimPlot(seu.int, reduction = "phate", group.by = "tissue_group") + 
  ggtitle("Tissue Group") + coord_fixed()

p2 <- DimPlot(seu.int, reduction = "phate", group.by = "patient") + 
  ggtitle("Patient") + coord_fixed()

p3 <- ggplot(phate_df, aes(x = PHATE_1, y = PHATE_2, color = CCAT)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_viridis_c(option = "plasma") +
  theme_classic() +
  ggtitle("CCAT Score") +
  coord_fixed()

library(patchwork)
p1 | p2 | p3

################################################
# --- Compare strict vs flexible assignments ---
################################################

# Method 1: Side-by-side violin plots
library(patchwork)

p1 <- VlnPlot(seu.int.annot, 
              features = "CCAT", 
              group.by = "celltype_auto_strict",
              pt.size = 0) +
  ggtitle("CCAT by Strict Assignment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- VlnPlot(seu.int.annot, 
              features = "CCAT", 
              group.by = "celltype_auto_flexible",
              pt.size = 0) +
  ggtitle("CCAT by Flexible Assignment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1 / p2

# --- Method 2: Summary statistics for both ---
ccat_summary_strict <- seu.int.annot@meta.data %>%
  group_by(celltype_auto_strict) %>%
  summarise(
    n_cells = n(),
    mean_CCAT = mean(CCAT, na.rm = TRUE),
    median_CCAT = median(CCAT, na.rm = TRUE),
    sd_CCAT = sd(CCAT, na.rm = TRUE)
  ) %>%
  arrange(desc(median_CCAT)) %>%
  mutate(assignment = "Strict")

ccat_summary_flexible <- seu.int.annot@meta.data %>%
  group_by(celltype_auto_flexible) %>%
  summarise(
    n_cells = n(),
    mean_CCAT = mean(CCAT, na.rm = TRUE),
    median_CCAT = median(CCAT, na.rm = TRUE),
    sd_CCAT = sd(CCAT, na.rm = TRUE)
  ) %>%
  arrange(desc(median_CCAT)) %>%
  mutate(assignment = "Flexible")

# Combine and view
ccat_comparison <- rbind(ccat_summary_strict, ccat_summary_flexible)
print(ccat_comparison)

# Save
write.csv(ccat_comparison, "CCAT_strict_vs_flexible_comparison.csv", row.names = FALSE)

# --- Method 3: PHATE comparison ---
phate_df$celltype_strict <- seu.int.annot$celltype_auto_strict
phate_df$celltype_flexible <- seu.int.annot$celltype_auto_flexible

p3 <- DimPlot(seu.int.annot, reduction = "phate", 
              group.by = "celltype_auto_strict", 
              label = TRUE, repel = TRUE, pt.size = 0.3) +
  ggtitle("Strict Assignment")

p4 <- DimPlot(seu.int.annot, reduction = "phate", 
              group.by = "celltype_auto_flexible", 
              label = TRUE, repel = TRUE, pt.size = 0.3) +
  ggtitle("Flexible Assignment")

p5 <- ggplot(phate_df, aes(x = PHATE_1, y = PHATE_2, color = CCAT)) +
  geom_point(size = 0.5, alpha = 0.8) +
  scale_color_viridis_c(option = "plasma") +
  theme_classic() +
  coord_fixed() +
  ggtitle("CCAT Score")

(p3 | p4) / p5

# --- Method 4: Identify mixed/transitional states in flexible ---
seu.int.annot$is_mixed <- grepl("_", seu.int.annot$celltype_auto_flexible)

# Compare CCAT in pure vs mixed states
ggplot(seu.int.annot@meta.data, aes(x = is_mixed, y = CCAT, fill = is_mixed)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values = c("FALSE" = "lightblue", "TRUE" = "coral"),
                    labels = c("FALSE" = "Pure Identity", "TRUE" = "Mixed Identity")) +
  theme_classic() +
  labs(x = "", y = "CCAT Score", 
       title = "CCAT in Pure vs Mixed Cell States",
       fill = "Cell State") +
  scale_x_discrete(labels = c("FALSE" = "Pure", "TRUE" = "Mixed"))

# Statistical test
wilcox.test(CCAT ~ is_mixed, data = seu.int.annot@meta.data)

# --- Method 5: Focus on specific mixed states ---
mixed_states <- seu.int.annot@meta.data %>%
  filter(grepl("_", celltype_auto_flexible)) %>%
  group_by(celltype_auto_flexible) %>%
  summarise(
    n_cells = n(),
    mean_CCAT = mean(CCAT, na.rm = TRUE),
    median_CCAT = median(CCAT, na.rm = TRUE)
  ) %>%
  arrange(desc(n_cells))

print("Mixed/Transitional States:")
print(mixed_states)

# --- Method 6: Heatmap showing CCAT across both annotation types ---
library(pheatmap)

# Create matrix: rows = strict types, columns = flexible types
annotation_matrix <- seu.int.annot@meta.data %>%
  group_by(celltype_auto_strict, celltype_auto_flexible) %>%
  summarise(mean_CCAT = mean(CCAT, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = celltype_auto_flexible, 
              values_from = mean_CCAT)

# Convert to matrix for heatmap
mat <- as.matrix(annotation_matrix[,-1])
rownames(mat) <- annotation_matrix$celltype_auto_strict

pheatmap(mat, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         color = viridis::plasma(100),
         main = "Mean CCAT: Strict (rows) vs Flexible (cols)",
         na_col = "grey90")


# Detach plyr (it's interfering with dplyr)
detach("package:plyr", unload = TRUE)

# Reload dplyr
library(dplyr)

# Now this should work
mixed_states <- seu.int.annot@meta.data %>%
  filter(grepl("_", celltype_auto_flexible)) %>%
  group_by(celltype_auto_flexible) %>%
  summarise(
    n_cells = n(),
    mean_CCAT = mean(CCAT, na.rm = TRUE),
    median_CCAT = median(CCAT, na.rm = TRUE)
  ) %>%
  arrange(desc(n_cells))

print("Mixed/Transitional States:")
print(mixed_states)

# Method 1: Bar plot of cell counts coloured by CCAT
library(ggplot2)

ggplot(mixed_states, aes(x = reorder(celltype_auto_flexible, n_cells), 
                         y = n_cells, 
                         fill = median_CCAT)) +
  geom_col() +
  coord_flip() +
  scale_fill_viridis_c(option = "plasma", name = "Median\nCCAT") +
  theme_classic() +
  labs(x = "Mixed Cell State", 
       y = "Number of Cells",
       title = "Mixed/Transitional Cell States") +
  theme(axis.text.y = element_text(size = 10))

# Method 2: Scatter plot - cell count vs CCAT
ggplot(mixed_states, aes(x = n_cells, y = median_CCAT)) +
  geom_point(aes(size = n_cells, color = median_CCAT), alpha = 0.7) +
  geom_text(aes(label = celltype_auto_flexible), 
            hjust = -0.1, size = 3, check_overlap = TRUE) +
  scale_color_viridis_c(option = "plasma") +
  scale_size_continuous(range = c(3, 10)) +
  theme_classic() +
  labs(x = "Number of Cells", 
       y = "Median CCAT Score",
       title = "Mixed States: Abundance vs Potency") +
  theme(legend.position = "none")

# Method 3: Lollipop plot ordered by CCAT
ggplot(mixed_states, aes(x = reorder(celltype_auto_flexible, median_CCAT), 
                         y = median_CCAT)) +
  geom_segment(aes(xend = celltype_auto_flexible, y = 0, yend = median_CCAT),
               color = "grey50") +
  geom_point(aes(size = n_cells, color = median_CCAT)) +
  coord_flip() +
  scale_color_viridis_c(option = "plasma") +
  scale_size_continuous(range = c(3, 8), name = "# Cells") +
  theme_classic() +
  labs(x = "", 
       y = "Median CCAT Score",
       title = "Differentiation Potential in Mixed States",
       color = "CCAT")

# Method 4: Visualize these cells on PHATE
# Extract cells with mixed identity
mixed_cell_ids <- rownames(seu.int.annot@meta.data[grepl("_", seu.int.annot@meta.data$celltype_auto_flexible), ])

phate_df$is_mixed <- rownames(phate_df) %in% mixed_cell_ids
phate_df$mixed_type <- seu.int.annot$celltype_auto_flexible[match(rownames(phate_df), colnames(seu.int.annot))]

# Highlight mixed states
ggplot(phate_df, aes(x = PHATE_1, y = PHATE_2)) +
  geom_point(data = subset(phate_df, !is_mixed), 
             color = "grey80", size = 0.3, alpha = 0.3) +
  geom_point(data = subset(phate_df, is_mixed), 
             aes(color = mixed_type), size = 0.5, alpha = 0.7) +
  theme_classic() +
  coord_fixed() +
  labs(title = "Mixed Cell States on PHATE",
       color = "Mixed State") +
  theme(legend.position = "right")

# Method 5: Compare to pure states - create comparison dataframe
pure_states <- seu.int.annot@meta.data %>%
  filter(!grepl("_", celltype_auto_flexible)) %>%
  group_by(celltype_auto_flexible) %>%
  summarise(
    n_cells = n(),
    mean_CCAT = mean(CCAT, na.rm = TRUE),
    median_CCAT = median(CCAT, na.rm = TRUE)
  ) %>%
  mutate(state_type = "Pure")

mixed_states_labeled <- mixed_states %>%
  mutate(state_type = "Mixed")

combined_states <- bind_rows(pure_states, mixed_states_labeled)

# Plot comparison
ggplot(combined_states, aes(x = state_type, y = median_CCAT, fill = state_type)) +
  geom_violin(alpha = 0.6) +
  geom_jitter(aes(size = n_cells), width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = c("Pure" = "lightblue", "Mixed" = "coral")) +
  theme_classic() +
  labs(x = "Cell State Type", 
       y = "Median CCAT Score",
       title = "CCAT: Pure vs Mixed Cell States",
       size = "# Cells") +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black")

# Method 6: Just a clean table view
library(knitr)
kable(mixed_states, 
      col.names = c("Mixed Cell State", "# Cells", "Mean CCAT", "Median CCAT"),
      digits = 3,
      caption = "Mixed/Transitional Cell States Summary")


