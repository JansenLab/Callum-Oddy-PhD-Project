if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("anndataR")
BiocManager::install("hdf5r")
BiocManager::install("rhdf5")




#Preparing Object to be Exported for Cell Typist
# Load packages
library(Seurat)
library(SeuratObject)
library(hdf5r)
library(rhdf5)
library(harmony)
library(ggplot2)
library(dplyr)
library(EnhancedVolcano)
library(pheatmap)
library(anndataR)


# Set working directory
setwd("/Users/callu/OneDrive - University College London/Protocol/scRNA/")

set.seed(1234)

# Load matrices (your existing code)
mtx.6Gas <- Read10X_h5(filename = "./CO_06_GAS_GEX/filtered_feature_bc_matrix.h5")
mtx.6IM1 <- Read10X_h5(filename = "./CO_06_IM1_GEX/filtered_feature_bc_matrix.h5")
mtx.6IM2 <- Read10X_h5(filename = "./CO_06_IM2_GEX/filtered_feature_bc_matrix.h5")
mtx.6Duo <- Read10X_h5(filename = "./CO_06_DUO_GEX/filtered_feature_bc_matrix.h5")
mtx.7Gas <- Read10X_h5(filename = "./CO_07_GAS_GEX/filtered_feature_bc_matrix.h5")
mtx.7IM1 <- Read10X_h5(filename = "./CO_07_IM1_GEX/filtered_feature_bc_matrix.h5")
mtx.7IM2 <- Read10X_h5(filename = "./CO_07_IM2_GEX/filtered_feature_bc_matrix.h5")
mtx.7Duo <- Read10X_h5(filename = "./CO_07_DUO_GEX/filtered_feature_bc_matrix.h5")

# Create Seurat Objects (your existing code)
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

# Filter cells (your filters)
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



# ============================================
# Proper export using anndataR
# ============================================

library(Seurat)
library(anndataR)

# IMPORTANT: join layers in Seurat v5, combines all the split layers into single matrices
seu.int <- JoinLayers(seu.int)

# export to h5ad with a proper filename (not tempfile)
write_h5ad(seu.int, "seu_int.h5ad")

print("Export complete! File saved as: seu_int.h5ad")


# ---------------------------------------------------------------------------
# PRE-REQUISITE:
# You need 'SeuratDisk' installed.
remotes::install_github("mojaveazure/seurat-disk")
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# PRE-REQUISITE:
# You need 'SeuratDisk' installed.
# remotes::install_github("mojaveazure/seurat-disk")
# ---------------------------------------------------------------------------

library(Seurat)
library(SeuratDisk)

# ---------------------------------------------------------------------------
# FIX FOR SEURAT V5 COMPATIBILITY
# ---------------------------------------------------------------------------

# 1. Consolidate the split layers into one matrix
print("Joining layers...")
seu.int <- JoinLayers(seu.int)

# 2. Convert the Assay to v3 format
# SeuratDisk cannot yet handle the new Seurat v5 "Assay5" class. 
# Must convert it to the standard "Assay" class.
print("Converting Assay5 to v3 Assay for SeuratDisk compatibility...")

# Create a temporary v3 assay
seu.int[["RNA3"]] <- as(object = seu.int[["RNA"]], Class = "Assay")

# Set it as default and remove the old v5 assay
DefaultAssay(seu.int) <- "RNA3"
seu.int[["RNA"]] <- NULL
seu.int <- RenameAssays(seu.int, RNA3 = 'RNA')

# Check if data is present
if (is.null(seu.int@assays$RNA@data) || length(seu.int@assays$RNA@data) == 0) {
  print("Data slot appears empty. Running NormalizeData to populate it...")
  seu.int <- NormalizeData(seu.int)
}

# ---------------------------------------------------------------------------
# SAVE AND CONVERT
# ---------------------------------------------------------------------------

print("Saving h5Seurat...")
SaveH5Seurat(seu.int, filename = "gut_atlas_integrated.h5Seurat", overwrite = TRUE)

print("Converting to h5ad...")
Convert("gut_atlas_integrated.h5Seurat", dest = "h5ad", overwrite = TRUE)

print("Conversion complete: gut_atlas_integrated.h5ad ready for Myriad upload.")









# ============================================
# Manual export (UPDATED for Seurat v5)
# ============================================

library(Seurat)
library(Matrix)

# Create output directory
dir.create("seurat_export", showWarnings = FALSE)

print("Exporting Seurat object to files for Python...")

# 1. Export normalized data matrix (updated syntax)
print("Exporting normalized expression data...")
writeMM(
  obj = LayerData(seu.int, layer = "data", assay = "RNA"),
  file = "seurat_export/normalized.mtx"
)

# 2. Export raw counts (updated syntax)
print("Exporting raw counts...")
writeMM(
  obj = LayerData(seu.int, layer = "counts", assay = "RNA"),
  file = "seurat_export/counts.mtx"
)

# 3. Export gene names
print("Exporting gene names...")
write.table(
  rownames(seu.int),
  file = "seurat_export/genes.csv",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# 4. Export cell barcodes
print("Exporting cell barcodes...")
write.table(
  colnames(seu.int),
  file = "seurat_export/barcodes.csv",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# 5. Export metadata
print("Exporting metadata...")
write.csv(
  seu.int@meta.data,
  file = "seurat_export/metadata.csv",
  row.names = TRUE
)

# 6. Export UMAP coordinates
print("Exporting UMAP coordinates...")
write.csv(
  Embeddings(seu.int, reduction = "umap"),
  file = "seurat_export/umap_coords.csv",
  row.names = TRUE
)

# 7. Export PCA coordinates (optional but useful)
print("Exporting PCA coordinates...")
write.csv(
  Embeddings(seu.int, reduction = "pca"),
  file = "seurat_export/pca_coords.csv",
  row.names = TRUE
)

# 8. Export Harmony coordinates (optional but useful)
print("Exporting Harmony coordinates...")
write.csv(
  Embeddings(seu.int, reduction = "harmony"),
  file = "seurat_export/harmony_coords.csv",
  row.names = TRUE
)

# Summary
cat("\n")
cat(strrep("=", 60), "\n")
cat("Export complete!\n")
cat(strrep("=", 60), "\n")
cat("Number of cells:", ncol(seu.int), "\n")
cat("Number of genes:", nrow(seu.int), "\n")
cat("Tissue groups:", paste(unique(seu.int$tissue_group), collapse=", "), "\n")
cat("Patients:", paste(unique(seu.int$patient), collapse=", "), "\n")
cat("\n")
cat("Files exported to seurat_export/:\n")
cat("  - normalized.mtx (log-normalized expression)\n")
cat("  - counts.mtx (raw counts)\n")
cat("  - genes.csv (gene names)\n")
cat("  - barcodes.csv (cell barcodes)\n")
cat("  - metadata.csv (all cell metadata)\n")
cat("  - umap_coords.csv (your exact UMAP coordinates)\n")
cat("  - pca_coords.csv (PCA coordinates)\n")
cat("  - harmony_coords.csv (Harmony coordinates)\n")
cat("\n")
cat("Ready for Python CellTypist annotation!\n")
cat(strrep("=", 60), "\n")

# Optional: Quick verification
cat("\nVerification:\n")
cat("  Normalized data dimensions:", dim(LayerData(seu.int, layer = "data")), "\n")
cat("  Counts data dimensions:", dim(LayerData(seu.int, layer = "counts")), "\n")
cat("  UMAP dimensions:", dim(Embeddings(seu.int, reduction = "umap")), "\n")