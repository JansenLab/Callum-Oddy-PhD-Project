#!/usr/bin/env Rscript

# --- 1. FORCE CORRECT PYTHON PATH ---
# Must be set before loading anndata/reticulate
Sys.setenv(RETICULATE_PYTHON =
"/home/regmddy/ACFS/Programmes/miniconda/envs/sc_env/bin/python")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(ggplot2)
  library(dplyr)
  library(anndata)
  library(hdf5r)
  library(Matrix)
  library(reticulate)
})

# Verify python binding
use_python(Sys.getenv("RETICULATE_PYTHON"), required = TRUE)
set.seed(1234)

cat("==============================================\n")
cat("Starting Seurat v5 Pipeline\n")
cat("Python Path:", py_config()$python, "\n")
cat("==============================================\n")

# --- 2. LOAD DATA ---
cat("Loading 10X h5 files...\n")
data_dirs <- list(
  "CO_06_GAS" = "./Filtered_Features/CO_06_GAS_GEX/filtered_feature_bc_matrix.h5",
  "CO_06_IM1" = "./Filtered_Features/CO_06_IM1_GEX/filtered_feature_bc_matrix.h5",
  "CO_06_IM2" = "./Filtered_Features/CO_06_IM2_GEX/filtered_feature_bc_matrix.h5",
  "CO_06_DUO" = "./Filtered_Features/CO_06_DUO_GEX/filtered_feature_bc_matrix.h5",
  "CO_07_GAS" = "./Filtered_Features/CO_07_GAS_GEX/filtered_feature_bc_matrix.h5",
  "CO_07_IM1" = "./Filtered_Features/CO_07_IM1_GEX/filtered_feature_bc_matrix.h5",
  "CO_07_IM2" = "./Filtered_Features/CO_07_IM2_GEX/filtered_feature_bc_matrix.h5",
  "CO_07_DUO" = "./Filtered_Features/CO_07_DUO_GEX/filtered_feature_bc_matrix.h5"
)

# Loop to create objects and perform basic QC
obj_list <- lapply(names(data_dirs), function(name) {
  cat("Processing:", name, "\n")
  mtx <- Read10X_h5(data_dirs[[name]])
  obj <- CreateSeuratObject(counts = mtx, project = name)
  obj[["percent.mito"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  return(subset(obj, nFeature_RNA > 400 & percent.mito < 10))
})

# --- 3. MERGE & JOIN LAYERS ---
cat("Merging and joining layers...\n")
seu.int <- merge(obj_list[[1]], y = obj_list[-1])
seu.int <- JoinLayers(seu.int)
rm(obj_list); gc()

# --- 4. WORKFLOW (CORRECTED ORDER - MATCHING SCRIPT 1) ---
cat("Normalization & Finding Variable Features...\n")
seu.int <- NormalizeData(seu.int)
seu.int <- FindVariableFeatures(seu.int)

# Scale WITHOUT regression first (for initial PCA)
cat("Initial scaling (no regression)...\n")
seu.int <- ScaleData(seu.int, features = VariableFeatures(seu.int))

cat("Running PCA...\n")
seu.int <- RunPCA(seu.int, features = VariableFeatures(seu.int), verbose = FALSE)

# Add metadata for Harmony
cat("Adding metadata...\n")
seu.int$patient <- sapply(strsplit(seu.int$orig.ident, "_"), function(x) paste(x[1],
x[2], sep = "_"))
seu.int$tissue_type <- sapply(strsplit(seu.int$orig.ident, "_"), function(x) x[3])
seu.int$tissue_group <- ifelse(seu.int$tissue_type %in% c("IM1", "IM2"), "IM", seu.int$tissu$

# Run Harmony BEFORE cell cycle regression
cat("Running Harmony batch correction...\n")
seu.int <- RunHarmony(seu.int, group.by.vars = "patient", verbose = FALSE)

# Cell cycle scoring AFTER Harmony
cat("Cell cycle scoring...\n")
seu.int <- CellCycleScoring(
  object = seu.int,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident = TRUE
)

# Regress out cell cycle AFTER Harmony
cat("Scaling with cell cycle regression...\n")
seu.int <- ScaleData(
  seu.int,
  vars.to.regress = c("S.Score", "G2M.Score"),
  features = VariableFeatures(seu.int)
)

# UMAP and clustering on Harmony embeddings
cat("Running UMAP and clustering...\n")
seu.int <- RunUMAP(seu.int, reduction = "harmony", dims = 1:30, verbose = FALSE)
seu.int <- FindNeighbors(seu.int, reduction = "harmony", dims = 1:30, verbose = FALSE)
seu.int <- FindClusters(seu.int, resolution = 2, verbose = FALSE)

# --- 5. VISUALIZATION (MUTED COLORS) ---
dir.create("celltypist_input", showWarnings = FALSE)

# Custom Muted Colors
# GAS = Muted Red, IM = Muted Purple, DUO = Muted Blue
muted_cols <- c("GAS" = "#CD5C5C", "IM" = "#9370DB", "DUO" = "#4682B4")

cat("Saving plots...\n")
pdf("celltypist_input/umap_visualizations.pdf", width = 12, height = 10)
# Tissue Plot with Muted Colors
p1 <- DimPlot(seu.int, reduction = "umap", group.by = "tissue_group", pt.size = 0.5) +
      scale_color_manual(values = muted_cols) +
      ggtitle("Tissue Groups (Muted Palette)")
print(p1)

# Cluster Plot
p2 <- DimPlot(seu.int, reduction = "umap", label = TRUE) + ggtitle("Clusters (Res
2.0)")
print(p2)

# Cell Cycle Phase Plot
p3 <- DimPlot(seu.int, reduction = "umap", group.by = "Phase", pt.size = 0.5) + ggtitle("Cel$
print(p3)

# Patient Plot
p4 <- DimPlot(seu.int, reduction = "umap", group.by = "patient", pt.size = 0.5) + ggtitle("P$
print(p4)

dev.off()

# --- 6. EXPORT TO ANNDATA ---
cat("Exporting to h5ad...\n")
# Ensure we use the correct sparse matrix format
data_mat <- GetAssayData(seu.int, assay = "RNA", layer = "data")

adata <- anndata::AnnData(
  X = Matrix::t(data_mat),
  obs = seu.int@meta.data,
  obsm = list(
    X_umap = Embeddings(seu.int, "umap"),
    X_pca = Embeddings(seu.int, "pca"),
    X_harmony = Embeddings(seu.int, "harmony")
  )
)

anndata::write_h5ad(adata, "celltypist_input/seu_int.h5ad")
saveRDS(seu.int, "celltypist_input/seu_int.rds")

cat("\n==============================================\n")
cat("Processing complete!\n")
cat("Output files in celltypist_input/:\n")
cat("  - seu_int.h5ad (for CellTypist)\n")
cat("  - seu_int.rds (Seurat backup)\n")
cat("  - umap_visualizations.pdf (all plots)\n")
cat("==============================================\n")
cat("Job Finished Successfully!\n")


