library(Seurat)
library(anndataR)

# 1. Load your object
seu <- readRDS("celltypist_input/seu_cytotrace_complete.rds")

# 2. Prepare for Seurat v5 (Merges layers if they were split by patient)
seu <- JoinLayers(seu)

# 3. Convert to AnnData with precise mapping
# We map the 'harmony' reduction to 'X_harmony' for Python/Scanpy compatibility
adata <- as_AnnData(
  seu,
  assay_name = "RNA",
  x_mapping = "data", # Use normalized data
  obsm_mapping = list(
    X_pca = "pca",
    X_harmony = "harmony",
    X_umap = "umap"
  ),
  obs_mapping = TRUE # Keeps all metadata like CytoTRACE2_Score and tissue_group
)

# 4. Save the file
adata$write_h5ad("seu_integrated_final.h5ad")

message("--- SUCCESS: seu_integrated_final.h5ad is ready for PHATE ---")
