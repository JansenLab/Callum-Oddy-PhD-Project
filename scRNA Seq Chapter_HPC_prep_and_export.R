#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

cat("--- R PHASE START ---\n")
seu <- readRDS("celltypist_input/seu_lineage_labeled.rds")

# 1. Update Annotations
neck_map <- c(
  "18" = "neck_mucous_proliferating",
  "14" = "neck_mucous_IM", "4" = "neck_mucous_IM", "19" = "neck_mucous_IM",
  "2"  = "neck_mucous_transitional",
  "0"  = "neck_mucous_gastric", "1"  = "neck_mucous_gastric", "6"  =
"neck_mucous_gastric",
  "13" = "neck_mucous_gastric", "20" = "neck_mucous_gastric", "10" =
"neck_mucous_gastric"
)

seu$final_annotation <- as.character(seu$paper_lineage)
cluster_ids <- as.character(seu$seurat_clusters)
for (num in names(neck_map)) {
  seu$final_annotation[cluster_ids == num] <- neck_map[num]
}

# 2. Setup Directories
bridge_dir <- "celltypist_input/bridge_files"
dir.create(bridge_dir, showWarnings = FALSE, recursive = TRUE)

# 3. Explicitly extract and name coordinates
cat("Extracting UMAP and PCA...\n")
umap_coords <- as.data.frame(Embeddings(seu, "umap"))
colnames(umap_coords) <- c("UMAP_1", "UMAP_2") # Force naming

pca_coords <- as.data.frame(Embeddings(seu, "pca"))
# Take the first 50 PCs
pca_coords <- pca_coords[, 1:min(50, ncol(pca_coords))]

# 4. Export
cat("Exporting Matrix (using layer='counts')...\n")
# Using layer instead of slot to satisfy Seurat V5
counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
writeMM(counts, file.path(bridge_dir, "matrix.mtx"))

metadata <- seu@meta.data
metadata$barcode <- rownames(metadata)
metadata <- cbind(metadata, umap_coords, pca_coords)
write.csv(metadata, file.path(bridge_dir, "metadata.csv"), row.names = FALSE)

write.table(rownames(seu), file.path(bridge_dir, "genes.tsv"), row.names = FALSE,
col.names = FALSE, quote = FALSE)
cat("--- R PHASE COMPLETE ---\n")

