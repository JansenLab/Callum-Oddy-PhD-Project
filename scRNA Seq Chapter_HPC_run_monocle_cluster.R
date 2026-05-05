library(monocle3)
library(Matrix)

path <- "celltypist_input/final_outputs/"
print("Loading data...")

# 1. Load the matrix and IMMEDIATELY transpose it
# This flips it from (Cells x Genes) to (Genes x Cells)
expression_matrix <- t(readMM(paste0(path, "counts.mtx")))

# 2. Load the metadata and gene names
cell_metadata <- read.csv(paste0(path, "metadata.csv"), row.names = 1)
gene_annotation <- read.csv(paste0(path, "gene_names.csv"))
umap_coords <- read.csv(paste0(path, "umap_coords.csv"), row.names = 1)

# 3. Now the dimensions will match perfectly
# Rows = Genes (38606), Columns = Cells (80133)
rownames(expression_matrix) <- gene_annotation$gene_short_name
colnames(expression_matrix) <- rownames(cell_metadata)
rownames(gene_annotation) <- gene_annotation$gene_short_name

print(paste("Matrix dimensions:", nrow(expression_matrix), "genes by",
ncol(expression_matrix), "cells"))

# 4. Continue with CDS creation as before...
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

# 5. Pre-process
# Even though we are using your UMAP coordinates, Monocle needs to
# run PCA internally to understand the global structure of the data.
print("Running preprocessing (PCA)...")
cds <- preprocess_cds(cds, num_dim = 50)

# 6. Inject your existing UMAP coordinates
# We convert the umap_coords dataframe to a matrix and slot it in.
# This ensures the trajectory lines match your known plot perfectly.
print("Injecting existing UMAP coordinates...")
reducedDims(cds)[["UMAP"]] <- as.matrix(umap_coords[colnames(cds), ])

# 7. Cluster the cells
# Monocle's graph-learning algorithm requires cells to be partitioned into clusters.
print("Clustering cells...")
cds <- cluster_cells(cds)

# 8. Learn the Trajectory Graph
# 'use_partition = FALSE' is critical here. It tells Monocle to try and
# connect all cells into one tree, even if they are in different UMAP clusters.
# This is usually best for seeing the transition into Intestinal Metaplasia.
print("Learning trajectory graph (this may take a few minutes)...")
cds <- learn_graph(cds, use_partition = FALSE)

# 9. Save the object for your laptop
# This is the file you will download to do the interactive plotting.
output_file <- paste0(path, "monocle3_umap_trajectory.rds")
print(paste("Saving final object to:", output_file))
saveRDS(cds, output_file)

print("Process complete! You can now transfer the .rds file to your local machine.")




