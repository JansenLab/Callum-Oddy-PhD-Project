if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'ggrastr'))

remotes::install_github("bnprks/BPCells/r")

devtools::install_github('cole-trapnell-lab/monocle3')

library(monocle3)

library(monocle3)

setwd("setdirectory")


# Load the processed object from Myriad
cds <- readRDS("monocle3_umap_trajectory.rds")

# Check it loaded correctly
print(cds)

# Plot to see UMAP with cell types
plot_cells(cds, 
           color_cells_by = "final_annotation",
           label_groups_by_cluster = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE) +
  scale_color_manual(values = annot_cols)

# Order cells in pseudotime - first visualise to pick your root
plot_cells(cds,
           label_cell_groups = FALSE,
           label_leaves = TRUE,
           label_branch_points = TRUE,
           graph_label_size = 1.5)

# Then order cells - this will open an interactive window to click root node
cds <- order_cells(cds)

# Plot pseudotime
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           graph_label_size = 1.5)

# Save with pseudotime added
saveRDS(cds, "monocle3_pseudotime_final_cytotrace2node.rds")
save_monocle_objects(cds=cds, directory_path = "setdirectory")




