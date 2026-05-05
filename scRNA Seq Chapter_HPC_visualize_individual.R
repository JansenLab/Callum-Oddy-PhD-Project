#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

# 1. Load Data
seu <- readRDS("celltypist_input/seu_int.rds")
preds <- read.csv("celltypist_input/hybrid_labels.csv", row.names = 1)
seu <- AddMetaData(seu, metadata = preds)

# 2. List of things we want to plot individually
plot_columns <- c("tissue_group", "seurat_clusters", "stomach_majority",
"intestine_majority", "hybrid_label")

# 3. Plotting Loop
pdf("celltypist_input/plots/individual_umaps_full.pdf", width = 12, height = 10)

for (col in plot_columns) {
  cat("Plotting:", col, "\n")

  # Dynamic title cleaning
  clean_title <- gsub("_", " ", col) %>% toupper()

  p <- DimPlot(seu, reduction = "umap", group.by = col, label = TRUE,
               label.size = 4, repel = TRUE, pt.size = 0.6) +
       theme_minimal() +
       labs(title = paste("UMAP Colored by:", clean_title)) +
       theme(legend.text = element_text(size = 8),
             legend.position = "right")

  print(p)
}

dev.off()

