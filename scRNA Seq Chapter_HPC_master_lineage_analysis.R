#!/usr/bin/Rscript --no-save
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(viridis)
})

# ============================================================================
# 1. SETUP & DATA LOADING
# ============================================================================
cat("Loading labeled Seurat object...\n")
seu <- readRDS("celltypist_input/seu_lineage_labeled.rds")

# Create Output Directories
base_dir <- "celltypist_input/lineage_results/deep_dive"
csv_dir  <- file.path(base_dir, "csv_tables")
plot_dir <- file.path(base_dir, "plots")

dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Set max cells to prevent HPC memory crashes
max_cells <- 400

# ============================================================================
# 2. LEVEL 1: GLOBAL LINEAGE MARKERS
# ============================================================================
cat("\n--- LEVEL 1: Finding Global Lineage Markers ---\n")
Idents(seu) <- "paper_lineage"

lineage_markers <- FindAllMarkers(
  seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
  max.cells.per.ident = max_cells
)
write.csv(lineage_markers, file.path(csv_dir, "Level1_Global_Lineage_Markers.csv"), row.name$

# ============================================================================
# 3. LEVEL 2: GLOBAL CLUSTER MARKERS
# ============================================================================
cat("\n--- LEVEL 2: Finding Global Cluster Markers ---\n")
Idents(seu) <- "seurat_clusters"

cluster_markers <- FindAllMarkers(
  seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
  max.cells.per.ident = max_cells
)
write.csv(cluster_markers, file.path(csv_dir, "Level2_Global_Cluster_Markers.csv"), row.name$

# ============================================================================
# 4. LEVEL 3: INTRA-LINEAGE HETEROGENEITY
# ============================================================================
cat("\n--- LEVEL 3: Finding Intra-Lineage Sub-cluster Markers ---\n")
lineages <- sort(unique(as.character(seu$paper_lineage)))
intra_results <- list()

for (lin in lineages) {
  seu_sub <- subset(seu, subset = paper_lineage == lin)
  sub_clusters <- unique(as.character(seu_sub$seurat_clusters))

  if (length(sub_clusters) > 1) {
    cat(paste0("  -> Zooming in on ", lin, " (Composed of clusters: ", paste(sub_clusters, c$
    Idents(seu_sub) <- "seurat_clusters"

    intra_diffs <- tryCatch({
      FindAllMarkers(
        seu_sub, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.15,
        max.cells.per.ident = max_cells
      )
    }, error = function(e) return(NULL))

    if (!is.null(intra_diffs) && nrow(intra_diffs) > 0) {
      intra_diffs$lineage_context <- lin
      intra_results[[lin]] <- intra_diffs
      write.csv(intra_diffs, file.path(csv_dir, paste0("Level3_Intra_", gsub(" ", "_", lin),$
    }
  }
}


# ============================================================================
# 5. GENERATE MASTER PDF REPORT
# ============================================================================
cat("\n--- Generating Master PDF Report ---\n")
pdf(file.path(plot_dir, "Master_Cluster_and_Lineage_DeepDive.pdf"), width = 14, height = 10)

# --- Page 1: Lineage Defining Genes ---
cat("  Plotting Level 1: Lineages\n")
top_lin <- lineage_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
Idents(seu) <- "paper_lineage"
p1 <- DotPlot(seu, features = unique(top_lin$gene)) + coord_flip() + theme_minimal() +
  scale_color_viridis_c(option = "magma") +
  labs(title = "Level 1: What defines each Lineage?", x = "Top Marker Genes", y = "Lineage")$
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p1)

# --- Page 2: Cluster Defining Genes ---
cat("  Plotting Level 2: Clusters\n")
# Only taking top 2 genes per cluster so the plot doesn't become impossibly large
top_clust <- cluster_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
Idents(seu) <- "seurat_clusters"
p2 <- DotPlot(seu, features = unique(top_clust$gene)) + coord_flip() + theme_minimal() +
  scale_color_viridis_c(option = "viridis") +
  labs(title = "Level 2: What defines each Numerical Cluster?", x = "Top Marker Genes", y = $
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), axis.text.y = element_tex$
print(p2)
# --- Pages 3+: Intra-Lineage Deep Dives ---
cat("  Plotting Level 3: Intra-Lineage Zoom-ins\n")
for (lin in names(intra_results)) {
  seu_sub <- subset(seu, subset = paper_lineage == lin)
  Idents(seu_sub) <- "seurat_clusters"
  sub_markers <- intra_results[[lin]]

  # A. The UMAP of just this lineage
  p_umap <- DimPlot(seu_sub, group.by = "seurat_clusters", label = TRUE, pt.size = 1) +
    ggtitle(paste("Zoomed UMAP:", lin)) + theme_void()

  # B. The DotPlot of what separates the clusters inside this lineage
  top_sub <- sub_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  p_dot <- DotPlot(seu_sub, features = unique(top_sub$gene)) +
    theme_bw() + scale_color_viridis_c(option = "plasma") +
    labs(title = "What separates these sub-clusters?", x = "Genes", y = "Cluster #") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # C. Combine and print
  combined_plot <- (p_umap | p_dot) + plot_layout(widths = c(1, 2))
  print(combined_plot)
}

dev.off()
cat("\nAnalysis Finished! Review CSVs in 'csv_tables/' and the PDF in 'plots/'.\n")

