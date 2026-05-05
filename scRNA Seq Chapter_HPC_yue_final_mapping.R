#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

cat("Loading Seurat object...\n")
seu <- readRDS("celltypist_input/seu_int.rds")
Idents(seu) <- "seurat_clusters"

# ============================================================================
# METADATA CONFIGURATION
# ============================================================================
# Based on your previous output, the column is likely 'tissue_group'
organoid_col <- "tissue_type"
patient_col  <- "patient"

# ============================================================================
# 1. OFFICIAL MARKER LIST
# ============================================================================
marker_list <- list(
  Gastric_Lineage = c("CLDN18", "ANXA10"),
  Intestinal_Lineage = c("CDX1", "CDX2", "CLDN3"),
  Stem = c("TNFRSF19", "LGR5", "EPHB2", "OLFM4", "SMOC2", "ASCL2", "MEX3A"),
  IsthSC = c("STMN1"),
  Proliferating = c("HELLS", "PCNA", "TOP2A", "MKI67", "BIRC5"),
  Chief = c("PGC", "LIPF", "PGA3"),
  Neck_Mucous = c("MUC6"),
  Gastric_Mucous = c("TFF2", "TFF1", "MUC5AC"),
  Parietal = c("ATP4A"),
  Enterocyte = c("FABP1", "KRT20", "ANPEP"),
  Goblet_Precursor = c("HES6"),
  Goblet = c("SPINK4", "ATOH1"),
  Paneth = c("DEFA6", "PLA2G2A"),
  EE_Progenitor = c("PROX1", "BMI1"),
  EE = c("CHGA", "NEUROG3", "SYP"),
  Tuft = c("POU2F3", "AVIL")
)

all_genes <- unique(unlist(marker_list))
valid_genes <- intersect(all_genes, rownames(seu))
clusters <- as.character(sort(unique(Idents(seu))))

# ============================================================================
# 2. THE 20% FILTER & Z-SCORES (YUE METHOD)
# ============================================================================
cat("Filtering genes by 20% expression threshold per cluster...\n")
expr_matrix <- GetAssayData(seu, layer = "data")
avg_table <- matrix(0, nrow = length(clusters), ncol = length(valid_genes),
                    dimnames = list(clusters, valid_genes))

for (cl in clusters) {
  cells_in_cl <- WhichCells(seu, idents = cl)
  for (gene in valid_genes) {
    pct <- sum(expr_matrix[gene, cells_in_cl] > 0) / length(cells_in_cl)
    if (pct >= 0.20) {
      avg_table[cl, gene] <- mean(expr_matrix[gene, cells_in_cl])
    }
  }
}

z_table <- as.data.frame(scale(avg_table))
z_table[is.na(z_table)] <- 0

# ============================================================================
# 3. WINNER-TAKES-ALL SCORING
# ============================================================================
cat("Calculating relative lineage scores...\n")
lineage_scores <- data.frame(row.names = clusters)
for (lineage in names(marker_list)) {
  genes <- intersect(marker_list[[lineage]], colnames(z_table))
  if (length(genes) > 1) {
    lineage_scores[[lineage]] <- rowMeans(z_table[, genes])
  } else if (length(genes) == 1) {
    lineage_scores[[lineage]] <- z_table[[genes]]
  } else {
    lineage_scores[[lineage]] <- -10
  }
}

cat("Assigning final identities...\n")
final_labels <- sapply(clusters, function(cl) {
  scores <- lineage_scores[cl, ]
  if (scores$Enterocyte > 0 & scores$Gastric_Mucous > 0) {
    k20 <- if ("KRT20" %in% colnames(z_table)) z_table[cl, "KRT20"] else 0
    if (k20 > 1.0) return("Late PitEntero")
    if (k20 > 0.2) return("Mid PitEntero")
    return("Early PitEntero")
  }
  winner <- names(which.max(scores))
  if (max(scores) < -0.5) return("Quiescent/Low Quality")
  return(winner)
})

mapping_table <- data.frame(Seurat_Cluster = names(final_labels), Assigned_Identity = unname$
seu$paper_lineage <- unname(setNames(final_labels, clusters)[as.character(seu$seurat_cluster$
Idents(seu) <- "paper_lineage"

# ============================================================================
# 4. PREPARE HEATMAP DATA
# ============================================================================
cat("Finding Top 30 Variable Genes...\n")
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
top30_hvg <- head(VariableFeatures(seu), 30)
seu <- ScaleData(seu, features = unique(c(top30_hvg, valid_genes)), verbose = FALSE)
seu_sub <- subset(seu, downsample = 100)

# ============================================================================
# 5. GENERATE COMPREHENSIVE PDF
# ============================================================================
out_dir <- "celltypist_input/lineage_results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Generating expanded PDF plots...\n")
pdf(file.path(out_dir, "Full_Biological_Validation_Report.pdf"), width = 16, height = 12)

# PAGE 1: Master UMAP
print(DimPlot(seu, group.by = "paper_lineage", label = TRUE, repel = TRUE) +
      ggtitle("Final Lineage Assignments (Yue et al. Logic)") + theme_minimal())

# PAGE 2: Organoid Contribution to Lineage (With Custom Muted Colors)
cat("Plotting tissue composition...\n")

# Logic to group IM1/IM2 into just "IM" for the color mapping
plot_meta <- seu@meta.data
plot_meta$tissue_color_group <- as.character(plot_meta[[organoid_col]])
plot_meta$tissue_color_group[grepl("IM", plot_meta$tissue_color_group)] <- "IM"
plot_meta$tissue_color_group <- factor(plot_meta$tissue_color_group, levels = c("GAS", "IM",$
custom_colors <- c("GAS" = "#C96156", "IM" = "#8C77A6", "DUO" = "#608AB5")

print(ggplot(plot_meta, aes(x = paper_lineage, fill = tissue_color_group)) +
  geom_bar(position = "fill", color = "black", size = 0.2) +
  scale_fill_manual(values = custom_colors, name = "Tissue Source") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Contribution of Tissue Type to Lineages", y = "Percent", x = "Lineage"))

# NEW PAGE 3: Cluster Contribution with Numerical Labels
cat("Plotting cluster breakdown with labels...\n")

# Create summary data for labeling
label_data <- plot_meta %>%
  group_by(paper_lineage, seurat_clusters) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(paper_lineage) %>%
  mutate(pct = count / sum(count)) %>%
  filter(pct > 0.03)
  print(ggplot(plot_meta, aes(x = paper_lineage, fill = seurat_clusters)) +
  geom_bar(position = "fill", color = "black", size = 0.2) +
  geom_text(data = label_data, aes(y = pct, label = seurat_clusters),
            position = position_stack(vjust = 0.5), size = 3, fontface = "bold") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Which Numerical Clusters form each Lineage?",
       subtitle = "Labels shown for clusters > 3% contribution",
       y = "Percent of Cluster Contribution", x = "Lineage", fill = "Cluster #"))

# PAGE 4: Patient Contribution to Clusters
cat("Plotting patient breakdown of clusters...\n")
print(ggplot(seu@meta.data, aes(x = seurat_clusters, fill = !!sym(patient_col))) +
  geom_bar(position = "fill", color = "black", size = 0.2) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Is there a Patient Batch Effect per Cluster?",
       y = "Percent of Patient Contribution", x = "Seurat Cluster Number", fill = "Patient"))

# PAGE 5: Top 30 HVG Heatmap
print(DoHeatmap(seu_sub, features = top30_hvg, group.by = "paper_lineage") +
      theme(axis.text.y = element_text(size = 9, face = "bold.italic")) +
      scale_fill_gradientn(colors = c("navy", "white", "firebrick3")) +
      ggtitle("Top 30 Most Variable Genes"))

# PAGE 6: Supervised Marker Heatmap
print(DoHeatmap(seu_sub, features = valid_genes, group.by = "paper_lineage") +
      theme(axis.text.y = element_text(size = 7, face = "italic")) +
      scale_fill_gradientn(colors = c("navy", "white", "firebrick3")) +
      ggtitle("Canonical Lineage Markers"))
# PAGE 7+: Individual Highlighted UMAPs
lineages <- sort(unique(seu$paper_lineage))
plot_list <- list()
for (lin in lineages) {
  plot_list[[lin]] <- DimPlot(seu, cells.highlight = WhichCells(seu, idents = lin),
               cols.highlight = "firebrick3", cols = "lightgrey") +
       ggtitle(lin) + theme_void() + theme(legend.position = "none", plot.title = element_te$
}
for (i in seq(1, length(plot_list), by = 4)) {
  indices <- i:min(i+3, length(plot_list))
  print(wrap_plots(plot_list[indices], ncol = 2))
}

dev.off()

# EXPORTS
write.csv(mapping_table, file.path(out_dir, "Cluster_to_Identity_Mapping.csv"), row.names = $
saveRDS(seu, "celltypist_input/seu_lineage_labeled.rds")
cat("Finished! Check 'celltypist_input/lineage_results/Full_Biological_Validation_Report.pdf$










