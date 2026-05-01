# ==============================================================================
# CyteType Results Visualization: tables, heatmaps, and plots
# ==============================================================================

# Load libraries
library(tidyverse)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(gridExtra)
library(gt)
library(gtExtras)
library(viridis)
library(scales)

# Set working directory
setwd("setdirectory")

# Create output directories
dir.create("figures/cytetype_analysis", recursive = TRUE, showWarnings = FALSE)
dir.create("tables/cytetype_analysis", recursive = TRUE, showWarnings = FALSE)

cat("================================================================================\n")
cat("CyteType Results Visualization\n")
cat("================================================================================\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
cat("Loading data...\n")

celltype_organoid_pct <- read.csv("celltype_by_organoid_percentage.csv", row.names = 1)
celltype_organoid_counts <- read.csv("celltype_by_organoid_counts.csv", row.names = 1)
celltype_patient_pct <- read.csv("celltype_by_patient_percentage.csv", row.names = 1)
celltype_patient_counts <- read.csv("celltype_by_patient_counts.csv", row.names = 1)
overall_summary <- read.csv("celltype_overall_summary.csv")
cluster_mapping <- read.csv("cluster_to_celltype_mapping.csv")
celltype_sample <- read.csv("celltype_by_sample_counts.csv", row.names = 1)
dataset_stats <- read.csv("dataset_summary_statistics.csv")

cat("  ✓ All files loaded successfully\n\n")

# ==============================================================================
# 1B. COMBINE IM1 AND IM2, SET ORDER
# ==============================================================================
cat("Combining IM1 and IM2 into IM category...\n")

# Function to combine IM1 and IM2 rows
combine_IM <- function(df) {
  if ("IM1" %in% rownames(df) && "IM2" %in% rownames(df)) {
    # Sum IM1 and IM2
    im_combined <- colSums(df[c("IM1", "IM2"), , drop = FALSE])
    # Remove IM1 and IM2
    df <- df[!rownames(df) %in% c("IM1", "IM2"), , drop = FALSE]
    # Add combined IM
    df <- rbind(df, IM = im_combined)
  }
  return(df)
}

# Apply to count data
celltype_organoid_counts <- combine_IM(celltype_organoid_counts)
celltype_patient_counts <- combine_IM(celltype_patient_counts)

# Recalculate percentages after combining
celltype_organoid_pct <- sweep(celltype_organoid_counts, 1, 
                               rowSums(celltype_organoid_counts), FUN = "/") * 100
celltype_patient_pct <- sweep(celltype_patient_counts, 1, 
                              rowSums(celltype_patient_counts), FUN = "/") * 100

# Set factor order: GAS -> IM -> DUO
organoid_order <- c("GAS", "IM", "DUO")

# Reorder rows
celltype_organoid_pct <- celltype_organoid_pct[organoid_order, ]
celltype_organoid_counts <- celltype_organoid_counts[organoid_order, ]

# Define organoid colors (updated for 3 categories)
organoid_colors <- c("GAS" = "#C25759", "IM" = "#8E7DBE", "DUO" = "#5B8FA8")

cat("  ✓ IM1 and IM2 combined into IM\n")
cat("  ✓ Order set to: GAS -> IM -> DUO\n\n")

# ==============================================================================
# 2. PUBLICATION-QUALITY TABLES
# ==============================================================================
cat("Creating publication-quality tables...\n")

# Table 1: Overall cell type summary with formatting
overall_table <- overall_summary %>%
  arrange(desc(Count)) %>%
  gt() %>%
  tab_header(
    title = "Cell Type Distribution",
    subtitle = "Overall dataset summary"
  ) %>%
  cols_label(
    Cell_Type = "Cell Type",
    Count = "Count",
    Percentage = "Percentage (%)"
  ) %>%
  fmt_number(
    columns = Count,
    decimals = 0,
    use_seps = TRUE
  ) %>%
  fmt_number(
    columns = Percentage,
    decimals = 2
  ) %>%
  data_color(
    columns = Percentage,
    colors = scales::col_numeric(
      palette = c("white", "steelblue"),
      domain = NULL
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  opt_row_striping()

gtsave(overall_table, "tables/cytetype_analysis/table1_overall_summary.html")
gtsave(overall_table, "tables/cytetype_analysis/table1_overall_summary.png")

# Table 2: Cell type by organoid type
organoid_table <- celltype_organoid_pct %>%
  rownames_to_column("Organoid_Type") %>%
  gt() %>%
  tab_header(
    title = "Cell Type Distribution by Organoid Type",
    subtitle = "Percentage of cells in each category"
  ) %>%
  fmt_number(
    columns = -Organoid_Type,
    decimals = 2
  ) %>%
  data_color(
    columns = -Organoid_Type,
    colors = scales::col_numeric(
      palette = viridis(100),
      domain = c(0, 100)
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  )

gtsave(organoid_table, "tables/cytetype_analysis/table2_organoid_distribution.html")
gtsave(organoid_table, "tables/cytetype_analysis/table2_organoid_distribution.png")

# Table 3: Cluster mapping
cluster_table <- cluster_mapping %>%
  arrange(Cluster) %>%
  gt() %>%
  tab_header(
    title = "Cluster to Cell Type Mapping",
    subtitle = "CyteType annotations with ontology terms"
  ) %>%
  cols_label(
    Cluster = "Cluster",
    Cell_Type = "Cell Type",
    Ontology_Term = "Cell Ontology ID",
    Cell_Count = "Cell Count"
  ) %>%
  fmt_number(
    columns = Cell_Count,
    decimals = 0,
    use_seps = TRUE
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  opt_row_striping()

gtsave(cluster_table, "tables/cytetype_analysis/table3_cluster_mapping.html")
gtsave(cluster_table, "tables/cytetype_analysis/table3_cluster_mapping.png")

cat("  ✓ Tables saved to tables/cytetype_analysis/\n\n")

# ==============================================================================
# 3. HEATMAPS
# ==============================================================================
cat("Creating heatmaps...\n")

# Heatmap 1: Cell type by organoid type (percentage)
pdf("figures/cytetype_analysis/heatmap1_organoid_percentage.pdf", width = 12, height = 8)
pheatmap(
  t(celltype_organoid_pct),
  color = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.1f",
  fontsize_number = 10,
  main = "Cell Type Distribution by Organoid Type (%)",
  angle_col = "0",
  cellwidth = 40,
  cellheight = 15,
  border_color = "grey60"
)
dev.off()

# Heatmap 2: Cell type by organoid type (log counts for better visualization)
log_counts <- log10(celltype_organoid_counts + 1)
pdf("figures/cytetype_analysis/heatmap2_organoid_log_counts.pdf", width = 12, height = 8)
pheatmap(
  t(log_counts),
  color = colorRampPalette(c("white", "orange", "red", "darkred"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "Cell Type Distribution by Organoid Type (log10 counts)",
  angle_col = "0",
  cellwidth = 40,
  cellheight = 15,
  border_color = "grey60"
)
dev.off()

# Heatmap 3: ComplexHeatmap version with better annotations
library(ComplexHeatmap)

# Prepare data
heatmap_data <- as.matrix(t(celltype_organoid_pct))

# Create color functions
col_fun <- colorRamp2(
  c(0, 25, 50, 75, 100),
  c("white", "lightblue", "steelblue", "darkblue", "black")
)

# Column annotations (organoid types)
organoid_colors <- c("GAS" = "#C25759", "IM" = "#8E7DBE", "DUO" = "#5B8FA8")
col_anno <- HeatmapAnnotation(
  Organoid = factor(colnames(heatmap_data), levels = organoid_order),
  col = list(Organoid = organoid_colors),
  show_legend = TRUE,
  annotation_name_side = "left"
)

# Create heatmap
pdf("figures/cytetype_analysis/heatmap3_complex_heatmap.pdf", width = 12, height = 10)
ht <- Heatmap(
  heatmap_data,
  name = "Percentage",
  col = col_fun,
  top_annotation = col_anno,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 12),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.1f", heatmap_data[i, j]), x, y, gp = gpar(fontsize = 8))
  },
  column_title = "Cell Type Distribution by Organoid Type (%)",
  row_title = "Cell Types",
  heatmap_legend_param = list(
    title = "Percentage",
    at = c(0, 25, 50, 75, 100),
    labels = c("0", "25", "50", "75", "100")
  )
)
draw(ht)
dev.off()

# Heatmap 4: Sample-level heatmap (more detailed)
# Only show top 15 most abundant cell types for clarity
top_celltypes <- overall_summary %>%
  arrange(desc(Count)) %>%
  head(15) %>%
  pull(Cell_Type)

sample_subset <- celltype_sample[, colnames(celltype_sample) %in% top_celltypes, drop = FALSE]

# Remove samples with zero counts
sample_subset <- sample_subset[rowSums(sample_subset) > 0, , drop = FALSE]

# Calculate percentages
sample_pct <- sweep(sample_subset, 1, rowSums(sample_subset), FUN = "/") * 100

# Only plot if we have valid data
if (nrow(sample_pct) > 0 && ncol(sample_pct) > 0 && !all(is.na(sample_pct))) {
  pdf("figures/cytetype_analysis/heatmap4_sample_level.pdf", width = 14, height = 10)
  pheatmap(
    t(sample_pct),
    color = viridis(100),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = FALSE,
    main = "Top 15 Cell Types by Sample",
    fontsize = 9,
    border_color = "grey60"
  )
  dev.off()
  cat("  ✓ Sample-level heatmap created\n")
} else {
  cat("  ⚠ Skipping sample-level heatmap (insufficient data)\n")
}

cat("  ✓ Heatmaps saved to figures/cytetype_analysis/\n\n")

# ==============================================================================
# 4. BAR PLOTS
# ==============================================================================
cat("Creating bar plots...\n")

# Plot 1: Overall cell type distribution
p1 <- ggplot(overall_summary %>% arrange(desc(Count)) %>% head(15), 
             aes(x = reorder(Cell_Type, Count), y = Count, fill = Cell_Type)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label = scales::comma(Count)), hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_fill_viridis_d(option = "turbo") +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Top 15 Cell Types by Abundance",
    x = "Cell Type",
    y = "Number of Cells"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

ggsave("figures/cytetype_analysis/barplot1_overall_abundance.pdf", 
       p1, width = 10, height = 8)

# Plot 2: Stacked bar chart by organoid type
celltype_long <- celltype_organoid_pct %>%
  rownames_to_column("Organoid_Type") %>%
  pivot_longer(-Organoid_Type, names_to = "Cell_Type", values_to = "Percentage") %>%
  filter(Percentage > 0) %>%
  mutate(Organoid_Type = factor(Organoid_Type, levels = organoid_order))

p2 <- ggplot(celltype_long, aes(x = Organoid_Type, y = Percentage, fill = Cell_Type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(length(unique(celltype_long$Cell_Type)))) +
  scale_x_discrete(limits = organoid_order) +
  labs(
    title = "Cell Type Composition by Organoid Type",
    x = "Organoid Type",
    y = "Percentage (%)",
    fill = "Cell Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right",
    legend.text = element_text(size = 8)
  )

ggsave("figures/cytetype_analysis/barplot2_stacked_organoid.pdf", 
       p2, width = 12, height = 8)

# Plot 3: Grouped bar chart showing comparison
celltype_comparison <- celltype_organoid_pct %>%
  rownames_to_column("Organoid_Type") %>%
  pivot_longer(-Organoid_Type, names_to = "Cell_Type", values_to = "Percentage") %>%
  group_by(Cell_Type) %>%
  filter(max(Percentage) > 5) %>%  # Only show cell types >5% somewhere
  ungroup() %>%
  mutate(Organoid_Type = factor(Organoid_Type, levels = organoid_order))

p3 <- ggplot(celltype_comparison, aes(x = reorder(Cell_Type, -Percentage), 
                                      y = Percentage, fill = Organoid_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = organoid_colors, limits = organoid_order) +
  coord_flip() +
  labs(
    title = "Cell Type Distribution Across Organoid Types",
    subtitle = "Cell types with >5% abundance in at least one organoid type",
    x = "Cell Type",
    y = "Percentage (%)",
    fill = "Organoid Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    legend.position = "bottom"
  )

ggsave("figures/cytetype_analysis/barplot3_grouped_comparison.pdf", 
       p3, width = 12, height = 10)

cat("  ✓ Bar plots saved to figures/cytetype_analysis/\n\n")

# ==============================================================================
# 5. DOT PLOTS / BUBBLE PLOTS
# ==============================================================================
cat("Creating dot plots...\n")

# Bubble plot showing abundance and prevalence
bubble_data <- celltype_organoid_counts %>%
  rownames_to_column("Organoid_Type") %>%
  pivot_longer(-Organoid_Type, names_to = "Cell_Type", values_to = "Count") %>%
  left_join(
    celltype_organoid_pct %>%
      rownames_to_column("Organoid_Type") %>%
      pivot_longer(-Organoid_Type, names_to = "Cell_Type", values_to = "Percentage"),
    by = c("Organoid_Type", "Cell_Type")
  ) %>%
  filter(Count > 0) %>%
  mutate(Organoid_Type = factor(Organoid_Type, levels = organoid_order))

p4 <- ggplot(bubble_data, aes(x = Organoid_Type, y = Cell_Type, 
                              size = Count, color = Percentage)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(2, 15), labels = scales::comma) +
  scale_color_viridis_c(option = "plasma") +
  scale_x_discrete(limits = organoid_order) +
  labs(
    title = "Cell Type Abundance Across Organoid Types",
    x = "Organoid Type",
    y = "Cell Type",
    size = "Cell Count",
    color = "Percentage (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 9),
    panel.grid.major = element_line(color = "grey90"),
    legend.position = "right"
  )

ggsave("figures/cytetype_analysis/dotplot1_bubble_abundance.pdf", 
       p4, width = 10, height = 12)

cat("  ✓ Dot plots saved to figures/cytetype_analysis/\n\n")

# ==============================================================================
# 6. PROPORTION PLOTS
# ==============================================================================
cat("Creating proportion plots...\n")

# Alluvial/Sankey-style plot showing cell type flow
# This requires counts data
flow_data <- celltype_organoid_counts %>%
  rownames_to_column("Organoid_Type") %>%
  pivot_longer(-Organoid_Type, names_to = "Cell_Type", values_to = "Count") %>%
  filter(Count > 100) %>%  
  mutate(Organoid_Type = factor(Organoid_Type, levels = organoid_order))

p5 <- ggplot(flow_data, aes(x = Organoid_Type, y = Count, fill = Cell_Type)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(flow_data$Cell_Type)))) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(limits = organoid_order) +
  labs(
    title = "Cell Type Proportions by Organoid Type",
    x = "Organoid Type",
    y = "Proportion",
    fill = "Cell Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right",
    legend.text = element_text(size = 8)
  )

ggsave("figures/cytetype_analysis/proportion1_filled_bars.pdf", 
       p5, width = 12, height = 8)

cat("  ✓ Proportion plots saved to figures/cytetype_analysis/\n\n")

# ==============================================================================
# 7. SUMMARY FIGURE (MULTI-PANEL)
# ==============================================================================
cat("Creating summary figure...\n")

# Create a multi-panel summary figure
library(patchwork)

# Reuse some plots but make them smaller
summary_fig <- (p1 + p3) / (p2 + p5) +
  plot_annotation(
    title = "CyteType Analysis Summary",
    subtitle = "Cell type distributions across gastric and duodenal organoids",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

ggsave("figures/cytetype_analysis/summary_figure_multipanel.pdf", 
       summary_fig, width = 20, height = 16)

cat("  ✓ Summary figure saved\n\n")

# ==============================================================================
# 8. EXPORT SESSION INFO
# ==============================================================================
sink("tables/cytetype_analysis/session_info.txt")
cat("R Session Information\n")
cat("====================\n\n")
sessionInfo()
sink()

# ==============================================================================
# COMPLETION MESSAGE
# ==============================================================================
cat("================================================================================\n")
cat("COMPLETE!\n")
cat("================================================================================\n\n")
cat("Generated files:\n")
cat("  Tables (HTML & PNG):\n")
cat("    - table1_overall_summary\n")
cat("    - table2_organoid_distribution\n")
cat("    - table3_cluster_mapping\n\n")
cat("  Heatmaps (PDF):\n")
cat("    - heatmap1_organoid_percentage\n")
cat("    - heatmap2_organoid_log_counts\n")
cat("    - heatmap3_complex_heatmap\n")
cat("    - heatmap4_sample_level\n\n")
cat("  Bar Plots (PDF):\n")
cat("    - barplot1_overall_abundance\n")
cat("    - barplot2_stacked_organoid\n")
cat("    - barplot3_grouped_comparison\n\n")
cat("  Dot Plots (PDF):\n")
cat("    - dotplot1_bubble_abundance\n\n")
cat("  Proportion Plots (PDF):\n")
cat("    - proportion1_filled_bars\n\n")
cat("  Summary:\n")
cat("    - summary_figure_multipanel (20x16 inches)\n\n")
cat("All files saved to:\n")
cat("  - figures/cytetype_analysis/\n")
cat("  - tables/cytetype_analysis/\n\n")
cat("================================================================================\n")

