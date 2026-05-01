library(Seurat)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(harmony)

setwd("/Users/callu/OneDrive - University College London/Protocol/scRNA/")

# 1. Load Data
message("Loading Seurat object...")
seu <- readRDS("seu_cytotrace_complete.rds")

# Define consistent color palette for tissue groups
tissue_cols <- c("GAS" = "#CC6666", 
                 "IM" = "#9966CC",  
                 "DUO" = "#6699CC")

annot_cols <-c("Chief" = "#2277b0",
               "IsthSC" = "#16bcce",
               "Early PitEntero" = "#fe7d08",
               "Paneth" = "#a9caeb", 
               "Enterocyte" = "#229f65", 
               "Parietal" = "#fdbd7e", 
               "Gastric_Lineage" = "#d22827", 
               "Proliferating" = "#93e184", 
               "Gastric_Mucous" = "#aa3ef9", 
               "neck_mucous_IM" = "#fc9999", 
               "Goblet" = "#8e554a", 
               "neck_mucous_gastric" = "#c7addb",
               "Goblet_Precursor" = "#db7ac2", 
               "neck_mucous_proliferating" = "#cd988d", 
               "Intestinal_Lineage" = "#b3bc65",
               "neck_mucous_transitional" = "#fbb3d2")

# --- Dimensionality Reduction Plots ---

# UMAP grouped by clusters
DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + 
  ggtitle("UMAP by Seurat Clusters")

# UMAP grouped by tissue group using the custom palette
DimPlot(seu, reduction = "umap", group.by = "tissue_group", cols = tissue_cols) + 
  ggtitle("UMAP by Tissue Group")

# UMAP grouped by tissue group using the custom palette
DimPlot(seu, reduction = "umap", group.by = "final_annotation", cols = annot_cols) + 
  ggtitle("UMAP by Cell Annotation")

# --- Heatmaps ---
library(dplyr)
library(Seurat)
library(ggplot2)

# 1. Set the active identity to your final annotations
Idents(seu) <- "final_annotation"

# 2. Find distinguishing markers for each cell type
celltype_markers <- FindAllMarkers(
  seu, 
  only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25
)

# 3. Extract the top 10 genes per cell type, sorted by average log2 fold change
top10_celltype <- celltype_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

# 4. Generate the Heatmap
# We reduce the y-axis text size because 10 genes * (e.g., 15 cell types) = 150 rows!
DoHeatmap(seu, features = top10_celltype$gene, size = 3,
          group.colors = annot_cols) +        
  theme(axis.text.y = element_text(size = 5)) +
  ggtitle("Top 10 Defining Genes per Cell Type")

# Unsupervised: Top 10 Highly Variable Genes per cluster
Idents(seu) <- "seurat_clusters"
seu.markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- seu.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

DoHeatmap(seu, features = top10$gene, size = 4) + 
  theme(axis.text.y = element_text(size = 6))

# Extract metadata for ggplot
meta <- seu@meta.data

# --- Cell Composition by Tissue Group ---
ggplot(meta, aes(x = tissue_group, fill = final_annotation)) +
  geom_bar(position = "fill", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = annot_cols) +        
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  labs(y = "Percentage of Cells", x = "Tissue Group", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- Cell Composition by Patient ---
ggplot(meta, aes(x = patient, fill = final_annotation)) +
  geom_bar(position = "fill", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = annot_cols) +        
  scale_y_continuous(labels = scales::percent) +
  facet_grid(~tissue_group, scales = "free_x") + 
  theme_minimal() +
  labs(y = "Percentage of Cells", x = "Patient", fill = "Cell Type")



# 1. Define the lists from your input
stomach_genes <- c("PGA4", "GIF", "GKN1", "PGA3", "PGA5", "ATP4A", "GAST", "ATP4B", "LIPF", "CLDN18", "MUC5AC", "GKN2", "TFF1", "TFF2", "MUC6", "PGC", "CHIA", "VSIG1", "CTSE", "ANXA10", "AQP5", "PSCA", "BARX1", "HNF4A") # truncated for brevity in code, use your full list
colon_genes <- c("INSL5", "AQP8", "MEP1A", "PRAC1", "ISX", "CA1", "MUC12", "MS4A12", "NOX1", "GUCA2A", "PYY", "NAT2", "LYPD8", "CDH17", "CLCA1", "GPA33", "FABP1", "VIL1", "CDX2", "ZG16", "TFF3", "KRT20") # truncated for brevity

# 2. Filter genes to only those present in the Seurat Object
stomach_genes_filt <- stomach_genes[stomach_genes %in% rownames(seu)]
colon_genes_filt <- colon_genes[colon_genes %in% rownames(seu)]

# 3. Add Module Scores
seu <- AddModuleScore(
  object = seu,
  features = list(stomach_genes_filt),
  name = "Stomach_Score",
  ctrl = 5
)

seu <- AddModuleScore(
  object = seu,
  features = list(colon_genes_filt),
  name = "Colon_Score",
  ctrl = 5
)

# 4. Visualization: Scoring Distribution
VlnPlot(seu, features = c("Stomach_Score1", "Colon_Score1"), 
        group.by = "final_annotation", pt.size = 0, ncol = 1) & 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 5. Feature Plot: Spatial distribution of scores
FeaturePlot(seu, features = c("Stomach_Score1", "Colon_Score1"), blend = TRUE)

# Make sure your color palette is loaded from earlier
tissue_cols <- c("GAS" = "#CC6666",   
                 "IM" = "#9370B8",    
                 "DUO" = "#6B8DB8")   

# ---------------------------------------------------------
# 1. Feature Plot: Custom Colors for Gastric vs Colon
# ---------------------------------------------------------
FeaturePlot(seu, 
            features = c("Stomach_Score1", "Colon_Score1"), 
            blend = TRUE, 
            cols = c("lightgrey", "#CC6666", "#6B8DB8"), 
            pt.size = 0.5)

# ---------------------------------------------------------
# 2. Violin Plots: Split by Tissue Group, then Final Annotation
# ---------------------------------------------------------
library(ggplot2)

# Extract the metadata to plot with pure ggplot2 for better control
meta_df <- seu@meta.data

# -- Plot A: Gastric Score --
ggplot(meta_df, aes(x = final_annotation, y = Stomach_Score1, fill = tissue_group)) +
  geom_violin(scale = "width", trim = FALSE, color = "black", linewidth = 0.2) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  facet_wrap(~ tissue_group, scales = "free_x") + 
  scale_fill_manual(values = tissue_cols) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none", 
        strip.background = element_rect(fill = "lightgrey", color = "black")) +
  labs(title = "Gastric Gene Signature Expression",
       subtitle = "Divided by Organoid Phenotype",
       x = "Cell Annotation",
       y = "Module Score")

# -- Plot B: Colon Score --
ggplot(meta_df, aes(x = final_annotation, y = Colon_Score1, fill = tissue_group)) +
  geom_violin(scale = "width", trim = FALSE, color = "black", linewidth = 0.2) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  facet_wrap(~ tissue_group, scales = "free_x") + 
  scale_fill_manual(values = tissue_cols) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        strip.background = element_rect(fill = "lightgrey", color = "black")) +
  labs(title = "Intestinal/Colon Gene Signature Expression",
       subtitle = "Divided by Organoid Phenotype",
       x = "Cell Annotation",
       y = "Module Score")

###########
#Signature Plots
###########
library(ggplot2)

# Extract metadata from Seurat object
meta_df <- seu@meta.data

# Define the specific muted color palette for the tissue groups
tissue_cols <- c("GAS" = "#CC6666",  
                 "IM"  = "#9966CC",  
                 "DUO" = "#6699CC") 

# Generate the faceted scatter plot
ggplot(meta_df, aes(x = Stomach_Score1, y = Colon_Score1, color = tissue_group)) +
  geom_point(size = 0.8, alpha = 0.6) +
  
  # Split the plot into 3 panels based on the organoid type
  facet_wrap(~ tissue_group) +
  
  # Apply the muted color palette
  scale_color_manual(values = tissue_cols) +
  
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  ) +
  
  # Meaningful labels
  labs(
    title = "Lineage Program Co-expression",
    subtitle = "Comparing Gastric vs. Intestinal signatures at single-cell resolution",
    x = "Gastric Signature Score",
    y = "Intestinal Signature Score"
  )




intestinal_markers <- c(
  # Master intestinal transcription factors
  "CDX1", "CDX2",
  # Enterocyte / absorptive lineage
  "VIL1", "ALPI", "FABP1", "FABP2",
  "APOA1", "APOA4",
  "SI", "SLC26A3", "SLC5A1",
  # Goblet cell / secretory
  "MUC2", "SPDEF", "AGR2",
  # Paneth cell
  "LYZ", "DEFA5", "DEFA6",
  # Intestinal stem / Wnt
  "LGR5", "ASCL2", "OLFM4",
  # Intestinal epithelial adhesion
  "CDH17"
)

# ---------------------------
# GASTRIC MARKERS
# ---------------------------

gastric_markers <- c(
  # Gastric lineage transcription factors
  "SOX2", "GATA4", "GATA6",
  # Gastric pit / foveolar
  "MUC5AC", "TFF1", "GKN1", "GKN2",
  # Gastric gland / neck
  "MUC6", "TFF2", "BHLHA15",
  # Chief cell
  "PGC", "GIF", "CPA1",
  # Parietal cell
  "ATP4A", "ATP4B", "SLC26A7",
  # Gastric epithelial identity
  "CLDN18", "KRT18"
)

# ---------------------------
# INTESTINAL METAPLASIA (IM) MARKERS
# ---------------------------
IM_markers <- c(
  # Core IM transcriptional drivers
  "CDX2", "SOX9", "HNF4A",
  # Hybrid gastric–intestinal mucins
  "MUC2", "MUC1", "MUC5AC", "MUC6",
  # Intestinal stem / progenitor features 
  "LGR5", "OLFM4",
  # Intestinal epithelial identity 
  "CDH17", "KRT20",
  # IM-associated secretory / trefoil factors
  "TFF3",
  # Hep antigen (IM / Barrett’s / intestinal differentiation)
  "ANXA10",
  # Metaplasia-associated differentiation
  "AGR2", "SPDEF", "REG4",
  # Wnt / regeneration
  "AXIN2", "RNF43",
  # Tight junction / IM-associated claudins
  "CLDN7", "CLDN10"
)





# --- 1. Data Extraction from Seurat ---
# Using 'data' slot (normalized counts) instead of raw counts
marker_genes <- unique(c(gastric_markers, IM_markers, intestinal_markers))

# Only keep genes that actually exist in the dataset
valid_markers <- marker_genes[marker_genes %in% rownames(seu)]

# Pull data and scale it
exp_mat <- as.matrix(GetAssayData(seu, layer = "data")[valid_markers, ])
heatmap_mat_z <- t(scale(t(exp_mat)))

# --- 2. Improved Marker Categorization ---
gene_categories <- data.frame(gene = rownames(heatmap_mat_z)) %>%
  mutate(Set = case_when(
    gene %in% gastric_markers ~ "Gastric",
    gene %in% IM_markers      ~ "IM", 
    gene %in% intestinal_markers ~ "Intestinal",
    TRUE ~ "Other"
  )) %>%
  mutate(Set = factor(Set, levels = c("Gastric", "IM", "Intestinal")))

# Sort the matrix by these categories
heatmap_mat_z <- heatmap_mat_z[order(gene_categories$Set), ]
row_ann <- data.frame(Marker_Type = gene_categories$Set[order(gene_categories$Set)])
rownames(row_ann) <- rownames(heatmap_mat_z)

# --- 3. Column (Cell/Sample) Annotation ---
# If plotting every cell is too "busy", we can average by cluster/phenotype
# But for organoids, showing the cells usually looks great.
annotation_col <- seu@meta.data[, c("tissue_group", "patient")]
colnames(annotation_col) <- c("Phenotype", "Patient")

# --- 4. Plotting with specific Hex Codes ---
library(pheatmap)
ann_colors <- list(
  Phenotype = c(
    "GAS" = "#C75A5A",
    "IM"  = "#9370B8",
    "DUO" = "#6B8DB8"
  ),
  Patient = c(
    "CO_06" = "#E69F00", 
    "CO_07" = "#56B4E9"
  ),
  Marker_Type = c(
    "Gastric" = "#C75A5A",
    "IM" = "#9370B8",
    "Intestinal" = "#6B8DB8"
  )
)

# 1. Set Identity to your final annotations
Idents(seu) <- "final_annotation"

# 2. Calculate average expression for the markers
avg_exp <- AverageExpression(seu, 
                             features = valid_markers, 
                             return.seurat = FALSE, 
                             slot = "data")$RNA

# 3. Z-score the data
avg_z <- t(scale(t(avg_exp)))

# 4. Prepare Row Annotation (The Gene Groups: Gastric, IM, Intestinal)
gene_categories <- data.frame(gene = rownames(avg_z)) %>%
  mutate(Marker_Type = case_when(
    gene %in% gastric_markers ~ "Gastric",
    gene %in% IM_markers      ~ "IM",
    gene %in% intestinal_markers ~ "Intestinal"
  )) %>%
  mutate(Marker_Type = factor(Marker_Type, levels = c("Gastric", "IM", "Intestinal")))

row_ann <- data.frame(Marker_Type = gene_categories$Marker_Type)
rownames(row_ann) <- rownames(avg_z)
row_gaps <- cumsum(table(row_ann$Marker_Type))

# 5. Define annotation colors for the Row Markers
ann_colors <- list(
  Marker_Type = c(
    "Gastric" = "#C75A5A",
    "IM" = "#9370B8",
    "Intestinal" = "#6B8DB8"
  )
)

# 6. Final Heatmap Call
pheatmap(
  avg_z,
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100),
  annotation_row = row_ann,
  annotation_colors = ann_colors,
  cluster_rows = FALSE,
  cluster_cols = TRUE,    
  show_colnames = TRUE,   
  border_color = NA,
  breaks = seq(-2, 2, length.out = 101), 
  gaps_row = row_gaps,
  main = "Lineage Signature: Average Expression by Cell Type"
)




# 1. Set the identity
Idents(seu) <- "final_annotation"

# 2. Generate the plot
main_umap <- DimPlot(seu, 
                     reduction = "umap", 
                     label = TRUE,        
                     label.size = 3.5,    
                     repel = TRUE,        
                     pt.size = 0.5) +
  theme_void() + 
  theme(
    legend.position = "right",            
    legend.text = element_text(size = 8), 
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) + 
ggtitle("Global Cellular Manifold by Annotated Cell Type")

# Display the plot
main_umap






library(ggplot2)
library(dplyr)

# 1. Prepare the data: Calculate percentage of tissue groups within each cell type
comp_data <- seu@meta.data %>%
  group_by(final_annotation, tissue_group) %>%
  tally() %>%
  group_by(final_annotation) %>%
  mutate(pct = n / sum(n) * 100)

# 2. Define your specific muted color palette
tissue_cols <- c("GAS" = "#CC6666", 
                 "IM"  = "#9370B8", 
                 "DUO" = "#6B8DB8")

# 3. Generate the plot
ggplot(comp_data, aes(x = final_annotation, y = pct, fill = tissue_group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = tissue_cols, name = "Organoid Type") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "Cell Type Origin Distribution",
    subtitle = "Relative contribution of Gastric, IM, and Duodenal organoids to each cell type",
    x = "Annotated Cell Type",
    y = "Percentage of Cells (%)"
  )



########
#GSEA Enrichment Analysis
#######

Idents(seu) <- "final_annotation"

message("Re-generating markers for all cell types. This may take a moment...")

# Find markers for each of cluster (final_annotation)
celltype_markers <- FindAllMarkers(
  seu, 
  only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25
)

# Quick check to ensure it worked
message(paste("Markers found:", nrow(celltype_markers)))
head(celltype_markers)

# STEP 1: Relax the fold-change threshold to 0.25
go_markers <- celltype_markers %>%
  filter(avg_log2FC > 0.25) %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC)

# -------------------------------------------------------------------
# STEP 2: DEFINE YOUR CLUSTERS
# -------------------------------------------------------------------
all_clusters <- c("neck_mucous_gastric", "Goblet_Precursor", "neck_mucous_proliferating", 
                  "Chief", "neck_mucous_transitional", "IsthSC", "Parietal", 
                  "Proliferating", "Gastric_Mucous", "Gastric_Lineage", 
                  "neck_mucous_IM", "Goblet", "Intestinal_Lineage", 
                  "Enterocyte", "Paneth", "Early PitEntero")

# Create lists to store the raw data and the plots
go_results_list <- list()
plot_list <- list()

# -------------------------------------------------------------------
# STEP 3: RUN THE AUTOMATED ENRICHMENT LOOP
# -------------------------------------------------------------------
for (cell_type in all_clusters) {
  
  message(paste("Running GO enrichment for:", cell_type))
  
  # Extract and clean the genes for this specific cluster
  cluster_genes <- go_markers$gene[go_markers$cluster == cell_type]
  clean_genes <- unique(toupper(cluster_genes))
  
  # Only proceed if we have at least 5 genes
  if (length(clean_genes) >= 5) {
    
    # tryCatch prevents the entire loop from crashing if one cluster fails
    ego <- tryCatch({
      enrichGO(gene          = clean_genes,
               OrgDb         = org.Hs.eg.db,
               keyType       = 'SYMBOL',
               ont           = "BP", 
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.05)
    }, error = function(e) {
      message(paste("  -> Error running enrichGO for", cell_type, ":", e$message))
      return(NULL)
    })
    
    # If successful and found enriched terms, save the result and create a plot
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      
      # Save the raw data table
      go_results_list[[cell_type]] <- as.data.frame(ego)
      
      # Create the plot and save it to the plot list
      p <- dotplot(ego, showCategory = 10) + 
        ggtitle(paste("GO Biological Process:", cell_type)) +
        theme_minimal()
      plot_list[[cell_type]] <- p
      
    } else {
      message(paste("  -> No significant GO terms found for", cell_type))
    }
    
  } else {
    message(paste("  -> Skipping", cell_type, "- not enough marker genes (", length(clean_genes), ")"))
  }
}

# -------------------------------------------------------------------
# STEP 4: EXPORT RESULTS FOR THE THESIS
# -------------------------------------------------------------------
# Save all plots into a single PDF 
pdf("GO_Enrichment_All_Clusters.pdf", width = 10, height = 8)
for (p in plot_list) {
  print(p)
}
dev.off()
message("Saved all GO plots to 'GO_Enrichment_All_Clusters.pdf'")

# Combine all text results into one big table and save to CSV
if (length(go_results_list) > 0) {
  final_go_table <- bind_rows(go_results_list, .id = "Cell_Type")
  write.csv(final_go_table, "GO_Enrichment_Results_All.csv", row.names = FALSE)
  message("Saved all GO data to 'GO_Enrichment_Results_All.csv'")
}


#####################################
#KEGG
#####################################
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

# Create lists to store the KEGG data and plots
kegg_results_list <- list()
kegg_plot_list <- list()

message("Starting KEGG Enrichment Pipeline...")

# Loop through every single cell type
for (cell_type in all_clusters) {
  
  message(paste("Running KEGG enrichment for:", cell_type))
  
  # 1. Extract and clean the Gene Symbols for this cluster
  cluster_genes <- go_markers$gene[go_markers$cluster == cell_type]
  clean_symbols <- unique(toupper(cluster_genes))
  
  # Only proceed if we have at least 5 genes
  if (length(clean_symbols) >= 5) {
    
    # 2. Convert SYMBOL to ENTREZID 
    gene_translation <- suppressMessages(tryCatch({
      bitr(clean_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    }, error = function(e) { return(NULL) }))
    
    # If translation worked, proceed to KEGG
    if (!is.null(gene_translation) && nrow(gene_translation) > 0) {
      entrez_ids <- gene_translation$ENTREZID
      
      # 3. Run KEGG Enrichment for Human ('hsa')
      ekegg <- tryCatch({
        enrichKEGG(gene         = entrez_ids,
                   organism     = 'hsa', 
                   pvalueCutoff = 0.05)
      }, error = function(e) {
        message(paste("  -> Error running enrichKEGG for", cell_type, ":", e$message))
        return(NULL)
      })
      
      # If successful and found enriched pathways
      if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
        
        # 4. Convert Entrez IDs back to Gene Symbols for the final table
        ekegg_readable <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
        
        # Save the raw data table
        kegg_results_list[[cell_type]] <- as.data.frame(ekegg_readable)
        
        # Create the plot and save it
        p <- dotplot(ekegg_readable, showCategory = 10) + 
          ggtitle(paste("KEGG Pathways:", cell_type)) +
          theme_minimal()
        kegg_plot_list[[cell_type]] <- p
        
      } else {
        message(paste("  -> No significant KEGG pathways found for", cell_type))
      }
    } else {
      message(paste("  -> Skipping", cell_type, "- ID translation failed."))
    }
  } else {
    message(paste("  -> Skipping", cell_type, "- not enough marker genes (", length(clean_symbols), ")"))
  }
}

# -------------------------------------------------------------------
# EXPORT RESULTS FOR THE THESIS
# -------------------------------------------------------------------
# Save all plots into a single PDF 
if (length(kegg_plot_list) > 0) {
  pdf("KEGG_Enrichment_All_Clusters.pdf", width = 10, height = 8)
  for (p in kegg_plot_list) {
    print(p)
  }
  dev.off()
  message("Saved all KEGG plots to 'KEGG_Enrichment_All_Clusters.pdf'")
}

# Combine all text results into one big table and save to CSV
if (length(kegg_results_list) > 0) {
  final_kegg_table <- bind_rows(kegg_results_list, .id = "Cell_Type")
  write.csv(final_kegg_table, "KEGG_Enrichment_Results_All.csv", row.names = FALSE)
  message("Saved all KEGG data to 'KEGG_Enrichment_Results_All.csv'")
}















#####################
#Plotting out cytotrace2 scores
#####################
ggplot(meta_df, aes(x = final_annotation, y = CytoTRACE2_Score, fill = final_annotation)) +
  geom_jitter(aes(color = final_annotation), width = 0.2, size = 0.3, alpha = 0.4) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black", linewidth = 0.3, fill = NA) +
  scale_fill_manual(values = annot_cols) +
  scale_color_manual(values = annot_cols) +
  scale_y_continuous(
    limits = c(0, 1),
    sec.axis = sec_axis(~ .,
                        breaks = c(0.9, 0.7, 0.5, 0.3, 0.1),
                        labels = c("Pluripotent", "Multipotent", "Oligopotent",
                                   "Unipotent", "Differentiated"))
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "none",
        axis.text.y.right = element_text(size = 8, color = "grey40")) +
  labs(title = "Developmental potential by phenotype",
       x = "Phenotype",
       y = "Potency score")

# Use the existing CytoTRACE2_Potency column directly
potency_order <- c("Differentiated", "Unipotent", "Oligopotent", "Multipotent", "Pluripotent")
meta_df$CytoTRACE2_Potency <- factor(meta_df$CytoTRACE2_Potency, levels = potency_order)

prop_df <- meta_df %>%
  group_by(CytoTRACE2_Potency, final_annotation) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(CytoTRACE2_Potency) %>%
  mutate(pct = n / sum(n) * 100)

ggplot(prop_df, aes(x = CytoTRACE2_Potency, y = pct, fill = final_annotation)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = annot_cols) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     breaks = seq(0, 100, 25)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm")) +
  guides(fill = guide_legend(ncol = 3)) +
  labs(title = "Lineage Composition within Potency Categories",
       x = "CytoTRACE 2 Potency",
       y = "Percentage of Cells",
       fill = "Lineage")