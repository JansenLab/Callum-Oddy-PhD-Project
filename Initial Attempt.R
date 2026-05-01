install.packages(EnhandedVolcano)


# Load packages
library(Seurat)
library(SeuratObject)
library(harmony)
library(ggplot2)
library(dplyr)
library(EnhancedVolcano)
library(pheatmap)

# Set working directory
setwd("/Users/callu/OneDrive - University College London/Protocol/scRNA/")

set.seed(1234)

# Load matrices (your existing code)
mtx.6Gas <- Read10X_h5(filename = "./CO_06_GAS_GEX/filtered_feature_bc_matrix.h5")
mtx.6IM1 <- Read10X_h5(filename = "./CO_06_IM1_GEX/filtered_feature_bc_matrix.h5")
mtx.6IM2 <- Read10X_h5(filename = "./CO_06_IM2_GEX/filtered_feature_bc_matrix.h5")
mtx.6Duo <- Read10X_h5(filename = "./CO_06_DUO_GEX/filtered_feature_bc_matrix.h5")
mtx.7Gas <- Read10X_h5(filename = "./CO_07_GAS_GEX/filtered_feature_bc_matrix.h5")
mtx.7IM1 <- Read10X_h5(filename = "./CO_07_IM1_GEX/filtered_feature_bc_matrix.h5")
mtx.7IM2 <- Read10X_h5(filename = "./CO_07_IM2_GEX/filtered_feature_bc_matrix.h5")
mtx.7Duo <- Read10X_h5(filename = "./CO_07_DUO_GEX/filtered_feature_bc_matrix.h5")

# Create Seurat Objects (your existing code)
seu.6Gas <- CreateSeuratObject(mtx.6Gas, project = "CO_06_GAS")
seu.6IM1 <- CreateSeuratObject(mtx.6IM1, project = "CO_06_IM1")
seu.6IM2 <- CreateSeuratObject(mtx.6IM2, project = "CO_06_IM2")
seu.6Duo <- CreateSeuratObject(mtx.6Duo, project = "CO_06_DUO")
seu.7Gas <- CreateSeuratObject(mtx.7Gas, project = "CO_07_GAS")
seu.7IM1 <- CreateSeuratObject(mtx.7IM1, project = "CO_07_IM1")
seu.7IM2 <- CreateSeuratObject(mtx.7IM2, project = "CO_07_IM2")
seu.7Duo <- CreateSeuratObject(mtx.7Duo, project = "CO_07_DUO")

# Calculate mitochondrial %
seu.6Gas <- PercentageFeatureSet(seu.6Gas, pattern = "^MT-", col.name = "percent.mito")
seu.6IM1 <- PercentageFeatureSet(seu.6IM1, pattern = "^MT-", col.name = "percent.mito")
seu.6IM2 <- PercentageFeatureSet(seu.6IM2, pattern = "^MT-", col.name = "percent.mito")
seu.6Duo <- PercentageFeatureSet(seu.6Duo, pattern = "^MT-", col.name = "percent.mito")
seu.7Gas <- PercentageFeatureSet(seu.7Gas, pattern = "^MT-", col.name = "percent.mito")
seu.7IM1 <- PercentageFeatureSet(seu.7IM1, pattern = "^MT-", col.name = "percent.mito")
seu.7IM2 <- PercentageFeatureSet(seu.7IM2, pattern = "^MT-", col.name = "percent.mito")
seu.7Duo <- PercentageFeatureSet(seu.7Duo, pattern = "^MT-", col.name = "percent.mito")

# Filter cells (your filters)
seu.6Gas <- subset(seu.6Gas, nFeature_RNA > 400 & percent.mito < 10)
seu.6IM1 <- subset(seu.6IM1, nFeature_RNA > 400 & percent.mito < 10)
seu.6IM2 <- subset(seu.6IM2, nFeature_RNA > 400 & percent.mito < 10)
seu.6Duo <- subset(seu.6Duo, nFeature_RNA > 400 & percent.mito < 10)
seu.7Gas <- subset(seu.7Gas, nFeature_RNA > 400 & percent.mito < 10)
seu.7IM1 <- subset(seu.7IM1, nFeature_RNA > 400 & percent.mito < 10)
seu.7IM2 <- subset(seu.7IM2, nFeature_RNA > 400 & percent.mito < 10)
seu.7Duo <- subset(seu.7Duo, nFeature_RNA > 400 & percent.mito < 10)

# Merge samples
seu.int <- merge(seu.6Gas, y = c(seu.6IM1, seu.6IM2, seu.6Duo, seu.7Gas, seu.7IM1, seu.7IM2, seu.7Duo))

# Normalize and find variable features BEFORE Harmony
seu.int <- NormalizeData(seu.int)
seu.int <- FindVariableFeatures(seu.int)

# Scale data WITHOUT regression (we regress later after Harmony)
seu.int <- ScaleData(seu.int, features = VariableFeatures(seu.int))

# Run PCA
seu.int <- RunPCA(seu.int, features = VariableFeatures(seu.int))

# Add metadata for patient and tissue info
seu.int$patient <- sapply(strsplit(seu.int$orig.ident, "_"), function(x) paste(x[1], x[2], sep = "_"))
seu.int$tissue_type <- sapply(strsplit(seu.int$orig.ident, "_"), function(x) x[3])
seu.int$tissue_group <- ifelse(seu.int$tissue_type %in% c("IM1", "IM2"), "IM", seu.int$tissue_type)

# Run Harmony for batch correction by patient
seu.int <- RunHarmony(seu.int, group.by.vars = "patient")

# Cell cycle scoring BEFORE regression (using default Seurat gene lists)
seu.int <- CellCycleScoring(
  object = seu.int,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident = TRUE
)

# Regress out cell cycle scores during scaling using Harmony embeddings
seu.int <- ScaleData(
  seu.int,
  vars.to.regress = c("S.Score", "G2M.Score"),
  features = VariableFeatures(seu.int)
)

# Run UMAP and clustering on Harmony embeddings
seu.int <- RunUMAP(seu.int, reduction = "harmony", dims = 1:30)
seu.int <- FindNeighbors(seu.int, reduction = "harmony", dims = 1:30)
seu.int <- FindClusters(seu.int, resolution = 2)


# Visualize clusters and cell cycle phase
DimPlot(seu.int, reduction = "umap", label = TRUE, pt.size = 0.5) + ggtitle("Clusters (Harmony Batch Corrected)")
DimPlot(seu.int, reduction = "umap", group.by = "Phase", pt.size = 0.5) + ggtitle("Cell Cycle Phase")
DimPlot(seu.int, reduction = "umap", group.by = "patient", pt.size = 0.5) + ggtitle("Patient")
DimPlot(seu.int, reduction = "umap", group.by = "tissue_group", pt.size = 0.5) + ggtitle("Organoid Identity")


##ClusTree
install.packages("clustree")
library("clustree")
head(seu.int[[]])
# Ensure you have your graph built
seu.int <- FindNeighbors(seu.int, dims = 1:30)

# Run clustering across a range of resolutions
seu.int <- FindClusters(seu.int, resolution = c(0.2, 0.4, 0.6, 0.8, 1, 1.5, 2))

clustree(seu.int, prefix = "RNA_snn_res.")




library(Seurat)
library(dplyr)

# --- STEP 1: Prepare your marker gene list ---
lineage_markers <- list(
  "Gastric" = c("CLDN18","ANXA10"),
  "Intestinal" = c("CDX1","CDX2","CLDN3"),
  "Stem" = c("TNFRSF19","LGR5","EPHB2","OLFM4","SMOC2","ASCL2","MEX3A"),
  "IsthmusStem" = c("STMN1"),
  "Proliferating" = c("HELLS", "PCNA", "TOP2A","MKI67","BIRC5"),
  "Chief" = c("PGC","LIPF","PGA3"),
  "NeckMucous" = c("MUC6"),
  "GastricMucous" = c("TFF2","TFF1","MUC5AC"),
  "Parietal" = c("ATP4A"),
  "Enterocytes" = c("FABP1","KRT20","ANPEP"),
  "GobletPre" = c("HES6"),
  "Goblet" = c("SPINK4","ATOH1","MUC2"),
  "Paneth" = c("DEFA6","PLA2G2A"),
  "EE.com.pro" = c("PROX1","BMI1"),
  "Enteroendo" = c("CHGA","NEUROG3","SYP"),
  "Tuft" = c("POU2F3","AVIL")
)

# Load required libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(RColorBrewer)

# Function to create violin plots for each lineage
create_lineage_violin_plots <- function(seu_obj, lineage_markers, 
                                        group.by = "seurat_clusters",
                                        pt.size = 0.1,
                                        ncol = 2,
                                        save_plots = TRUE,
                                        output_dir = "violin_plots") {
  
  # Create output directory if saving plots
  if (save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Get all unique genes
  all_genes <- unique(unlist(lineage_markers))
  
  # Check which genes are present in the dataset
  genes_present <- intersect(all_genes, rownames(seu_obj))
  genes_missing <- setdiff(all_genes, rownames(seu_obj))
  
  if (length(genes_missing) > 0) {
    cat("Warning: The following genes are not present in the dataset:\n")
    cat(paste(genes_missing, collapse = ", "), "\n\n")
  }
  
  # Filter lineage_markers to only include present genes
  lineage_markers_filtered <- lapply(lineage_markers, function(genes) {
    intersect(genes, genes_present)
  })
  
  # Remove empty lineages
  lineage_markers_filtered <- lineage_markers_filtered[lengths(lineage_markers_filtered) > 0]
  
  # Create plots for each lineage
  lineage_plots <- list()
  
  for (lineage in names(lineage_markers_filtered)) {
    genes <- lineage_markers_filtered[[lineage]]
    
    cat(paste("Creating plots for", lineage, "lineage with", length(genes), "genes\n"))
    
    # Create individual violin plots for each gene
    gene_plots <- list()
    
    for (gene in genes) {
      p <- VlnPlot(seu_obj, 
                   features = gene,
                   group.by = group.by,
                   pt.size = pt.size) +
        ggtitle(paste(lineage, "-", gene)) +
        theme(axis.title.x = element_blank(),
              plot.title = element_text(size = 12, hjust = 0.5),
              legend.position = "none") +
        xlab("Cluster")
      
      gene_plots[[gene]] <- p
    }
    
    # Combine plots for this lineage
    if (length(gene_plots) == 1) {
      combined_plot <- gene_plots[[1]]
    } else {
      combined_plot <- wrap_plots(gene_plots, ncol = ncol)
    }
    
    # Add overall title
    lineage_plot <- combined_plot + 
      plot_annotation(title = paste(lineage, "Lineage Markers"),
                      theme = theme(plot.title = element_text(size = 16, hjust = 0.5)))
    
    lineage_plots[[lineage]] <- lineage_plot
    
    # Save individual lineage plots
    if (save_plots) {
      filename <- paste0(output_dir, "/", lineage, "_violin_plots.pdf")
      ggsave(filename, lineage_plot, width = 12, height = 8)
      cat(paste("Saved:", filename, "\n"))
    }
  }
  
  return(lineage_plots)
}

# Function to create a summary plot with all genes
create_summary_violin_plot <- function(seu_obj, lineage_markers,
                                       group.by = "seurat_clusters",
                                       genes_per_row = 4,
                                       pt.size = 0.05) {
  
  # Get all present genes
  all_genes <- unique(unlist(lineage_markers))
  genes_present <- intersect(all_genes, rownames(seu_obj))
  
  # Create a single plot with all genes
  if (length(genes_present) > 0) {
    summary_plot <- VlnPlot(seu_obj, 
                            features = genes_present,
                            group.by = group.by,
                            pt.size = pt.size,
                            ncol = genes_per_row) +
      plot_annotation(title = "All Lineage Markers - Expression by Cluster",
                      theme = theme(plot.title = element_text(size = 16, hjust = 0.5)))
    
    return(summary_plot)
  } else {
    stop("No genes found in the dataset")
  }
}

# Main execution
# Assuming your Seurat object is named 'seu.int'
if (exists("seu.int")) {
  
  cat("Creating violin plots for lineage markers...\n\n")
  
  # Create plots by lineage
  lineage_plots <- create_lineage_violin_plots(
    seu_obj = seu.int,
    lineage_markers = lineage_markers,
    group.by = "seurat_clusters", 
    pt.size = 0.1,
    ncol = 3,  
    save_plots = TRUE,
    output_dir = "violin_plots"
  )
  
  # Display plots, such as: lineage_plots$Gastric, lineage_plots$Stem, etc.
  
  # Create and save a summary plot with all genes
  cat("\nCreating summary plot with all genes...\n")
  summary_plot <- create_summary_violin_plot(
    seu_obj = seu.int,
    lineage_markers = lineage_markers,
    group.by = "seurat_clusters",
    genes_per_row = 4,
    pt.size = 0.05
  )
  
  # Save summary plot
  ggsave("violin_plots/all_markers_summary.pdf", summary_plot, 
         width = 16, height = 20)
  
  cat("All plots have been created and saved!\n")
  
} else {
  cat("Please make sure your Seurat object is named 'seu.int' or modify the script accordingly.\n")
}

###Heatmap Expression

create_expression_heatmap <- function(seu_obj, lineage_markers, 
                                      group.by = "seurat_clusters") {
  
  # Get average expression per cluster
  all_genes <- unique(unlist(lineage_markers))
  genes_present <- intersect(all_genes, rownames(seu_obj))
  
  if (length(genes_present) > 0) {
    
    # Check if genes are in scale.data, if not, scale them
    scaled_genes <- intersect(genes_present, rownames(GetAssayData(seu_obj, layer = "scale.data")))
    
    if (length(scaled_genes) < length(genes_present)) {
      cat("Scaling missing genes for heatmap...\n")
      # Scale the genes that aren't already scaled
      genes_to_scale <- setdiff(genes_present, scaled_genes)
      seu_obj <- ScaleData(seu_obj, features = genes_to_scale, verbose = FALSE)
    }
    
    # Create heatmap using DoHeatmap
    tryCatch({
      heatmap_plot <- DoHeatmap(seu_obj, 
                                features = genes_present,
                                group.by = group.by,
                                size = 3) +
        ggtitle("Average Expression Heatmap - All Lineage Markers") +
        theme(plot.title = element_text(hjust = 0.5))
      
      return(heatmap_plot)
      
    }, error = function(e) {
      
      # If DoHeatmap fails, create alternative heatmap using AggregateExpression
      cat("DoHeatmap failed, creating alternative heatmap using AggregateExpression...\n")
      
      # Use AggregateExpression (recommended in Seurat v5)
      avg_exp <- AggregateExpression(seu_obj, 
                                     features = genes_present,
                                     group.by = group.by,
                                     return.seurat = FALSE)
      
      # Convert to matrix and create heatmap
      library(pheatmap)
      library(RColorBrewer)
      
      # Get the RNA assay data
      exp_matrix <- as.matrix(avg_exp$RNA)
      
      # Create annotation for lineages
      gene_annotation <- data.frame(
        Gene = rownames(exp_matrix),
        Lineage = NA
      )
      
      for (lineage in names(lineage_markers)) {
        lineage_genes <- intersect(lineage_markers[[lineage]], rownames(exp_matrix))
        gene_annotation$Lineage[gene_annotation$Gene %in% lineage_genes] <- lineage
      }
      
      rownames(gene_annotation) <- gene_annotation$Gene
      gene_annotation$Gene <- NULL
      
      # Create color palette
      lineage_colors <- rainbow(length(unique(gene_annotation$Lineage)))
      names(lineage_colors) <- unique(gene_annotation$Lineage)
      
      annotation_colors <- list(Lineage = lineage_colors)
      
      # Create heatmap
      heatmap_plot <- pheatmap(exp_matrix,
                               scale = "row",
                               clustering_distance_rows = "euclidean",
                               clustering_distance_cols = "euclidean",
                               color = colorRampPalette(c("blue", "white", "red"))(100),
                               annotation_row = gene_annotation,
                               annotation_colors = annotation_colors,
                               main = "Average Expression Heatmap - All Lineage Markers",
                               fontsize = 10,
                               fontsize_row = 8,
                               fontsize_col = 10)
      
      return(heatmap_plot)
    })
    
  } else {
    stop("No genes found in the dataset")
  }
}

# Alternative function using only ggplot2 (no additional dependencies)
create_ggplot_heatmap <- function(seu_obj, lineage_markers, 
                                  group.by = "seurat_clusters") {
  
  all_genes <- unique(unlist(lineage_markers))
  genes_present <- intersect(all_genes, rownames(seu_obj))
  
  if (length(genes_present) > 0) {
    
    # Use AggregateExpression
    avg_exp <- AggregateExpression(seu_obj, 
                                   features = genes_present,
                                   group.by = group.by,
                                   return.seurat = FALSE)
    
    # Convert to long format for ggplot
    exp_matrix <- as.matrix(avg_exp$RNA)
    
    # Create long format data
    library(reshape2)
    heatmap_data <- melt(exp_matrix)
    colnames(heatmap_data) <- c("Gene", "Cluster", "Expression")
    
    # Add lineage information
    heatmap_data$Lineage <- NA
    for (lineage in names(lineage_markers)) {
      lineage_genes <- intersect(lineage_markers[[lineage]], genes_present)
      heatmap_data$Lineage[heatmap_data$Gene %in% lineage_genes] <- lineage
    }
    
    # Scale expression values
    heatmap_data <- heatmap_data %>%
      group_by(Gene) %>%
      mutate(Scaled_Expression = scale(Expression)[,1]) %>%
      ungroup()
    
    # Create heatmap
    heatmap_plot <- ggplot(heatmap_data, aes(x = Cluster, y = Gene, fill = Scaled_Expression)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                           midpoint = 0, name = "Scaled\nExpression") +
      facet_grid(Lineage ~ ., scales = "free_y", space = "free_y") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 8),
            strip.text.y = element_text(angle = 0),
            plot.title = element_text(hjust = 0.5)) +
      ggtitle("Average Expression Heatmap - All Lineage Markers") +
      xlab("Cluster") +
      ylab("Gene")
    
    return(heatmap_plot)
    
  } else {
    stop("No genes found in the dataset")
  }
}

heatmap_plot <- create_expression_heatmap(seu.int, lineage_markers)

heatmap_plot



####Feature plots of those lineage markers
# Load required libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# Define lineage markers
lineage_markers <- list(
  "Gastric" = c("CLDN18","ANXA10"),
  "Intestinal" = c("CDX1","CDX2","CLDN3"),
  "Stem" = c("TNFRSF19","LGR5","EPHB2","OLFM4","SMOC2","ASCL2","MEX3A"),
  "IsthmusStem" = c("STMN1"),
  "Proliferating" = c("HELLS", "PCNA", "TOP2A","MKI67","BIRC5"),
  "Chief" = c("PGC","LIPF","PGA3"),
  "NeckMucous" = c("MUC6"),
  "GastricMucous" = c("TFF2","TFF1","MUC5AC"),
  "Parietal" = c("ATP4A"),
  "Enterocytes" = c("FABP1","KRT20","ANPEP"),
  "GobletPre" = c("HES6"),
  "Goblet" = c("SPINK4","ATOH1","MUC2"),
  "Paneth" = c("DEFA6","PLA2G2A"),
  "EE.com.pro" = c("PROX1","BMI1"),
  "Enteroendo" = c("CHGA","NEUROG3","SYP"),
  "Tuft" = c("POU2F3","AVIL")
)

# Function to create feature plots for all lineage markers
create_lineage_feature_plots <- function(seurat_obj, lineage_markers, ncol = 4) {
  
  # Create individual plots for each lineage
  lineage_plots <- list()
  
  for (lineage_name in names(lineage_markers)) {
    genes <- lineage_markers[[lineage_name]]
    
    # Check which genes are present in the dataset
    available_genes <- genes[genes %in% rownames(seurat_obj)]
    missing_genes <- genes[!genes %in% rownames(seurat_obj)]
    
    if (length(missing_genes) > 0) {
      message(paste("Missing genes in", lineage_name, ":", paste(missing_genes, collapse = ", ")))
    }
    
    if (length(available_genes) > 0) {
      # Create feature plot for available genes
      p <- FeaturePlot(seurat_obj, 
                       features = available_genes,
                       ncol = min(length(available_genes), ncol),
                       reduction = "umap", # Change to "tsne" if using t-SNE
                       pt.size = 0.5,
                       order = TRUE) +
        plot_annotation(title = paste(lineage_name, "Markers"),
                        theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")))
      
      lineage_plots[[lineage_name]] <- p
    }
  }
  
  return(lineage_plots)
}

# Function to create a comprehensive feature plot with all genes
create_comprehensive_feature_plot <- function(seurat_obj, lineage_markers, ncol = 6) {
  
  # Flatten all genes into a single vector
  all_genes <- unlist(lineage_markers)
  
  # Check which genes are present
  available_genes <- all_genes[all_genes %in% rownames(seurat_obj)]
  missing_genes <- all_genes[!all_genes %in% rownames(seurat_obj)]
  
  if (length(missing_genes) > 0) {
    message(paste("Missing genes:", paste(missing_genes, collapse = ", ")))
  }
  
  if (length(available_genes) > 0) {
    # Create comprehensive feature plot
    p <- FeaturePlot(seurat_obj, 
                     features = available_genes,
                     ncol = ncol,
                     reduction = "umap", # Change to "tsne" if using t-SNE
                     pt.size = 0.3,
                     order = TRUE,
                     combine = TRUE) +
      plot_annotation(title = "All Lineage Markers",
                      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
    
    return(p)
  }
}

# Function to create violin plots for lineage markers
create_lineage_violin_plots <- function(seurat_obj, lineage_markers) {
  
  violin_plots <- list()
  
  for (lineage_name in names(lineage_markers)) {
    genes <- lineage_markers[[lineage_name]]
    available_genes <- genes[genes %in% rownames(seurat_obj)]
    
    if (length(available_genes) > 0) {
      p <- VlnPlot(seurat_obj, 
                   features = available_genes,
                   ncol = min(length(available_genes), 3),
                   pt.size = 0) +
        plot_annotation(title = paste(lineage_name, "Expression"),
                        theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")))
      
      violin_plots[[lineage_name]] <- p
    }
  }
  
  return(violin_plots)
}

# Function to create folders and save feature plots for each lineage
save_lineage_feature_plots <- function(seurat_obj, lineage_markers, 
                                       base_dir = "lineage_feature_plots",
                                       ncol = 3, plot_width = 12, plot_height = 8,
                                       image_format = "png", dpi = 300) {
  
  # Create base directory if it doesn't exist
  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE)
    message(paste("Created directory:", base_dir))
  }
  
  # Process each lineage
  for (lineage_name in names(lineage_markers)) {
    genes <- lineage_markers[[lineage_name]]
    
    # Create lineage-specific folder
    lineage_dir <- file.path(base_dir, lineage_name)
    if (!dir.exists(lineage_dir)) {
      dir.create(lineage_dir, recursive = TRUE)
      message(paste("Created directory:", lineage_dir))
    }
    
    # Check which genes are present
    available_genes <- genes[genes %in% rownames(seurat_obj)]
    missing_genes <- genes[!genes %in% rownames(seurat_obj)]
    
    if (length(missing_genes) > 0) {
      message(paste("Missing genes in", lineage_name, ":", paste(missing_genes, collapse = ", ")))
    }
    
    if (length(available_genes) > 0) {
      # Create feature plot for available genes
      p <- FeaturePlot(seurat_obj, 
                       features = available_genes,
                       ncol = ncol,
                       reduction = "umap", 
                       pt.size = 0.5,
                       order = TRUE) +
        plot_annotation(title = paste(lineage_name, "Markers"),
                        theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")))
      
      # Calculate dynamic plot dimensions based on number of genes
      n_genes <- length(available_genes)
      n_rows <- ceiling(n_genes / ncol)
      dynamic_height <- max(plot_height, n_rows * 3)
      dynamic_width <- max(plot_width, ncol * 4)  
      
      # Save the combined plot as PDF (for overview)
      pdf_path <- file.path(lineage_dir, paste0(lineage_name, "_markers_combined.pdf"))
      ggsave(pdf_path, p, width = dynamic_width, height = dynamic_height, device = "pdf")
      message(paste("Saved combined plot:", pdf_path))
      
      # Save individual gene plots as images
      for (gene in available_genes) {
        individual_plot <- FeaturePlot(seurat_obj, 
                                       features = gene,
                                       reduction = "umap",
                                       pt.size = 0.5,
                                       order = TRUE) +
          ggtitle(paste(lineage_name, "-", gene)) +
          theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
        
        individual_path <- file.path(lineage_dir, paste0(gene, "_individual.", image_format))
        ggsave(individual_path, individual_plot, 
               width = 6, height = 5, 
               device = image_format, dpi = dpi)
      }
      
      message(paste("Saved", length(available_genes), "individual", image_format, "plots for", lineage_name))
    }
  }
  
  message(paste("All plots saved in:", base_dir))
}

# Function to also create violin plots in the same folder structure
save_lineage_violin_plots <- function(seurat_obj, lineage_markers, 
                                      base_dir = "lineage_feature_plots",
                                      plot_width = 12, plot_height = 8,
                                      image_format = "png", dpi = 300) {
  
  # Process each lineage
  for (lineage_name in names(lineage_markers)) {
    genes <- lineage_markers[[lineage_name]]
    lineage_dir <- file.path(base_dir, lineage_name)
    
    # Check which genes are present
    available_genes <- genes[genes %in% rownames(seurat_obj)]
    
    if (length(available_genes) > 0) {
      # Create violin plot
      v <- VlnPlot(seurat_obj, 
                   features = available_genes,
                   ncol = min(length(available_genes), 3),
                   pt.size = 0) +
        plot_annotation(title = paste(lineage_name, "Expression Distribution"),
                        theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")))
      
      # Calculate dimensions
      n_genes <- length(available_genes)
      n_cols <- min(n_genes, 3)
      n_rows <- ceiling(n_genes / n_cols)
      dynamic_height <- max(plot_height, n_rows * 4)
      dynamic_width <- max(plot_width, n_cols * 5)
      
      # Save violin plot as image
      violin_path <- file.path(lineage_dir, paste0(lineage_name, "_violin.", image_format))
      ggsave(violin_path, v, 
             width = dynamic_width, height = dynamic_height, 
             device = image_format, dpi = dpi)
      message(paste("Saved violin plot:", violin_path))
    }
  }
}

# Set working directory and run analysis
setwd("/Users/callu/OneDrive - University College London/Protocol/scRNA/")

# Usage with your specific setup:

# 1. Create folders and save individual plots as PNG images (default)
save_lineage_feature_plots(seu.int, lineage_markers)

# 2. Also save violin plots as PNG images
save_lineage_violin_plots(seu.int, lineage_markers)

# 3. Alternative: Save as JPEG images with custom DPI
# save_lineage_feature_plots(seu.int, lineage_markers, image_format = "jpeg", dpi = 600)
# save_lineage_violin_plots(seu.int, lineage_markers, image_format = "jpeg", dpi = 600)

# 4. Alternative: Save as TIFF images (high quality)
# save_lineage_feature_plots(seu.int, lineage_markers, image_format = "tiff", dpi = 300)
# save_lineage_violin_plots(seu.int, lineage_markers, image_format = "tiff", dpi = 300)

# 5. Create comprehensive feature plot and save as image
comprehensive_plot <- create_comprehensive_feature_plot(seu.int, lineage_markers)
if (!dir.exists("lineage_feature_plots")) dir.create("lineage_feature_plots")
ggsave("lineage_feature_plots/comprehensive_all_markers.png", comprehensive_plot, 
       width = 20, height = 16, device = "png", dpi = 300)

 # Alternative: Create a simple feature plot for specific genes
# Example for just a few key markers:
# key_markers <- c("LGR5", "MKI67", "MUC2", "CHGA", "ATP4A")
# FeaturePlot(seurat_obj, features = key_markers, ncol = 3)

# Generate t-SNE for  integrated object

# 1. Run t-SNE (this will add tsne coordinates to object)
seu.int <- RunTSNE(seu.int, 
                   reduction = "pca",     
                   dims = 1:30,       
                   perplexity = 30,     
                   verbose = TRUE)


# 2. Basic t-SNE plot colored by clusters
tsne_clusters <- DimPlot(seu.int, 
                         reduction = "tsne", 
                         group.by = "seurat_clusters",
                         label = TRUE, 
                         label.size = 4) +
  ggtitle("t-SNE: Clusters") +
  theme(plot.title = element_text(hjust = 0.5))

# 3. t-SNE plot colored by sample (assuming sample info is in metadata)
tsne_samples <- DimPlot(seu.int, 
                        reduction = "tsne", 
                        group.by = "orig.ident",
                        label = FALSE) +
  ggtitle("t-SNE: Samples") +
  theme(plot.title = element_text(hjust = 0.5))

tsne_patient <- DimPlot(seu.int, 
                        reduction = "tsne", 
                        group.by = "patient",  
                        label = FALSE) +
  ggtitle("t-SNE: Samples") +
  theme(plot.title = element_text(hjust = 0.5))

tsne_phenotype <- DimPlot(seu.int, 
                          reduction = "tsne", 
                          group.by = "tissue_group",
                          label = FALSE) +
  ggtitle("t-SNE: Samples") +
  theme(plot.title = element_text(hjust = 0.5))

# 4. Display plots side by side
combined_plot1 <- tsne_clusters | tsne_samples
print(combined_plot1)

combined_plot2 <- tsne_patient | tsne_phenotype
print(combined_plot2)



# 5. Save the plots
ggsave("tsne_clusters.png", tsne_clusters, width = 8, height = 6, dpi = 300)
ggsave("tsne_samples.png", tsne_samples, width = 10, height = 6, dpi = 300)
ggsave("tsne_combined.png", combined_plot, width = 16, height = 6, dpi = 300)
