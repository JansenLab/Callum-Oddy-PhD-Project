############################################
# Gene Promoter Methylation Analysis - PAIRED STATISTICS
# Completely rebuilt from scratch with correct paired testing
############################################


library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggsignif)
library(dplyr)

############################################
# FUNCTION: Analyze single gene with CORRECT paired statistics
############################################

analyze_gene_promoter <- function(gene_name, 
                                  beta_matrix,
                                  annotation,
                                  pheno_data,
                                  promoter_regions = c("TSS1500", "TSS200"),
                                  output_dir = "Gene_Methylation_Plots",
                                  save_plot = TRUE,
                                  show_stats = TRUE) {
  
  if (save_plot && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("\n========================================\n")
  cat("Analyzing:", gene_name, "\n")
  cat("========================================\n")
  
  # 1. Find gene probes in annotation
  gene_indices <- grep(paste0("\\b", gene_name, "\\b"), annotation$UCSC_RefGene_Name)
  
  if (length(gene_indices) == 0) {
    cat("⚠️  No probes found for", gene_name, "\n")
    return(NULL)
  }
  
  # 2. Filter for promoter probes
  gene_anno <- annotation[gene_indices, ]
  is_promoter <- grepl(paste(promoter_regions, collapse="|"), gene_anno$UCSC_RefGene_Group)
  gene_promoter_probes <- gene_anno$Name[is_promoter]
  
  if (length(gene_promoter_probes) == 0) {
    cat("⚠️  No promoter probes found for", gene_name, "\n")
    return(NULL)
  }
  
  # 3. Check which probes exist in beta matrix
  available_probes <- gene_promoter_probes[gene_promoter_probes %in% rownames(beta_matrix)]
  
  if (length(available_probes) == 0) {
    cat("⚠️  No promoter probes in beta matrix for", gene_name, "\n")
    return(NULL)
  }
  
  cat("Found", length(available_probes), "promoter probes\n")
  
  # 4. Calculate mean beta across promoter probes for each sample
  beta_gene_promoter <- beta_matrix[available_probes, , drop = FALSE]
  gene_promoter_means <- colMeans(beta_gene_promoter, na.rm = TRUE)
  
  # 5. Create data frame with proper ordering
  plot_data <- pheno_data %>%
    mutate(
      Gene_Promoter_Mean = gene_promoter_means,
      Sample_Type = factor(Sample_Type, levels = c("Gastric", "IM", "Duodenal"))
    )
  
  # 6. Mixed-effects model accounting for patient pairing
  cat("\nFitting mixed-effects model...\n")
  model <- lmer(Gene_Promoter_Mean ~ Sample_Type + (1|Patient), data = plot_data)
  
  if (show_stats) {
    cat("\n--- Model Summary ---\n")
    print(summary(model))
  }
  
  # 7. Pairwise comparisons with Tukey adjustment
  cat("\n--- Pairwise Comparisons (Tukey-adjusted) ---\n")
  emm <- emmeans(model, specs = "Sample_Type")
  pairs_result <- pairs(emm, adjust = "tukey")
  comparisons_df <- as.data.frame(pairs_result)
  print(comparisons_df)
  
  # 8. Extract p-values in the order we need for plotting
  # emmeans gives us: Gastric-IM, Gastric-Duodenal, IM-Duodenal
  # But we need to flip signs for plotting (group1 vs group2)
  
  # Find each comparison
  gastric_vs_im <- comparisons_df[comparisons_df$contrast == "Gastric - IM", ]
  gastric_vs_duo <- comparisons_df[comparisons_df$contrast == "Gastric - Duodenal", ]
  im_vs_duo <- comparisons_df[comparisons_df$contrast == "IM - Duodenal", ]
  
  # Create plotting annotations
  max_y <- max(plot_data$Gene_Promoter_Mean, na.rm = TRUE)
  
  plot_annotations <- data.frame(
    group1 = c("Gastric", "Gastric", "IM"),
    group2 = c("IM", "Duodenal", "Duodenal"),
    p_value = c(
      gastric_vs_im$p.value,
      gastric_vs_duo$p.value,
      im_vs_duo$p.value
    ),
    y_position = c(
      max_y * 1.08,
      max_y * 1.18,
      max_y * 1.13
    )
  )
  
  # Add significance labels
  plot_annotations$label <- sapply(plot_annotations$p_value, function(p) {
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    return("ns")
  })
  
  cat("\n--- Plotting Annotations ---\n")
  print(plot_annotations)
  
  # 9. Define colors
  tissue_colors <- c(
    "Gastric"   = "#c85a5a",
    "IM"        = "#8e7cc3",
    "Duodenal"  = "#487cac"
  )
  
  # 10. Create plot
  p <- ggplot(plot_data, aes(x = Sample_Type, y = Gene_Promoter_Mean)) +
    
    # Boxplots
    geom_boxplot(
      aes(fill = Sample_Type),
      outlier.shape = NA,
      alpha = 0.4,
      width = 0.5
    ) +
    
    # Individual points
    geom_point(
      aes(fill = Patient),
      size = 3,
      shape = 21,
      alpha = 0.8
    ) +
    
    # Lines connecting paired samples
    geom_line(
      aes(group = Patient, color = Patient),
      alpha = 0.6,
      linewidth = 0.8
    ) +
    
    # Significance brackets
    geom_signif(
      data = plot_annotations,
      aes(
        xmin = group1,
        xmax = group2,
        annotations = label,
        y_position = y_position
      ),
      manual = TRUE,
      tip_length = 0.02,
      textsize = 5,
      vjust = -0.2,
      size = 0.5
    ) +
    
    # Colors
    scale_fill_manual(
      values = tissue_colors,
      aesthetics = "fill",
      breaks = c("Gastric", "IM", "Duodenal")
    ) +
    scale_color_brewer(palette = "Set3") +
    
    # Theme and labels
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
      axis.title = element_text(face = "bold"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste0(gene_name, " Promoter Methylation"),
      subtitle = "Lines connect matched samples from same patient",
      y = "Mean Beta Value (Promoter)",
      x = "Tissue Type",
      fill = "Tissue Type",
      color = "Patient ID"
    ) +
    scale_y_continuous(
      limits = c(0, max_y * 1.25),
      breaks = seq(0, 1, 0.2)
    )
  
  # 11. Save outputs
  if (save_plot) {
    # PNG
    png_file <- file.path(output_dir, paste0(gene_name, "_Promoter_Methylation.png"))
    ggsave(png_file, p, width = 10, height = 7, dpi = 300)
    cat("✅ Plot saved:", png_file, "\n")
    
    # PDF
    pdf_file <- file.path(output_dir, paste0(gene_name, "_Promoter_Methylation.pdf"))
    ggsave(pdf_file, p, width = 10, height = 7)
    
    # Statistics file
    stats_file <- file.path(output_dir, paste0(gene_name, "_Statistics.txt"))
    sink(stats_file)
    cat("=====================================\n")
    cat("Gene:", gene_name, "\n")
    cat("Number of promoter probes:", length(available_probes), "\n")
    cat("Promoter regions:", paste(promoter_regions, collapse = ", "), "\n")
    cat("=====================================\n\n")
    cat("--- Probe IDs ---\n")
    print(available_probes)
    cat("\n--- Mixed-Effects Model Summary ---\n")
    print(summary(model))
    cat("\n--- Pairwise Comparisons (Tukey-adjusted) ---\n")
    print(comparisons_df)
    cat("\n--- Plotting Annotations (matching plot) ---\n")
    print(plot_annotations)
    sink()
    cat("✅ Statistics saved:", stats_file, "\n")
  }
  
  # 12. Return results
  return(list(
    plot = p,
    model = model,
    comparisons = comparisons_df,
    plot_annotations = plot_annotations,
    data = plot_data,
    n_probes = length(available_probes),
    probe_ids = available_probes
  ))
}

############################################
# FUNCTION: Batch process multiple genes
############################################

analyze_multiple_genes <- function(gene_list,
                                   beta_matrix,
                                   annotation,
                                   pheno_data,
                                   promoter_regions = c("TSS1500", "TSS200"),
                                   output_dir = "Gene_Methylation_Plots",
                                   create_summary = TRUE) {
  
  cat("\n🔬 BATCH ANALYSIS: Processing", length(gene_list), "genes\n")
  cat("Output directory:", output_dir, "\n\n")
  
  results <- list()
  successful <- c()
  failed <- c()
  
  for (gene in gene_list) {
    result <- tryCatch({
      analyze_gene_promoter(
        gene_name = gene,
        beta_matrix = beta_matrix,
        annotation = annotation,
        pheno_data = pheno_data,
        promoter_regions = promoter_regions,
        output_dir = output_dir,
        save_plot = TRUE,
        show_stats = FALSE
      )
    }, error = function(e) {
      cat("❌ Error processing", gene, ":", conditionMessage(e), "\n")
      return(NULL)
    })
    
    if (!is.null(result)) {
      results[[gene]] <- result
      successful <- c(successful, gene)
    } else {
      failed <- c(failed, gene)
    }
  }
  
  # Summary report
  if (create_summary) {
    summary_file <- file.path(output_dir, "BATCH_SUMMARY.txt")
    sink(summary_file)
    cat("=====================================\n")
    cat("BATCH ANALYSIS SUMMARY\n")
    cat("=====================================\n")
    cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    cat("Total genes:", length(gene_list), "\n")
    cat("Successful:", length(successful), "\n")
    cat("Failed:", length(failed), "\n\n")
    
    if (length(successful) > 0) {
      cat("--- Successfully Processed ---\n")
      for (gene in successful) {
        cat(sprintf("  %s (%d probes)\n", gene, results[[gene]]$n_probes))
      }
    }
    
    if (length(failed) > 0) {
      cat("\n--- Failed ---\n")
      for (gene in failed) {
        cat(" ", gene, "\n")
      }
    }
    
    cat("\n--- Pairwise Comparisons Summary ---\n\n")
    for (gene in successful) {
      cat("=====================================\n")
      cat("Gene:", gene, "\n")
      print(results[[gene]]$plot_annotations)
      cat("\n")
    }
    sink()
    
    cat("\n✅ Summary saved:", summary_file, "\n")
  }
  
  cat("\n=====================================\n")
  cat("BATCH COMPLETE\n")
  cat("=====================================\n")
  cat("Successful:", length(successful), "genes\n")
  cat("Failed:", length(failed), "genes\n")
  cat("Output:", output_dir, "\n\n")
  
  return(results)
}

############################################
# USAGE
############################################

# Gene lists
genes_of_interest <- c(
  "ATP4A","BHLHA15","CLDN18","GATA4","GATA6","GIF","GKN1","GKN2","KRT18",
  "MUC5AC","MUC6","PGC","SLC26A7","SOX2","TFF1","TFF2","AGR2","ANXA10","AXIN2",
  "CDH17","CDX2","CLDN10","CLDN7","HNF4A","KRT20","LGR5","MUC1","MUC2","OLFM4",
  "REG4","RNF43","SOX9","SPDEF","TFF3","ALPI","APOA1","APOA4","ASCL2","CDX1",
  "FABP1","FABP2","LYZ","SI","SLC26A3","SLC5A1","VIL1"
)

genes_of_interest_from_RNAseq <- c(
  "A1BG","BGN","BLOC1S5-TXNDC5","CGB5","EPHA2","LIF",
  "LOC102723901", "MAP3K5", "PIK3CG", "SLC2A6","TMX2-CTNND1A",
  "WNT7A", "ZNF625-ZNF20","ZNF664-FAM101A", "HOXA1","ITGB2-AS1",
  "NKX6-3","PCSK1N","ZNF462","CA12","ADAMTSL4","C15orf52","C1orf21","CASC9",
  "COL4A2","FGF19","FLG","GABRP","GRHL3","KRT80","MXRA5","MYEOV","NLRP6","REG1A",
  "SERPINB7","SLC4A3","STAC","TRPM2","ADAMTS14","ALDH1A3","ALDH1A3","AZGP1",
  "CLDN10","GABRE","GCT6","PLXNB3","SHISA6","SIX1","TM6SF2"
)

# Test single gene first
test_result <- analyze_gene_promoter(
  gene_name = "CDX2",
  beta_matrix = myNorm,
  annotation = annEPICv2,
  pheno_data = Sample_Data
)
print(test_result$plot)

# Batch process all genes
results1 <- analyze_multiple_genes(
  gene_list = genes_of_interest,
  beta_matrix = myNorm,
  annotation = annEPICv2,
  pheno_data = Sample_Data,
  output_dir = "Gene_Methylation_Plots/GenesOfInterest"
)

results2 <- analyze_multiple_genes(
  gene_list = genes_of_interest_from_RNAseq,
  beta_matrix = myNorm,
  annotation = annEPICv2,
  pheno_data = Sample_Data,
  output_dir = "Gene_Methylation_Plots/RNAseq_Genes"
)


############################################
# Gene Promoter Methylation Heatmaps
# Creates sample-level heatmaps with unsupervised clustering
############################################

library(pheatmap)
library(RColorBrewer)

############################################
# FUNCTION: Create sample-level heatmap
############################################

create_sample_heatmap <- function(gene_list,
                                  beta_matrix,
                                  annotation,
                                  pheno_data,
                                  promoter_regions = c("TSS1500", "TSS200"),
                                  output_file = "Gene_Methylation_Heatmap.pdf",
                                  output_dir = "Gene_Methylation_Plots",
                                  title = "Promoter Methylation Heatmap",
                                  width = 10,
                                  height = NULL) {
  
  cat("\n🔥 Creating heatmap for", length(gene_list), "genes...\n")
  
  # Create matrix to store mean beta values
  heatmap_matrix <- matrix(NA, nrow = length(gene_list), ncol = nrow(pheno_data))
  rownames(heatmap_matrix) <- gene_list
  colnames(heatmap_matrix) <- rownames(pheno_data)
  
  genes_found <- c()
  genes_missing <- c()
  
  # Calculate mean promoter methylation for each gene
  for (gene in gene_list) {
    # Find gene probes
    gene_indices <- grep(paste0("\\b", gene, "\\b"), annotation$UCSC_RefGene_Name)
    
    if (length(gene_indices) == 0) {
      cat("  ⚠️  No probes found for", gene, "\n")
      genes_missing <- c(genes_missing, gene)
      next
    }
    
    gene_anno <- annotation[gene_indices, ]
    
    # Filter for promoter probes
    is_promoter <- grepl(paste(promoter_regions, collapse="|"), gene_anno$UCSC_RefGene_Group)
    gene_promoter_probes <- gene_anno$Name[is_promoter]
    
    if (length(gene_promoter_probes) == 0) {
      cat("  ⚠️  No promoter probes found for", gene, "\n")
      genes_missing <- c(genes_missing, gene)
      next
    }
    
    # Calculate mean beta values across promoter probes
    # Filter for probes that exist in beta matrix
    available_probes <- gene_promoter_probes[gene_promoter_probes %in% rownames(beta_matrix)]
    
    if (length(available_probes) == 0) {
      cat("  ⚠️  No promoter probes found in beta matrix for", gene, "\n")
      genes_missing <- c(genes_missing, gene)
      next
    }
    
    beta_gene_promoter <- beta_matrix[available_probes, , drop = FALSE]
    gene_promoter_means <- colMeans(beta_gene_promoter, na.rm = TRUE)
    
    # Add to matrix
    heatmap_matrix[gene, ] <- gene_promoter_means
    genes_found <- c(genes_found, gene)
    
    cat("  ✓", gene, "(", length(available_probes), "probes )\n")
  }
  
  # Remove genes with no data
  heatmap_matrix <- heatmap_matrix[genes_found, , drop = FALSE]
  
  cat("\n📊 Successfully found", length(genes_found), "genes\n")
  if (length(genes_missing) > 0) {
    cat("⚠️  Missing data for", length(genes_missing), "genes:", 
        paste(genes_missing, collapse = ", "), "\n")
  }
  
  # Create annotation for samples
  annotation_col <- data.frame(
    Patient = pheno_data$Patient,
    Tissue = pheno_data$Sample_Type,
    row.names = rownames(pheno_data)
  )
  
  # Define colors
  tissue_colors <- c(
    "Gastric"   = "#c85a5a",
    "IM"        = "#8e7cc3",
    "Duodenal"  = "#487cac"
  )
  
  # Get patient colors (using Set3 palette)
  n_patients <- length(unique(pheno_data$Patient))
  patient_colors <- setNames(
    RColorBrewer::brewer.pal(max(3, n_patients), "Set3")[1:n_patients],
    unique(pheno_data$Patient)
  )
  
  annotation_colors <- list(
    Tissue = tissue_colors,
    Patient = patient_colors
  )
  
  # Auto-calculate height if not provided
  if (is.null(height)) {
    height <- max(8, nrow(heatmap_matrix) * 0.25)
  }
  
  # Create heatmap with pheatmap
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  output_path <- file.path(output_dir, output_file)
  
  pheatmap(
    heatmap_matrix,
    cluster_rows = TRUE,           
    cluster_cols = TRUE,          
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = seq(0, 1, length.out = 101),
    main = title,
    fontsize = 10,
    fontsize_row = 8,
    fontsize_col = 8,
    show_colnames = TRUE,
    show_rownames = TRUE,
    border_color = NA,
    filename = output_path,
    width = width,
    height = height
  )
  
  cat("✅ Heatmap saved:", output_path, "\n")
  
  # Also create PNG version
  output_path_png <- sub("\\.pdf$", ".png", output_path)
  pheatmap(
    heatmap_matrix,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = seq(0, 1, length.out = 101),
    main = title,
    fontsize = 10,
    fontsize_row = 8,
    fontsize_col = 8,
    show_colnames = TRUE,
    show_rownames = TRUE,
    border_color = NA,
    filename = output_path_png,
    width = width,
    height = height
  )
  
  cat("✅ PNG version saved:", output_path_png, "\n")
  
  # Return the matrix for further analysis
  return(list(
    matrix = heatmap_matrix,
    genes_found = genes_found,
    genes_missing = genes_missing
  ))
}

############################################
# USAGE EXAMPLES
############################################

# Gene lists
genes_of_interest <- c(
  "ATP4A","BHLHA15","CLDN18","GATA4","GATA6","GIF","GKN1","GKN2","KRT18",
  "MUC5AC","MUC6","PGC","SLC26A7","SOX2","TFF1","TFF2","AGR2","ANXA10","AXIN2",
  "CDH17","CDX2","CLDN10","CLDN7","HNF4A","KRT20","LGR5","MUC1","MUC2","OLFM4",
  "REG4","RNF43","SOX9","SPDEF","TFF3","ALPI","APOA1","APOA4","ASCL2","CDX1",
  "FABP1","FABP2","LYZ","SI","SLC26A3","SLC5A1","VIL1"
)

genes_of_interest_from_RNAseq <- c(
  "A1BG","BGN","BLOC1S5-TXNDC5","CGB5","EPHA2","LIF",
  "LOC102723901", "MAP3K5", "PIK3CG", "SLC2A6","TMX2-CTNND1A",
  "WNT7A", "ZNF625-ZNF20","ZNF664-FAM101A", "HOXA1","ITGB2-AS1",
  "NKX6-3","PCSK1N","ZNF462","CA12","ADAMTSL4","C15orf52","C1orf21","CASC9",
  "COL4A2","FGF19","FLG","GABRP","GRHL3","KRT80","MXRA5","MYEOV","NLRP6","REG1A",
  "SERPINB7","SLC4A3","STAC","TRPM2","ADAMTS14","ALDH1A3","ALDH1A3","AZGP1",
  "CLDN10","GABRE","GCT6","PLXNB3","SHISA6","SIX1","TM6SF2"
)

# Create heatmap 1: Main genes of interest
heatmap1 <- create_sample_heatmap(
  gene_list = genes_of_interest,
  beta_matrix = myNorm,
  annotation = annEPICv2,
  pheno_data = Sample_Data,
  output_file = "Heatmap_GenesOfInterest.pdf",
  title = "Promoter Methylation - Genes of Interest",
  width = 10,
  height = 12
)

# Create heatmap 2: RNA-seq genes
heatmap2 <- create_sample_heatmap(
  gene_list = genes_of_interest_from_RNAseq,
  beta_matrix = myNorm,
  annotation = annEPICv2,
  pheno_data = Sample_Data,
  output_file = "Heatmap_RNAseq_Genes.pdf",
  title = "Promoter Methylation - RNA-seq Genes",
  width = 10,
  height = 14
)



############################################
#Alternative: Differential Methylation Analysis
############################################

library(limma)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

############################################
# STEP 1: Prepare sample data and design
############################################

# Use your Sample_Data (not myLoad$pd)
pd <- Sample_Data

# Ensure factors are set correctly
pd$Sample_Type <- factor(
  pd$Sample_Type,
  levels = c("Gastric", "IM", "Duodenal")
)
pd$Patient <- factor(pd$Patient)

# Verify
print(head(pd))
print(table(pd$Sample_Type))
print(table(pd$Patient))

############################################
# STEP 2: Convert beta to M-values (CORRECTED)
############################################

# Clamp beta values to avoid log(0) or log(Inf)
beta <- pmin(pmax(myNorm, 1e-6), 1 - 1e-6)

# Convert to M-values
Mvals <- log2(beta / (1 - beta))

# Check for any NA or Inf values
sum(is.na(Mvals))
sum(is.infinite(Mvals))

# If you have any NA/Inf, remove those probes
if (sum(is.na(Mvals)) > 0 | sum(is.infinite(Mvals)) > 0) {
  bad_probes <- apply(Mvals, 1, function(x) any(is.na(x) | is.infinite(x)))
  cat("Removing", sum(bad_probes), "probes with NA/Inf values\n")
  Mvals <- Mvals[!bad_probes, ]
}

cat("M-values matrix dimensions:", dim(Mvals), "\n")

############################################
# STEP 3: Get and align annotation
############################################

annEPICv2 <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

# Match probes
match_idx <- match(rownames(Mvals), annEPICv2$Name)
cat("Unmatched probes:", sum(is.na(match_idx)), "\n")

# Subset to matched probes BEFORE modeling
Mvals_sub <- Mvals[!is.na(match_idx), ]
annEPICv2Sub <- annEPICv2[match_idx[!is.na(match_idx)], ]

# Verify alignment
stopifnot(all(rownames(Mvals_sub) == annEPICv2Sub$Name))
cat("Final M-values matrix for modeling:", dim(Mvals_sub), "\n")

############################################
# STEP 4: Create design matrix
############################################

design <- model.matrix(~ Sample_Type, data = pd)
colnames(design)
cat("\nDesign matrix:\n")
print(design)

############################################
# STEP 5: Estimate within-patient correlation
############################################

cat("\nEstimating within-patient correlation...\n")
corfit <- duplicateCorrelation(
  Mvals_sub,
  design,
  block = pd$Patient
)

cat("Consensus correlation:", corfit$consensus, "\n")

# Check if correlation is reasonable
if (corfit$consensus < 0) {
  warning("Negative correlation detected. Patient blocking may not be necessary.")
} else if (corfit$consensus > 0.9) {
  warning("Very high correlation (>0.9). Check for potential issues.")
} else {
  cat("Correlation looks reasonable (0.3-0.8 is typical)\n")
}

############################################
# STEP 6: Fit linear model
############################################

cat("\nFitting linear model...\n")
fit <- lmFit(
  Mvals_sub,
  design,
  block = pd$Patient,
  correlation = corfit$consensus
)

fit <- eBayes(fit)

############################################
# STEP 7: Extract annotated DMPs
############################################

cat("\nExtracting DMPs for all comparisons...\n")

# 1. IM vs Gastric
dmp_IM_vs_Gas <- topTable(
  fit,
  coef = "Sample_TypeIM",
  number = Inf,
  adjust.method = "BH",
  genelist = annEPICv2Sub
)

# 2. Duodenal vs Gastric
dmp_Duo_vs_Gas <- topTable(
  fit,
  coef = "Sample_TypeDuodenal",
  number = Inf,
  adjust.method = "BH",
  genelist = annEPICv2Sub
)

# 3. IM vs Duodenal (via contrast)
contrast.matrix <- makeContrasts(
  IM_vs_Duo = Sample_TypeIM - Sample_TypeDuodenal,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

dmp_IM_vs_Duo <- topTable(
  fit2,
  coef = "IM_vs_Duo",
  number = Inf,
  adjust.method = "BH",
  genelist = annEPICv2Sub
)

############################################
# STEP 8: Verify results
############################################

cat("\n=== VERIFICATION ===\n")

# Check dimensions
cat("\nDimensions:\n")
cat("IM vs Gastric:", dim(dmp_IM_vs_Gas), "\n")
cat("Duodenal vs Gastric:", dim(dmp_Duo_vs_Gas), "\n")
cat("IM vs Duodenal:", dim(dmp_IM_vs_Duo), "\n")

# Check significant DMPs (FDR < 0.05)
cat("\nSignificant DMPs (adj.P < 0.05):\n")
cat("IM vs Gastric:", sum(dmp_IM_vs_Gas$adj.P.Val < 0.05), "\n")
cat("Duodenal vs Gastric:", sum(dmp_Duo_vs_Gas$adj.P.Val < 0.05), "\n")
cat("IM vs Duodenal:", sum(dmp_IM_vs_Duo$adj.P.Val < 0.05), "\n")

# Check significant DMPs with effect size (adj.P < 0.05 & |logFC| > 0.2)
cat("\nSignificant DMPs (adj.P < 0.05 & |logFC| > 0.2):\n")
cat("IM vs Gastric:", sum(dmp_IM_vs_Gas$adj.P.Val < 0.05 & abs(dmp_IM_vs_Gas$logFC) > 0.2), "\n")
cat("Duodenal vs Gastric:", sum(dmp_Duo_vs_Gas$adj.P.Val < 0.05 & abs(dmp_Duo_vs_Gas$logFC) > 0.2), "\n")
cat("IM vs Duodenal:", sum(dmp_IM_vs_Duo$adj.P.Val < 0.05 & abs(dmp_IM_vs_Duo$logFC) > 0.2), "\n")

# Check annotation columns
cat("\nAnnotation columns in dmp_IM_vs_Gas:\n")
print(colnames(dmp_IM_vs_Gas))

# View top results
cat("\nTop 10 DMPs (IM vs Gastric):\n")
print(head(dmp_IM_vs_Gas, 10))

############################################
# STEP 9: Save results
############################################

cat("\nSaving results...\n")
write.csv(dmp_IM_vs_Gas, "DMPs_IM_vs_Gastric.csv", row.names = FALSE)
write.csv(dmp_Duo_vs_Gas, "DMPs_Duodenal_vs_Gastric.csv", row.names = FALSE)
write.csv(dmp_IM_vs_Duo, "DMPs_IM_vs_Duodenal.csv", row.names = FALSE)

cat("\n✅ COMPLETE! All DMP tables saved.\n")


############################################
# STEP 7b: Calculate Beta Differences (Delta Beta)
############################################

dmp_IM_vs_Gas <- read.csv("setdirectory/DMPs_IM_vs_Gastric.csv")
dmp_Duo_vs_Gas <- read.csv("setdirectory/DMPs_Duodenal_vs_Gastric.csv")
dmp_IM_vs_Duodenal <- read.csv("setdirectory/DMPs_IM_vs_Duodenal.csv")


# 1. Get the original Beta values for the matched probes
# (Use the subset we created in Step 3 so rows match the fit object)
Beta_sub <- myNorm[rownames(Mvals_sub), ]

# 2. Calculate Group Means (using Beta values)
# Note: simple rowMeans is usually sufficient for interpretation
mean_Gastric <- rowMeans(Beta_sub[, pd$Sample_Type == "Gastric"], na.rm=TRUE)
mean_IM      <- rowMeans(Beta_sub[, pd$Sample_Type == "IM"], na.rm=TRUE)
mean_Duo     <- rowMeans(Beta_sub[, pd$Sample_Type == "Duodenal"], na.rm=TRUE)

# 3. Add these to your topTables

# --- Process IM vs Gastric ---
dmp_IM_vs_Gas <- topTable(fit, coef="Sample_TypeIM", number=Inf, genelist=annEPICv2Sub, adjust.method="BH")

# Add Beta columns
dmp_IM_vs_Gas$Mean_Gastric <- mean_Gastric[rownames(dmp_IM_vs_Gas)]
dmp_IM_vs_Gas$Mean_IM      <- mean_IM[rownames(dmp_IM_vs_Gas)]
dmp_IM_vs_Gas$Delta_Beta   <- dmp_IM_vs_Gas$Mean_IM - dmp_IM_vs_Gas$Mean_Gastric

# --- Process Duodenal vs Gastric ---
dmp_Duo_vs_Gas <- topTable(fit, coef="Sample_TypeDuodenal", number=Inf, genelist=annEPICv2Sub, adjust.method="BH")

# Add Beta columns
dmp_Duo_vs_Gas$Mean_Gastric <- mean_Gastric[rownames(dmp_Duo_vs_Gas)]
dmp_Duo_vs_Gas$Mean_Duodenal<- mean_Duo[rownames(dmp_Duo_vs_Gas)]
dmp_Duo_vs_Gas$Delta_Beta   <- dmp_Duo_vs_Gas$Mean_Duodenal - dmp_Duo_vs_Gas$Mean_Gastric

# --- Process IM vs Duodenal ---
dmp_IM_vs_Duo <- topTable(fit2, coef="IM_vs_Duo", number=Inf, genelist=annEPICv2Sub, adjust.method="BH")

# Add Beta columns
dmp_IM_vs_Duo$Mean_IM       <- mean_IM[rownames(dmp_IM_vs_Duo)]
dmp_IM_vs_Duo$Mean_Duodenal <- mean_Duo[rownames(dmp_IM_vs_Duo)]
dmp_IM_vs_Duo$Delta_Beta    <- dmp_IM_vs_Duo$Mean_IM - dmp_IM_vs_Duo$Mean_Duodenal

###Visualizations of DMPs
############################################
# A. Prepare volcano plot function
############################################

library(ggplot2)
library(ggrepel)

# Define thresholds
pval_thresh <- 0.05
logfc_thresh <- 1

# Function to create volcano plot
create_volcano <- function(dmp_data, comparison_name, 
                           top_n_genes = 10,
                           adj_pval_thresh = 0.05,  
                           logfc_thresh = 1) {
  
  # Prepare data
  plot_data <- dmp_data %>%
    mutate(
      log10P = -log10(P.Value),
      log10adjP = -log10(adj.P.Val),  
      Significant = case_when(
        adj.P.Val < adj_pval_thresh & logFC > logfc_thresh ~ "Hypermethylated",
        adj.P.Val < adj_pval_thresh & logFC < -logfc_thresh ~ "Hypomethylated",
        TRUE ~ "Not Significant"
      )
    )
  
  # Get top genes to label
  top_hyper <- plot_data %>%
    filter(Significant == "Hypermethylated",
           grepl("TSS", UCSC_RefGene_Group)) %>%
    arrange(adj.P.Val) %>%
    head(top_n_genes)
  
  top_hypo <- plot_data %>%
    filter(Significant == "Hypomethylated",
           grepl("TSS", UCSC_RefGene_Group)) %>%
    arrange(adj.P.Val) %>%
    head(top_n_genes)
  
  genes_to_label <- bind_rows(top_hyper, top_hypo) %>%
    mutate(gene_label = sapply(strsplit(UCSC_RefGene_Name, ";"), `[`, 1)) %>%
    filter(gene_label != "")
  
  # Create plot, using adj.P.Val for y-axis
  p <- ggplot(plot_data, aes(x = logFC, y = log10adjP, color = Significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(
      values = c(
        "Hypermethylated" = "red",
        "Hypomethylated" = "blue",
        "Not Significant" = "gray"
      )
    ) +
    geom_vline(xintercept = c(-logfc_thresh, logfc_thresh), 
               linetype = "dashed", color = "black", alpha = 0.5) +
    geom_hline(yintercept = -log10(adj_pval_thresh),
               linetype = "dashed", color = "black", alpha = 0.5) +
    theme_minimal() +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste("Volcano Plot:", comparison_name),
      x = "log2 Fold Change (M-value difference)",
      y = "-log10(Adjusted P-value)", 
      color = "Methylation Status"
    ) +
    geom_text_repel(
      data = genes_to_label,
      aes(label = gene_label),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "gray50",
      segment.size = 0.3
    )
  
  # Add count annotations
  n_hyper <- sum(plot_data$Significant == "Hypermethylated")
  n_hypo <- sum(plot_data$Significant == "Hypomethylated")
  
  p <- p + annotate("text", 
                    x = max(plot_data$logFC) * 0.7, 
                    y = max(plot_data$log10adjP) * 0.95,
                    label = paste0("Hypermethylated: ", n_hyper),
                    color = "red", size = 4, fontface = "bold") +
    annotate("text",
             x = min(plot_data$logFC) * 0.7,
             y = max(plot_data$log10adjP) * 0.95,
             label = paste0("Hypomethylated: ", n_hypo),
             color = "blue", size = 4, fontface = "bold")
  
  return(p)
}
############################################
# B. Create all three volcano plots
############################################

# IM vs Gastric
volcano_IM_vs_Gas <- create_volcano(
  dmp_IM_vs_Gas, 
  "IM vs Gastric",
  top_n_genes = 10
)
print(volcano_IM_vs_Gas)
ggsave("Volcano_IM_vs_Gastric.png", volcano_IM_vs_Gas, 
       width = 10, height = 8, dpi = 300)

# Duodenal vs Gastric
volcano_Duo_vs_Gas <- create_volcano(
  dmp_Duo_vs_Gas,
  "Duodenal vs Gastric",
  top_n_genes = 10
)
print(volcano_Duo_vs_Gas)
ggsave("Volcano_Duodenal_vs_Gastric.png", volcano_Duo_vs_Gas,
       width = 10, height = 8, dpi = 300)

# IM vs Duodenal
volcano_IM_vs_Duo <- create_volcano(
  dmp_IM_vs_Duo,
  "IM vs Duodenal",
  top_n_genes = 10
)
print(volcano_IM_vs_Duo)
ggsave("Volcano_IM_vs_Duodenal.png", volcano_IM_vs_Duo,
       width = 10, height = 8, dpi = 300)

############################################
# C. Combine all three volcano plots
############################################

library(patchwork)

# Create combined plot
combined_volcano <- volcano_IM_vs_Gas + volcano_Duo_vs_Gas + volcano_IM_vs_Duo +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

print(combined_volcano)
ggsave("Volcano_Combined_All_Comparisons.png", combined_volcano,
       width = 18, height = 6, dpi = 300)

# Alternative: stacked layout
combined_volcano_stacked <- volcano_IM_vs_Gas / volcano_Duo_vs_Gas / volcano_IM_vs_Duo +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(combined_volcano_stacked)
ggsave("Volcano_Combined_Stacked.png", combined_volcano_stacked,
       width = 10, height = 18, dpi = 300)


####################################################################
#Genomic Distrubution Bar Charts, coloured by hyper/hypomethylation
####################################################################
library(tidyverse)
library(ggplot2)
library(patchwork)

# Define the logical biological order for the X-axis
region_order <- c("TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "3'UTR", "IGR")

# ---------------------------------------------------------
# HELPER FUNCTION: Maps raw Illumina strings to Clean Groups
# ---------------------------------------------------------
clean_annotation_string <- function(raw_string) {
  # 1. Handle Empty/NA -> IGR
  if (is.na(raw_string) || raw_string == "") return("IGR")
  
  # 2. Split by semicolon and take the first entry
  # (Illumina often lists multiples like "exon_2;exon_2". We just need the first one.)
  first_item <- strsplit(raw_string, ";")[[1]][1]
  
  # 3. Categorize based on your specific file format
  if (first_item == "TSS1500") return("TSS1500")
  if (first_item == "TSS200")  return("TSS200")
  if (first_item == "5UTR")    return("5'UTR")
  if (first_item == "3UTR")    return("3'UTR")  
  if (first_item == "exon_1")  return("1stExon")
  
  # 4. Catch-all for other exons (exon_2, exon_3, exon_115...) -> "Body"
  if (grepl("exon_", first_item)) return("Body")
  
  # 5. Fallback for anything else
  return("IGR") 
}

# ---------------------------------------------------------
# PLOTTING FUNCTION
# ---------------------------------------------------------
plot_clean_proportions <- function(dmp_data, comparison_name) {
  
  # 1. Clean the Data
  # We apply the cleaning function to every row
  plot_data <- dmp_data %>%
    filter(adj.P.Val < 0.05) %>%
    rowwise() %>%
    mutate(
      Clean_Region = clean_annotation_string(as.character(UCSC_RefGene_Group)),
      Direction = ifelse(logFC > 0, "Hypermethylated", "Hypomethylated")
    ) %>%
    ungroup() %>%
    filter(Clean_Region %in% region_order) 
  
  # 2. Calculate Percentages within each Direction (Hyper vs Hypo)
  summary_data <- plot_data %>%
    group_by(Direction, Clean_Region) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    group_by(Direction) %>%
    mutate(Percentage = (Count / sum(Count)) * 100)
  
  # 3. Enforce X-axis Order
  summary_data$Clean_Region <- factor(summary_data$Clean_Region, levels = region_order)
  
  # 4. Generate Plot
  p <- ggplot(summary_data, aes(x = Clean_Region, y = Percentage, fill = Direction)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    scale_fill_manual(values = c("Hypermethylated" = "#D53E4F", 
                                 "Hypomethylated" = "#3288BD")) +
    labs(
      title = comparison_name,
      y = "% of DMPs",
      x = NULL,
      fill = ""
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11)
    ) +
    # Add percentage labels above bars
    geom_text(aes(label = sprintf("%.1f", Percentage)), 
              position = position_dodge(width = 0.7), 
              vjust = -0.5, size = 2.5) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) 
  
  return(p)
}

# ---------------------------------------------------------
# RUN AND SAVE
# ---------------------------------------------------------

# Create the 3 plots
p1 <- plot_clean_proportions(dmp_IM_vs_Gas, "IM vs Gastric")
p2 <- plot_clean_proportions(dmp_Duo_vs_Gas, "Duodenal vs Gastric")
p3 <- plot_clean_proportions(dmp_IM_vs_Duo, "IM vs Duodenal")

# Combine them side-by-side
combined_plot <- p1 + p2 + p3 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

print(combined_plot)

# Save high-res version
ggsave("Genomic_Region_Distribution_Fixed.png", combined_plot, width = 12, height = 5, dpi = 300)


######################################
#Do the same, but for CGI Context.
######################################

library(tidyverse)
library(ggplot2)
library(patchwork)

# Define CpG context order (Illumina standard)
cpg_order <- c("Island", "Shore", "Shelf", "OpenSea")

# Plotting Function
plot_cpg_proportions <- function(dmp_data, comparison_name) {
  
  plot_data <- dmp_data %>%
    filter(adj.P.Val < 0.05) %>%
    mutate(
      CpG_Context = as.character(Relation_to_Island),
      Direction = ifelse(logFC > 0,
                         "Hypermethylated",
                         "Hypomethylated")
    ) %>%
    filter(CpG_Context %in% cpg_order)
  
  # Calculate percentages within hyper / hypo
  summary_data <- plot_data %>%
    group_by(Direction, CpG_Context) %>%
    summarise(Count = n(), .groups = "drop") %>%
    group_by(Direction) %>%
    mutate(Percentage = (Count / sum(Count)) * 100)
  
  # enforce order
  summary_data$CpG_Context <- factor(
    summary_data$CpG_Context,
    levels = cpg_order
  )
  
  # -------------------------------------------------------
  # Plot
  # -------------------------------------------------------
  p <- ggplot(
    summary_data,
    aes(x = CpG_Context, y = Percentage, fill = Direction)
  ) +
    geom_bar(
      stat = "identity",
      position = "dodge",
      width = 0.7
    ) +
    scale_fill_manual(
      values = c(
        "Hypermethylated" = "#D53E4F",
        "Hypomethylated"  = "#3288BD"
      )
    ) +
    labs(
      title = comparison_name,
      y = "% of DMPs",
      x = NULL,
      fill = ""
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        color = "black"
      ),
      legend.position = "top",
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 11
      )
    ) +
    geom_text(
      aes(label = sprintf("%.1f", Percentage)),
      position = position_dodge(width = 0.7),
      vjust = -0.4,
      size = 2.7
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.15))
    )
  
  return(p)
}

p1 <- plot_cpg_proportions(dmp_IM_vs_Gas, "IM vs Gastric")
p2 <- plot_cpg_proportions(dmp_Duo_vs_Gas, "Duodenal vs Gastric")
p3 <- plot_cpg_proportions(dmp_IM_vs_Duo, "IM vs Duodenal")

combined_plot <- p1 + p2 + p3 +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(combined_plot)

ggsave(
  "CpG_Context_Distribution.png",
  combined_plot,
  width = 12,
  height = 5,
  dpi = 300
)





############################################
# Genome-wide Manhattan Plot
############################################

library(ggplot2)
library(dplyr)
library(ggrepel)

create_manhattan_plot <- function(dmp_data, comparison_name,
                                  fdr_thresh = 0.05,
                                  bonf_thresh = NULL,
                                  label_top_n = 20) {
  
  # Calculate Bonferroni threshold if not provided
  n_tests <- nrow(dmp_data)
  if (is.null(bonf_thresh)) {
    bonf_thresh <- 0.05 / n_tests
  }
  
  # Prepare data
  manhattan_data <- dmp_data %>%
    filter(!is.na(chr), !is.na(pos)) %>%
    mutate(
      # Clean chromosome names and convert to numeric
      chr_clean = gsub("chr", "", chr),
      chr_num = case_when(
        chr_clean == "X" ~ 23,
        chr_clean == "Y" ~ 24,
        chr_clean == "M" ~ 25,
        TRUE ~ as.numeric(chr_clean)
      ),
      # Calculate -log10 P-values
      log10P = -log10(P.Value),
      # Determine significance status
      Significance = case_when(
        P.Value < bonf_thresh ~ "Bonferroni",
        adj.P.Val < fdr_thresh ~ "FDR",
        TRUE ~ "Not Significant"
      )
    ) %>%
    filter(!is.na(chr_num)) %>%
    arrange(chr_num, pos)
  
  # Calculate cumulative positions for x-axis
  # Use as.numeric() to avoid integer overflow
  chr_lengths <- manhattan_data %>%
    group_by(chr_num) %>%
    summarise(
      chr_len = as.numeric(max(pos)),
      n_probes = n(),
      .groups = "drop"
    ) %>%
    arrange(chr_num) %>%
    mutate(
      # Use as.numeric to prevent overflow
      cumsum_len = cumsum(chr_len),
      start_pos = lag(cumsum_len, default = 0),
      mid_pos = start_pos + chr_len / 2
    )
  
  # Add cumulative positions to data
  # Explicitly use dplyr::select to avoid any AnnotationDbi conflict
  manhattan_data <- manhattan_data %>%
    left_join(
      chr_lengths %>% dplyr::select(chr_num, start_pos), 
      by = "chr_num"
    ) %>%
    mutate(pos_cumulative = start_pos + as.numeric(pos))
  
  # Identify probes to label (Bonferroni significant with gene names)
  probes_to_label <- manhattan_data %>%
    filter(
      Significance == "Bonferroni",
      grepl("TSS", UCSC_RefGene_Group),
      UCSC_RefGene_Name != ""
    ) %>%
    mutate(gene_label = sapply(strsplit(UCSC_RefGene_Name, ";"), `[`, 1)) %>%
    filter(gene_label != "") %>%
    arrange(P.Value) %>%
    head(label_top_n)
  
  # If no Bonferroni hits, use FDR instead
  if (nrow(probes_to_label) == 0) {
    cat("Note: No Bonferroni-significant probes found. Using FDR threshold for labeling.\n")
    probes_to_label <- manhattan_data %>%
      filter(
        Significance == "FDR",
        grepl("TSS", UCSC_RefGene_Group),
        UCSC_RefGene_Name != ""
      ) %>%
      mutate(gene_label = sapply(strsplit(UCSC_RefGene_Name, ";"), `[`, 1)) %>%
      filter(gene_label != "") %>%
      arrange(adj.P.Val) %>%
      head(label_top_n)
  }
  
  # Calculate threshold lines
  bonf_line <- -log10(bonf_thresh)
  # For FDR line, find the minimum P-value among FDR-significant probes
  fdr_pvals <- manhattan_data %>%
    filter(adj.P.Val < fdr_thresh) %>%
    pull(P.Value)
  
  if (length(fdr_pvals) > 0) {
    fdr_line <- -log10(max(fdr_pvals))
  } else {
    fdr_line <- -log10(0.05
  }
  
  # Create the plot
  p <- ggplot(manhattan_data, aes(x = pos_cumulative, y = log10P)) +
    # Plot points - colored by chromosome for alternating pattern
    geom_point(
      aes(color = factor(chr_num %% 2)),
      data = manhattan_data %>% filter(Significance == "Not Significant"),
      alpha = 0.5, size = 0.8
    ) +
    # FDR significant points (orange/yellow)
    geom_point(
      data = manhattan_data %>% filter(Significance == "FDR"),
      color = "orange", alpha = 0.7, size = 1.2
    ) +
    # Bonferroni significant points (red)
    geom_point(
      data = manhattan_data %>% filter(Significance == "Bonferroni"),
      color = "red", alpha = 0.8, size = 1.5
    ) +
    # Color scheme for alternating chromosomes
    scale_color_manual(
      values = c("0" = "gray30", "1" = "gray60"),
      guide = "none"
    ) +
    # Threshold lines
    geom_hline(
      yintercept = bonf_line,
      linetype = "dashed",
      color = "red",
      linewidth = 0.8
    ) +
    geom_hline(
      yintercept = fdr_line,
      linetype = "dashed",
      color = "orange",
      linewidth = 0.8
    ) +
    # Add chromosome labels on x-axis
    scale_x_continuous(
      breaks = chr_lengths$mid_pos,
      labels = c(1:22, "X", "Y")[1:nrow(chr_lengths)],
      expand = c(0.01, 0.01)
    ) +
    # Labels for significant genes
    {if(nrow(probes_to_label) > 0) {
      geom_text_repel(
        data = probes_to_label,
        aes(label = gene_label),
        size = 3,
        max.overlaps = 30,
        box.padding = 0.3,
        point.padding = 0.2,
        segment.color = "gray50",
        segment.size = 0.3,
        min.segment.length = 0,
        force = 2
      )
    }} +
    # Theme
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
      axis.ticks.x = element_line(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    labs(
      title = paste("Manhattan Plot:", comparison_name),
      x = "Chromosome",
      y = expression(-log[10](P-value))
    )
  
  # Add threshold annotations conditionally
  max_y <- max(manhattan_data$log10P, na.rm = TRUE)
  max_x <- max(manhattan_data$pos_cumulative, na.rm = TRUE)
  
  if (bonf_line < max_y) {
    p <- p + annotate(
      "text",
      x = max_x * 0.85,
      y = bonf_line + max_y * 0.05,
      label = paste0("Bonferroni (p < ", formatC(bonf_thresh, format = "e", digits = 1), ")"),
      color = "red",
      size = 3.5,
      fontface = "bold"
    )
  }
  
  if (fdr_line < max_y) {
    p <- p + annotate(
      "text",
      x = max_x * 0.85,
      y = fdr_line + max_y * 0.05,
      label = paste0("FDR (q < ", fdr_thresh, ")"),
      color = "orange",
      size = 3.5,
      fontface = "bold"
    )
  }
  
  # Print summary statistics
  n_bonf <- sum(manhattan_data$Significance == "Bonferroni")
  n_fdr <- sum(manhattan_data$Significance == "FDR")
  
  cat("\n=== Manhattan Plot Summary ===\n")
  cat("Comparison:", comparison_name, "\n")
  cat("Total probes tested:", n_tests, "\n")
  cat("Bonferroni threshold:", formatC(bonf_thresh, format = "e", digits = 3), "\n")
  cat("FDR threshold:", fdr_thresh, "\n")
  cat("Bonferroni significant probes:", n_bonf, "\n")
  cat("FDR significant probes:", n_fdr, "\n")
  cat("Labeled genes:", nrow(probes_to_label), "\n\n")
  
  return(list(plot = p, data = manhattan_data, labeled = probes_to_label))
}


############################################
# Create Manhattan plots for all comparisons
############################################

# IM vs Gastric
manhattan_IM_Gas <- create_manhattan_plot(
  dmp_IM_vs_Gas,
  "IM vs Gastric",
  fdr_thresh = 0.05,
  label_top_n = 20
)
print(manhattan_IM_Gas$plot)
ggsave("Manhattan_IM_vs_Gastric.png", manhattan_IM_Gas$plot,
       width = 14, height = 6, dpi = 300)

# Duodenal vs Gastric
manhattan_Duo_Gas <- create_manhattan_plot(
  dmp_Duo_vs_Gas,
  "Duodenal vs Gastric",
  fdr_thresh = 0.05,
  label_top_n = 20
)
print(manhattan_Duo_Gas$plot)
ggsave("Manhattan_Duodenal_vs_Gastric.png", manhattan_Duo_Gas$plot,
       width = 14, height = 6, dpi = 300)

# IM vs Duodenal
manhattan_IM_Duo <- create_manhattan_plot(
  dmp_IM_vs_Duo,
  "IM vs Duodenal",
  fdr_thresh = 0.05,
  label_top_n = 20
)
print(manhattan_IM_Duo$plot)
ggsave("Manhattan_IM_vs_Duodenal.png", manhattan_IM_Duo$plot,
       width = 14, height = 6, dpi = 300)


##################################################
#Venn Diagrams for looking at what DMPs are shared
##################################################

library(VennDiagram)
library(dplyr)
library(openxlsx)
library(gridExtra)
library(grid)

# Set significance thresholds
pval_thresh <- 0.05
logfc_thresh <- 0.2

############################################
# STEP 1: Prepare DMP sets for each comparison
############################################

# All significant DMPs
dmps_IM_Gas_all <- dmp_IM_vs_Gas %>%
  filter(adj.P.Val < pval_thresh, abs(logFC) > logfc_thresh) %>%
  mutate(comparison = "IM_vs_Gastric")

dmps_Duo_Gas_all <- dmp_Duo_vs_Gas %>%
  filter(adj.P.Val < pval_thresh, abs(logFC) > logfc_thresh) %>%
  mutate(comparison = "Duodenal_vs_Gastric")

dmps_IM_Duo_all <- dmp_IM_vs_Duo %>%
  filter(adj.P.Val < pval_thresh, abs(logFC) > logfc_thresh) %>%
  mutate(comparison = "IM_vs_Duodenal")

# Promoter-only DMPs
dmps_IM_Gas_promoter <- dmps_IM_Gas_all %>%
  filter(grepl("TSS", UCSC_RefGene_Group))

dmps_Duo_Gas_promoter <- dmps_Duo_Gas_all %>%
  filter(grepl("TSS", UCSC_RefGene_Group))

dmps_IM_Duo_promoter <- dmps_IM_Duo_all %>%
  filter(grepl("TSS", UCSC_RefGene_Group))

# Hypermethylated DMPs
dmps_IM_Gas_hyper <- dmps_IM_Gas_all %>%
  filter(logFC > logfc_thresh)

dmps_Duo_Gas_hyper <- dmps_Duo_Gas_all %>%
  filter(logFC > logfc_thresh)

dmps_IM_Duo_hyper <- dmps_IM_Duo_all %>%
  filter(logFC > logfc_thresh)

# Hypomethylated DMPs
dmps_IM_Gas_hypo <- dmps_IM_Gas_all %>%
  filter(logFC < -logfc_thresh)

dmps_Duo_Gas_hypo <- dmps_Duo_Gas_all %>%
  filter(logFC < -logfc_thresh)

dmps_IM_Duo_hypo <- dmps_IM_Duo_all %>%
  filter(logFC < -logfc_thresh)

############################################
# STEP 2: Function to create Venn diagrams
############################################

create_clean_venn <- function(set1_data, set2_data, set3_data,
                              set1_name, set2_name, set3_name,
                              title, filename_prefix,
                              fill_colors = c("#c85a5a", "#8e7cc3", "#487cac")) {
  
  # Get probe name vectors
  set1_probes <- set1_data$Name
  set2_probes <- set2_data$Name
  set3_probes <- set3_data$Name
  
  # Create clean Venn diagram with pyramid layout
  venn_plot <- venn.diagram(
    x = list(
      set1_probes,
      set2_probes,
      set3_probes
    ),
    category.names = c(set1_name, set2_name, set3_name),
    filename = NULL,
    
    # Circle styling
    fill = fill_colors,
    alpha = 0.5,
    lwd = 2,
    col = "black",
    
    # Number styling
    cex = 2.5,
    fontfamily = "sans",
    fontface = "bold",
    
    # Category label styling
    cat.cex = 2,
    cat.fontface = "bold",
    cat.fontfamily = "sans",
    cat.default.pos = "outer",
    
    # Positioning for pyramid layout
    cat.pos = c(-27, 27, 180),  
    cat.dist = c(0.055, 0.055, 0.085),  
    
    # Title
    main = title,
    main.cex = 2.2,
    main.fontface = "bold",
    main.fontfamily = "sans",
    main.pos = c(0.5, 1.05),
    
    # Remove extra elements
    print.mode = "raw",
    sigdigs = 3,
    
    # Better spacing
    margin = 0.15
  )
  
  # Save individual Venn diagram
  png(paste0(filename_prefix, "_VennDiagram.png"),
      width = 1200, height = 1200, res = 150)
  grid.newpage()
  grid.draw(venn_plot)
  dev.off()
  
  # Calculate intersections
  unique_set1 <- setdiff(setdiff(set1_probes, set2_probes), set3_probes)
  unique_set2 <- setdiff(setdiff(set2_probes, set1_probes), set3_probes)
  unique_set3 <- setdiff(setdiff(set3_probes, set1_probes), set2_probes)
  
  set1_set2_only <- setdiff(intersect(set1_probes, set2_probes), set3_probes)
  set1_set3_only <- setdiff(intersect(set1_probes, set3_probes), set2_probes)
  set2_set3_only <- setdiff(intersect(set2_probes, set3_probes), set1_probes)
  
  all_three <- Reduce(intersect, list(set1_probes, set2_probes, set3_probes))
  
  # Print summary
  cat("\n========================================\n")
  cat(title, "\n")
  cat("========================================\n")
  cat(sprintf("%-40s: %d\n", paste("Unique to", set1_name), length(unique_set1)))
  cat(sprintf("%-40s: %d\n", paste("Unique to", set2_name), length(unique_set2)))
  cat(sprintf("%-40s: %d\n", paste("Unique to", set3_name), length(unique_set3)))
  cat(sprintf("%-40s: %d\n", paste(set1_name, "&", set2_name, "only"), length(set1_set2_only)))
  cat(sprintf("%-40s: %d\n", paste(set1_name, "&", set3_name, "only"), length(set1_set3_only)))
  cat(sprintf("%-40s: %d\n", paste(set2_name, "&", set3_name, "only"), length(set2_set3_only)))
  cat(sprintf("%-40s: %d\n", "All three comparisons", length(all_three)))
  cat(sprintf(
    "%-40s: %d\n",
    "TOTAL",
    length(unique(c(set1_probes, set2_probes, set3_probes)))
  ))
  cat("\n")
  
  # Combine all data
  all_data <- bind_rows(set1_data, set2_data, set3_data)
  
  get_annotated_probes <- function(probe_list) {
    if (length(probe_list) == 0) {
      return(data.frame(Message = "No probes in this intersection"))
    }
    
    all_data %>%
      filter(Name %in% probe_list) %>%
      dplyr::select(
        Name, chr, pos, logFC, P.Value, adj.P.Val,
        UCSC_RefGene_Name, UCSC_RefGene_Group,
        Relation_to_Island, comparison
      ) %>%
      arrange(adj.P.Val)
  }
  
  # Create Excel file
  wb <- createWorkbook()
  
  # Summary sheet
  summary_df <- data.frame(
    Intersection = c(
      paste("Unique to", set1_name),
      paste("Unique to", set2_name),
      paste("Unique to", set3_name),
      paste(set1_name, "&", set2_name, "only"),
      paste(set1_name, "&", set3_name, "only"),
      paste(set2_name, "&", set3_name, "only"),
      "All three comparisons",
      "TOTAL unique probes"
    ),
    Count = c(
      length(unique_set1),
      length(unique_set2),
      length(unique_set3),
      length(set1_set2_only),
      length(set1_set3_only),
      length(set2_set3_only),
      length(all_three),
      length(unique(c(set1_probes, set2_probes, set3_probes)))
    )
  )
  
  addWorksheet(wb, "Summary")
  writeData(wb, "Summary", summary_df)
  
  # Add detailed sheets for each intersection
  # Use shorter names to avoid Excel's 31-character limit
  addWorksheet(wb, "Unique_Set1")
  writeData(wb, "Unique_Set1", get_annotated_probes(unique_set1))
  
  addWorksheet(wb, "Unique_Set2")
  writeData(wb, "Unique_Set2", get_annotated_probes(unique_set2))
  
  addWorksheet(wb, "Unique_Set3")
  writeData(wb, "Unique_Set3", get_annotated_probes(unique_set3))
  
  addWorksheet(wb, "Set1_and_Set2_only")
  writeData(wb, "Set1_and_Set2_only", get_annotated_probes(set1_set2_only))
  
  addWorksheet(wb, "Set1_and_Set3_only")
  writeData(wb, "Set1_and_Set3_only", get_annotated_probes(set1_set3_only))
  
  addWorksheet(wb, "Set2_and_Set3_only")
  writeData(wb, "Set2_and_Set3_only", get_annotated_probes(set2_set3_only))
  
  addWorksheet(wb, "All_Three")
  writeData(wb, "All_Three", get_annotated_probes(all_three))
  
  # Add a key sheet explaining which set is which
  key_df <- data.frame(
    Sheet = c("Set1", "Set2", "Set3"),
    Comparison = c(set1_name, set2_name, set3_name)
  )
  addWorksheet(wb, "Sheet_Key")
  writeData(wb, "Sheet_Key", key_df)
  
  saveWorkbook(
    wb,
    paste0(filename_prefix, "_Intersections.xlsx"),
    overwrite = TRUE
  )
  
  cat(sprintf("✅ Saved: %s_VennDiagram.png\n", filename_prefix))
  cat(sprintf("✅ Saved: %s_Intersections.xlsx\n\n", filename_prefix))
  
  return(venn_plot)
}

############################################
# STEP 3: Define consistent color scheme
############################################

# Use your tissue colors consistently across all Venns
tissue_colors <- c("#c85a5a", "#8e7cc3", "#487cac")  # Gastric, IM, Duodenal

# For hyper/hypo, use color gradients
colors_all <- c("#c85a5a", "#8e7cc3", "#487cac")
colors_promoter <- c("#b04848", "#7a6aaf", "#3a6a94")
colors_hyper <- c("#d46a6a", "#9e8cd3", "#588cbc")  # Lighter shades
colors_hypo <- c("#a84a4a", "#7e6cb3", "#386c9c")  # Darker shades

############################################
# STEP 4: Create all 4 Venn diagrams
############################################

# 1. ALL DMPs
venn1 <- create_clean_venn(
  dmps_IM_Gas_all, dmps_Duo_Gas_all, dmps_IM_Duo_all,
  "IM vs Gastric", "Duodenal vs Gastric", "IM vs Duodenal",
  "Overlap of All Significant DMPs",
  "Venn1_All_DMPs",
  fill_colors = colors_all
)

# 2. PROMOTER DMPs only
venn2 <- create_clean_venn(
  dmps_IM_Gas_promoter, dmps_Duo_Gas_promoter, dmps_IM_Duo_promoter,
  "IM vs Gastric", "Duodenal vs Gastric", "IM vs Duodenal",
  "Overlap of Promoter DMPs (TSS200/TSS1500)",
  "Venn2_Promoter_DMPs",
  fill_colors = colors_promoter
)

# 3. HYPERMETHYLATED DMPs only
venn3 <- create_clean_venn(
  dmps_IM_Gas_hyper, dmps_Duo_Gas_hyper, dmps_IM_Duo_hyper,
  "IM vs Gastric", "Duodenal vs Gastric", "IM vs Duodenal",
  "Overlap of Hypermethylated DMPs",
  "Venn3_Hypermethylated_DMPs",
  fill_colors = colors_hyper
)

# 4. HYPOMETHYLATED DMPs only
venn4 <- create_clean_venn(
  dmps_IM_Gas_hypo, dmps_Duo_Gas_hypo, dmps_IM_Duo_hypo,
  "IM vs Gastric", "Duodenal vs Gastric", "IM vs Duodenal",
  "Overlap of Hypomethylated DMPs",
  "Venn4_Hypomethylated_DMPs",
  fill_colors = colors_hypo
)

############################################
# STEP 5: Create combined panel figure
############################################

png("Combined_Venn_Diagrams_All_Four.png", 
    width = 2400, height = 2400, res = 150)

# Create layout with proper spacing
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 2, 
                                           heights = unit(c(0.1, 1, 1), "null"),
                                           widths = unit(c(1, 1), "null"))))

# Main title
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:2))
grid.text("DMP Overlap Analysis Across Comparisons", 
          gp = gpar(fontsize = 28, fontface = "bold"))
popViewport()

# Plot 1: All DMPs (top left)
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
grid.draw(venn1)
popViewport()

# Plot 2: Promoter DMPs (top right)
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
grid.draw(venn2)
popViewport()

# Plot 3: Hypermethylated (bottom left)
pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
grid.draw(venn3)
popViewport()

# Plot 4: Hypomethylated (bottom right)
pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 2))
grid.draw(venn4)
popViewport()

dev.off()

############################################
# BONUS: Create a color legend file
############################################

png("Venn_Color_Legend.png", width = 600, height = 400, res = 150)
par(mar = c(2, 2, 3, 2))
plot.new()
title("Venn Diagram Color Scheme", cex.main = 1.5, font.main = 2)

legend("center", 
       legend = c("IM vs Gastric", 
                  "Duodenal vs Gastric", 
                  "IM vs Duodenal"),
       fill = tissue_colors,
       border = "black",
       bty = "n",
       cex = 1.5,
       pt.cex = 3,
       title = "Comparison Colors",
       title.adj = 0.5)

text(0.5, 0.2, 
     "These colors are used consistently\nacross all Venn diagrams", 
     cex = 1.2, col = "gray40")
dev.off()

cat("\n========================================\n")
cat("✅ COMPLETE! Generated:\n")
cat("========================================\n")
cat("  📊 4 individual Venn diagrams (PNG) - clean pyramid layout\n")
cat("  📁 4 Excel files with detailed probe lists\n")
cat("  🖼️  1 combined panel figure (2x2 grid)\n")
cat("  🎨 1 color legend reference\n\n")
cat("All Venn diagrams use:\n")
cat("  ✓ Consistent color scheme matching your tissue colors\n")
cat("  ✓ Pyramid arrangement for visual clarity\n")
cat("  ✓ Clear, non-overlapping labels\n")
cat("  ✓ Bold, readable text\n\n")


####################################
#Exploring DMP Significance & Burden
####################################
plot_top_promoter_genes <- function(
    dmp_data,
    comparison_name,
    top_n = 50,
    fdr = 0.05
) {
  
  promoter_terms <- c("TSS200", "TSS1500")
  
  promoter_dmps <- dmp_data %>%
    dplyr::filter(
      adj.P.Val < fdr,
      grepl(paste(promoter_terms, collapse = "|"),
            UCSC_RefGene_Group),
      UCSC_RefGene_Name != ""
    ) %>%
    dplyr::mutate(
      Gene = sapply(strsplit(UCSC_RefGene_Name, ";"), `[`, 1),
      Direction = ifelse(logFC > 0,
                         "Hypermethylated",
                         "Hypomethylated")
    )
  
  gene_summary <- promoter_dmps %>%
    dplyr::count(Gene, Direction) %>%
    tidyr::pivot_wider(
      names_from = Direction,
      values_from = n,
      values_fill = 0
    ) %>%
    dplyr::mutate(
      Total = Hypermethylated + Hypomethylated
    ) %>%
    dplyr::arrange(desc(Total)) %>%
    dplyr::slice_head(n = top_n)
  
  plot_df <- gene_summary %>%
    dplyr::select(Gene, Hypermethylated, Hypomethylated) %>%
    tidyr::pivot_longer(
      cols = c(Hypermethylated, Hypomethylated),
      names_to = "Direction",
      values_to = "Count"
    )
  
  plot_df$Gene <- factor(
    plot_df$Gene,
    levels = gene_summary$Gene
  )
  
  ggplot(
    plot_df,
    aes(x = Gene, y = Count, fill = Direction)
  ) +
    geom_bar(stat = "identity", width = 0.8) +
    coord_flip() +
    scale_fill_manual(
      values = c(
        "Hypermethylated" = "#D53E4F",
        "Hypomethylated"  = "#3288BD"
      )
    ) +
    labs(
      title = paste0(
        comparison_name,
        "\nTop ", top_n,
        " genes by promoter DMP burden"
      ),
      x = NULL,
      y = "Number of promoter DMPs",
      fill = ""
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      axis.text.y = element_text(size = 8),
      legend.position = "top"
    )
}

p_IM_Gas <- plot_top_promoter_genes(
  dmp_IM_vs_Gas,
  "IM vs Gastric",
  top_n = 50
)

p_Duo_Gas <- plot_top_promoter_genes(
  dmp_Duo_vs_Gas,
  "Duodenal vs Gastric",
  top_n = 50
)

p_IM_Duo <- plot_top_promoter_genes(
  dmp_IM_vs_Duo,
  "IM vs Duodenal",
  top_n = 50
)

ggsave(
  "Top50_Promoter_Genes_IM_vs_Gastric.png",
  p_IM_Gas,
  width = 8,
  height = 12,
  dpi = 300
)

ggsave(
  "Top50_Promoter_Genes_Duodenal_vs_Gastric.png",
  p_Duo_Gas,
  width = 8,
  height = 12,
  dpi = 300
)

ggsave(
  "Top50_Promoter_Genes_IM_vs_Duodenal.png",
  p_IM_Duo,
  width = 8,
  height = 12,
  dpi = 300
)


############################################
# Analyze and plot gene promoter methylation
############################################

analyze_gene_promoter <- function(gene_name, 
                                  beta_matrix = myNorm,
                                  annotation = annEPICv2Sub,
                                  pheno_data = pd,
                                  promoter_regions = c("TSS1500", "TSS200"),
                                  output_dir = "Gene_Methylation_Plots",
                                  save_plot = TRUE,
                                  show_stats = TRUE) {
  
  # Create output directory if it doesn't exist
  if (save_plot && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("\n========================================\n")
  cat("Analyzing:", gene_name, "\n")
  cat("========================================\n")
  
  # 1. Find gene probes
  gene_indices <- grep(paste0("\\b", gene_name, "\\b"), annotation$UCSC_RefGene_Name)
  
  if (length(gene_indices) == 0) {
    cat("⚠️  No probes found for", gene_name, "\n")
    return(NULL)
  }
  
  gene_anno <- annotation[gene_indices, ]
  
  # 2. Filter for PROMOTER probes
  is_promoter <- grepl(paste(promoter_regions, collapse="|"), gene_anno$UCSC_RefGene_Group)
  gene_promoter_probes <- gene_anno$Name[is_promoter]
  
  if (length(gene_promoter_probes) == 0) {
    cat("⚠️  No promoter probes found for", gene_name, "\n")
    return(NULL)
  }
  
  cat("Found", length(gene_promoter_probes), "promoter probes for", gene_name, "\n")
  
  # 3. Extract Beta values and calculate mean
  beta_gene_promoter <- beta_matrix[gene_promoter_probes, , drop = FALSE]
  gene_promoter_means <- colMeans(beta_gene_promoter, na.rm = TRUE)
  
  # 4. Add to phenotype data (create temporary copy to avoid modifying original)
  plot_data <- pheno_data
  plot_data$Gene_Promoter_Mean <- gene_promoter_means
  
  # 5. Statistical testing using mixed-effects model
  cat("\nRunning mixed-effects model...\n")
  model <- lmer(Gene_Promoter_Mean ~ Sample_Type + (1|Patient), data = plot_data)
  
  # Print model summary
  if (show_stats) {
    cat("\n--- Model Summary ---\n")
    print(summary(model))
  }
  
  # 6. Pairwise comparisons with emmeans
  cat("\nPerforming pairwise comparisons...\n")
  emm <- emmeans(model, pairwise ~ Sample_Type, adjust = "tukey")
  comparisons <- as.data.frame(emm$contrasts)
  
  cat("\n--- Pairwise Comparisons ---\n")
  print(comparisons)
  
  # 7. Format p-values for plotting
  # Create comparison labels for stat_compare_means
  comparison_list <- list(
    c("Gastric", "IM"),
    c("Gastric", "Duodenal"),
    c("IM", "Duodenal")
  )
  
  # Define custom colors with transparency
  tissue_colors <- c(
    "Gastric"   = "#c85a5a",
    "IM"        = "#8e7cc3",
    "Duodenal"  = "#487cac"
  )
  
  # 8. Create the plot
  p <- ggplot(plot_data, aes(x = Sample_Type, y = Gene_Promoter_Mean)) +
    # Boxplot with custom colors and transparency
    geom_boxplot(
      aes(fill = Sample_Type),
      outlier.shape = NA, 
      alpha = 0.4,
      width = 0.5
    ) +
    
    # Individual points with patient colors
    geom_point(aes(fill = Patient), size = 3, shape = 21, alpha = 0.8) +
    
    # Connect points by patient
    geom_line(aes(group = Patient, color = Patient), alpha = 0.6, linewidth = 0.8) +
    
    # Add pairwise comparisons with STARS ONLY (no text)
    stat_compare_means(
      comparisons = comparison_list,
      method = "t.test",
      label = "p.signif",  
      hide.ns = FALSE,     
      bracket.size = 0.5,
      tip.length = 0.02,
      step.increase = 0.1,
      size = 5  
    ) +
    
    # Custom color scale for Sample_Type (boxplot fill)
    scale_fill_manual(
      values = tissue_colors,
      aesthetics = "fill",
      breaks = c("Gastric", "IM", "Duodenal")  
    ) +
    
    # Patient colors for points and lines
    scale_color_brewer(palette = "Set3") +
    
    # Styling
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
      axis.title = element_text(face = "bold"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste0(gene_name, " Promoter Methylation"),
      subtitle = "Connected lines show matched samples from same patient",
      y = "Mean Beta Value (Promoter)",
      x = "Tissue Type",
      fill = "Tissue Type",
      color = "Patient ID"
    ) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))
  
  # 9. Save plot
  if (save_plot) {
    filename <- file.path(output_dir, paste0(gene_name, "_Promoter_Methylation.png"))
    ggsave(filename, p, width = 10, height = 7, dpi = 300)
    cat("✅ Plot saved:", filename, "\n")
    
    # Also save as PDF for publications
    filename_pdf <- file.path(output_dir, paste0(gene_name, "_Promoter_Methylation.pdf"))
    ggsave(filename_pdf, p, width = 10, height = 7)
  }
  
  # 10. Save statistics
  if (save_plot) {
    stats_file <- file.path(output_dir, paste0(gene_name, "_Statistics.txt"))
    sink(stats_file)
    cat("=====================================\n")
    cat("Gene:", gene_name, "\n")
    cat("Number of promoter probes:", length(gene_promoter_probes), "\n")
    cat("Promoter regions:", paste(promoter_regions, collapse = ", "), "\n")
    cat("=====================================\n\n")
    cat("--- Probe IDs ---\n")
    print(gene_promoter_probes)
    cat("\n--- Mixed-Effects Model Summary ---\n")
    print(summary(model))
    cat("\n--- Pairwise Comparisons (Tukey-adjusted) ---\n")
    print(comparisons)
    sink()
    cat("✅ Statistics saved:", stats_file, "\n")
  }
  
  # Return results
  return(list(
    plot = p,
    model = model,
    comparisons = comparisons,
    data = plot_data,
    n_probes = length(gene_promoter_probes),
    probe_ids = gene_promoter_probes
  ))
}

############################################
# FUNCTION: Batch process multiple genes
############################################

analyze_multiple_genes <- function(gene_list,
                                   beta_matrix = myNorm,
                                   annotation = annEPICv2Sub,
                                   pheno_data = pd,
                                   promoter_regions = c("TSS1500", "TSS200"),
                                   output_dir = "Gene_Methylation_Plots",
                                   create_summary = TRUE) {
  
  cat("\n🔬 BATCH ANALYSIS: Processing", length(gene_list), "genes\n")
  cat("Output directory:", output_dir, "\n\n")
  
  # Store results
  results <- list()
  successful <- c()
  failed <- c()
  
  # Process each gene
  for (gene in gene_list) {
    result <- tryCatch({
      analyze_gene_promoter(
        gene_name = gene,
        beta_matrix = beta_matrix,
        annotation = annotation,
        pheno_data = pheno_data,
        promoter_regions = promoter_regions,
        output_dir = output_dir,
        save_plot = TRUE,
        show_stats = FALSE
      )
    }, error = function(e) {
      cat("❌ Error processing", gene, ":", conditionMessage(e), "\n")
      return(NULL)
    })
    
    if (!is.null(result)) {
      results[[gene]] <- result
      successful <- c(successful, gene)
    } else {
      failed <- c(failed, gene)
    }
  }
  
  # Create summary report
  if (create_summary) {
    summary_file <- file.path(output_dir, "BATCH_SUMMARY.txt")
    sink(summary_file)
    cat("=====================================\n")
    cat("BATCH ANALYSIS SUMMARY\n")
    cat("=====================================\n")
    cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    cat("Total genes requested:", length(gene_list), "\n")
    cat("Successfully processed:", length(successful), "\n")
    cat("Failed:", length(failed), "\n\n")
    
    cat("--- Successfully Processed ---\n")
    for (gene in successful) {
      cat(sprintf("  %s (%d probes)\n", gene, results[[gene]]$n_probes))
    }
    
    if (length(failed) > 0) {
      cat("\n--- Failed ---\n")
      for (gene in failed) {
        cat(" ", gene, "\n")
      }
    }
    
    cat("\n--- Individual Gene Statistics ---\n\n")
    for (gene in successful) {
      cat("=====================================\n")
      cat("Gene:", gene, "\n")
      cat("Number of probes:", results[[gene]]$n_probes, "\n")
      cat("Pairwise comparisons:\n")
      print(results[[gene]]$comparisons)
      cat("\n")
    }
    sink()
    
    cat("\n✅ Summary report saved:", summary_file, "\n")
  }
  
  # Print summary to console
  cat("\n=====================================\n")
  cat("BATCH PROCESSING COMPLETE\n")
  cat("=====================================\n")
  cat("Successfully processed:", length(successful), "genes\n")
  cat("Failed:", length(failed), "genes\n")
  cat("All plots saved to:", output_dir, "\n")
  
  return(results)
}

############################################
# BONUS: Create comparison heatmap
############################################

create_gene_comparison_heatmap <- function(results, output_dir = "Gene_Methylation_Plots") {
  
  # Extract mean values per tissue type for each gene
  gene_means <- do.call(rbind, lapply(names(results), function(gene) {
    data <- results[[gene]]$data
    means <- aggregate(Gene_Promoter_Mean ~ Sample_Type, data = data, FUN = mean)
    values <- setNames(means$Gene_Promoter_Mean, means$Sample_Type)
    return(values)
  }))
  
  rownames(gene_means) <- names(results)
  
  # Create heatmap
  pheatmap(
    gene_means,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = seq(0, 1, length.out = 101),
    main = "Gene Promoter Methylation Across Tissue Types",
    fontsize = 10,
    cellwidth = 40,
    cellheight = 20,
    filename = file.path(output_dir, "Gene_Comparison_Heatmap.png"),
    width = 8,
    height = length(results) * 0.3 + 2
  )
  
  cat("✅ Heatmap saved to:", file.path(output_dir, "Gene_Comparison_Heatmap.png"), "\n")
}

############################################
# USAGE
############################################

# Test with CDX2
cdx2_result <- analyze_gene_promoter("NKX6-3")
print(cdx2_result$plot)

# Batch process multiple genes
genes_of_interest <- c(
  "ATP4A","BHLHA15","CLDN18","GATA4","GATA6","GIF","GKN1","GKN2","KRT18",
  "MUC5AC","MUC6","PGC","SLC26A7","SOX2","TFF1","TFF2","AGR2","ANXA10","AXIN2",
  "CDH17","CDX2","CLDN10","CLDN7","HNF4A","KRT20","LGR5","MUC1","MUC2","OLFM4",
  "REG4","RNF43","SOX9","SPDEF","TFF3","ALPI","APOA1","APOA4","ASCL2","CDX1",
  "FABP1","FABP2","LYZ","SI","SLC26A3","SLC5A1","VIL1"
  )

genes_of_interest_from_RNAseq <- c(
  "A1BG","BGN","BLOC1S5-TXNDC5","CGB5","EPHA2","LIF",
  "LOC102723901", "MAP3K5", "PIK3CG", "SLC2A6","TMX2-CTNND1A",
  "WNT7A", "ZNF625-ZNF20","ZNF664-FAM101A", "HOXA1","ITGB2-AS1",
  "NKX6-3","PCSK1N","ZNF462","CA12","ADAMTSL4","C15orf52","C1orf21","CASC9",
  "COL4A2","FGF19","FLG","GABRP","GRHL3","KRT80","MXRA5","MYEOV","NLRP6","REG1A",
  "SERPINB7","SLC4A3","STAC","TRPM2","ADAMTS14","ALDH1A3","ALDH1A3","AZGP1",
  "CLDN10","GABRE","GCT6","PLXNB3","SHISA6","SIX1","TM6SF2")

all_results <- analyze_multiple_genes(genes_of_interest)
all_results_from_RNAseq <- analyze_multiple_genes(genes_of_interest_from_RNAseq)

create_gene_comparison_heatmap(all_results)







library(missMethyl)
library(ggplot2)
library(dplyr)
library(stringr)

################################################################
# HELPER FUNCTION: Run Analysis & Create Plots
################################################################

run_pathway_analysis <- function(dmp_table,
                                 comparison_name,
                                 array_type = "EPIC_V2") {
  
  cat(paste0("\n=== Processing: ", comparison_name, " ===\n"))
  
  ##############################################################
  # 1. Define Significant & Background Probes
  ##############################################################
  
  sig_probes <- rownames(dmp_table)[dmp_table$adj.P.Val < 0.05]
  all_probes <- rownames(dmp_table)
  
  cat("Significant Probes:", length(sig_probes), "\n")
  
  if (length(sig_probes) < 50) {
    warning(paste("Skipping", comparison_name,
                  "- Not enough significant probes (<50)."))
    return(NULL)
  }
  
  ##############################################################
  # 2. GO Analysis (Biological Process)
  ##############################################################
  
  go_res <- gometh(
    sig.cpg    = sig_probes,
    all.cpg    = all_probes,
    collection = "GO",
    array.type = array_type
  )
  
  top_GO <- topGSA(go_res, number = 20)
  top_GO <- top_GO[top_GO$Ont == "BP", ]
  
  ##############################################################
  # 3. KEGG Analysis
  ##############################################################
  
  kegg_res <- gometh(
    sig.cpg    = sig_probes,
    all.cpg    = all_probes,
    collection = "KEGG",
    array.type = array_type
  )
  
  top_KEGG <- topGSA(kegg_res, number = 20)
  
  ##############################################################
  # 4. Plotting Function
  ##############################################################
  
  create_dotplot <- function(res_df, title) {
    
    res_df$Term <- str_trunc(res_df$Term, 50)
    
    res_df$Term <- factor(
      res_df$Term,
      levels = res_df$Term[order(res_df$FDR, decreasing = TRUE)]
    )
    
    ggplot(res_df, aes(x = DE, y = Term)) +
      geom_point(aes(color = FDR, size = DE)) +
      scale_color_gradient(low = "red", high = "blue") +
      labs(
        title = title,
        x = "Number of Genes",
        y = "Pathway / Term",
        color = "FDR",
        size = "Count"
      ) +
      theme_bw() +
      theme(axis.text.y = element_text(size = 10))
  }
  
  ##############################################################
  # 5. Generate plots
  ##############################################################
  
  p_GO <- create_dotplot(
    top_GO,
    paste("Top GO Terms (BP):", comparison_name)
  )
  
  p_KEGG <- create_dotplot(
    top_KEGG,
    paste("Top KEGG Pathways:", comparison_name)
  )
  
  ##############################################################
  # 6. SAVE PLOTS TO WORKING DIRECTORY
  ##############################################################
  
  safe_name <- gsub(" ", "_", comparison_name)
  
  ggsave(
    filename = paste0("GO_", safe_name, ".pdf"),
    plot     = p_GO,
    width    = 9,
    height   = 7
  )
  
  ggsave(
    filename = paste0("KEGG_", safe_name, ".pdf"),
    plot     = p_KEGG,
    width    = 9,
    height   = 7
  )
  
  ##############################################################
  # Output
  ##############################################################
  
  return(list(
    GO_Table   = top_GO,
    KEGG_Table = top_KEGG,
    Plot_GO    = p_GO,
    Plot_KEGG  = p_KEGG
  ))
}

################################################################
# EXECUTION: Run for all Comparisons
################################################################

results_list <- list()

if (exists("dmp_IM_vs_Gas")) {
  results_list[["IM_vs_Gastric"]] <-
    run_pathway_analysis(dmp_IM_vs_Gas, "IM vs Gastric")
}

if (exists("dmp_Duo_vs_Gas")) {
  results_list[["Duo_vs_Gastric"]] <-
    run_pathway_analysis(dmp_Duo_vs_Gas, "Duodenal vs Gastric")
}

if (exists("dmp_IM_vs_Duo")) {
  results_list[["IM_vs_Duodenal"]] <-
    run_pathway_analysis(dmp_IM_vs_Duo, "IM vs Duodenal")
}

################################################################
# OPTIONAL: SAVE TABLES
################################################################

if (!is.null(results_list$IM_vs_Gastric)) {
  
  write.csv(
    results_list$IM_vs_Gastric$GO_Table,
    "Pathways_GO_IM_vs_Gastric.csv",
    row.names = FALSE
  )
  
  write.csv(
    results_list$IM_vs_Gastric$KEGG_Table,
    "Pathways_KEGG_IM_vs_Gastric.csv",
    row.names = FALSE
  )
}
