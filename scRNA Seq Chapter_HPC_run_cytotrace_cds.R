 --- Environment Setup ---
.libPaths("/home/regmddy/ACFS/Programmes/Packages")
library(Seurat)
library(CytoTRACE2)
library(ggplot2)

# --- Load Data ---
input_file <- "celltypist_input/seu_lineage_labeled_v2.rds"
cat("  Loading Seurat object:", input_file, "\n")
seu <- readRDS(input_file)

# --- Run CytoTRACE2 ---
# We use a batch size of 10k to process your 80k cells safely
cat("  Running CytoTRACE2. Processing 80,000 cells in batches...\n")
ct2_result <- cytotrace2(
  seu,
  species = "human",
  is_seurat = TRUE,
  slot_type = "counts",
  batch_size = 10000,
  ncores = 8,
  seed = 42
)

# --- Add Scores to Metadata ---
seu$CytoTRACE2_Score <- ct2_result$CytoTRACE2_Score
seu$CytoTRACE2_Potency <- ct2_result$CytoTRACE2_Potency

# --- Save Updated Object ---
output_rds <- "celltypist_input/seu_cytotrace_complete.rds"
saveRDS(seu, output_rds)
cat("  Updated Seurat object saved to:", output_rds, "\n")

# --- Generate Diagnostic Plots ---
# --- 1. Extract Data into 'Simple' Formats ---
cat("  Bridging Seurat v5 to CytoTRACE2 dependencies...\n")

# Extract counts as a matrix (bypasses Seurat v5 'Layer' logic)
# Note: This may be memory intensive for 80k cells
exp_matrix <- as.matrix(GetAssayData(seu, assay = "RNA", layer = "counts"))

# Create the specific annotation dataframe required by plotData
anno_df <- data.frame(
  phenotype = as.character(seu$paper_lineage),
  row.names = colnames(seu)
)


# --- 2. Run plotData (The Official Way, but with is_seurat = FALSE) ---
# This forces the package to use the matrix we just made,
# skipping the broken Seurat-compatibility code.
cat("  Generating official CytoTRACE2 plots...\n")

plots <- plotData(
  cytotrace2_result = ct2_result,
  annotation = anno_df,
  expression_data = exp_matrix,
  is_seurat = TRUE
)

# --- 3. Save to PDF ---
pdf("CytoTRACE2_Official_Results.pdf", width = 12, height = 10)
  print(plots$CytoTRACE2_UMAP)
  print(plots$CytoTRACE2_Potency_UMAP)
  print(plots$CytoTRACE2_Relative_UMAP)
  print(plots$Phenotype_UMAP)
  print(plots$CytoTRACE2_Boxplot_byPheno)
dev.off()


