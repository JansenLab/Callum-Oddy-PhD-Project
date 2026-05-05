#!/usr/bin/Rscript --no-save
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db) # Human gene database
})

# ============================================================================
# 1. SETUP PATHS
# ============================================================================
cat("Starting Pathway Enrichment Analysis...\n")
csv_dir  <- "celltypist_input/lineage_results/deep_dive/csv_tables"
out_dir  <- "celltypist_input/lineage_results/deep_dive/pathways"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Find all Level 3 CSV files we generated in the last script
level3_files <- list.files(csv_dir, pattern = "Level3_Intra_.*_Diffs.csv", full.names = TRUE)

if (length(level3_files) == 0) {
  stop("No Level 3 CSV files found. Did the previous script finish successfully?")
}

# ============================================================================
# 2. RUN GO ENRICHMENT LOOP
# ============================================================================
for (file in level3_files) {

  # Extract the lineage name from the filename (e.g., "Neck_Mucous")
  lin_name <- gsub("Level3_Intra_|_Diffs.csv", "", basename(file))
  cat(paste0("\n--- Analyzing Pathways for: ", lin_name, " ---\n"))

  # Load the markers
  markers <- read.csv(file)

  # Filter for strong, significant markers (p < 0.05 and good fold change)
  sig_markers <- markers %>%
    filter(p_val_adj < 0.05 & avg_log2FC > 0.25)

  if (nrow(sig_markers) < 10) {
    cat("  -> Not enough significant genes for pathway analysis. Skipping.\n")
    next
  }

  cat("  -> Converting Gene Symbols to Entrez IDs...\n")
  # clusterProfiler requires 'Entrez IDs' instead of standard gene names (like MUC6)
  gene_conversion <- bitr(sig_markers$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb $

  # Merge the Entrez IDs back with the cluster information
  sig_markers <- sig_markers %>%
    inner_join(gene_conversion, by = c("gene" = "SYMBOL"))

  cat("  -> Running GO Biological Process Enrichment...\n")
  # compareCluster puts all clusters side-by-side for easy visual comparison
  tryCatch({
    go_compare <- compareCluster(
      ENTREZID ~ cluster,
      data = sig_markers,
      fun = "enrichGO",
      OrgDb = org.Hs.eg.db,
      ont = "BP",         # BP = Biological Process (e.g., "cell division", "metabolism")
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      readable = TRUE     # Translates Entrez IDs back to Gene Symbols in the output CSV
    )

    # Save the raw pathway data
    write.csv(as.data.frame(go_compare), file.path(out_dir, paste0("GO_Pathways_", lin_name,$

    # Plot and save the Dotplot
    cat("  -> Generating Plot...\n")
    p <- dotplot(go_compare, showCategory = 5) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
      labs(title = paste("Biological Pathways:", lin_name),
           subtitle = "Top 5 Biological Processes per Sub-cluster",
           x = "Seurat Cluster #")

    ggsave(file.path(out_dir, paste0("GO_Plot_", lin_name, ".pdf")), plot = p, width = 12, h$
    cat("  -> Success!\n")

  }, error = function(e) {
    cat("  -> Error calculating pathways (often caused by too few genes mapping to pathways)$
  })
}

cat("\nAll Pathway analyses finished! Check the 'pathways/' folder.\n")


