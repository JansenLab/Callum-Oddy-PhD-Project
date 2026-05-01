############################################
# DIRECTIONAL Gene Set Analysis
# Separates hyper- vs hypo-methylated genes
############################################

library(missMethyl)
library(org.Hs.eg.db)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

cat("\n=== DIRECTIONAL GENE SET ANALYSIS ===\n")

############################################
# Function for directional GSA
############################################

run_directional_gsa <- function(dmp_results, comparison_name, fdr_threshold = 0.05, logfc_threshold = 0) {
  
  cat("\n--- Processing:", comparison_name, "---\n")
  
  # Split into hyper and hypo based on logFC
  hyper_cpgs <- dmp_results$Name[dmp_results$adj.P.Val < fdr_threshold & dmp_results$logFC > logfc_threshold]
  hypo_cpgs <- dmp_results$Name[dmp_results$adj.P.Val < fdr_threshold & dmp_results$logFC < -logfc_threshold]
  all_cpgs <- dmp_results$Name
  
  cat("  Hypermethylated CpGs:", length(hyper_cpgs), "\n")
  cat("  Hypomethylated CpGs:", length(hypo_cpgs), "\n")
  cat("  Total CpGs tested:", length(all_cpgs), "\n")
  
  # Get annotation
  annEPICv2 <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
  
  # Helper function to process one direction
  process_direction <- function(sig_cpgs, direction_name) {
    
    if (length(sig_cpgs) < 5) {
      cat("  ⚠️  Too few", direction_name, "CpGs\n")
      return(list(GO = NULL, KEGG = NULL))
    }
    
    # Map to genes
    sig_genes <- annEPICv2$UCSC_RefGene_Name[match(sig_cpgs, annEPICv2$Name)]
    all_genes <- annEPICv2$UCSC_RefGene_Name[match(all_cpgs, annEPICv2$Name)]
    
    sig_genes <- unique(unlist(strsplit(sig_genes, ";")))
    all_genes <- unique(unlist(strsplit(all_genes, ";")))
    
    sig_genes <- sig_genes[!is.na(sig_genes) & sig_genes != ""]
    all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]
    
    # Convert to Entrez
    sig_entrez <- mapIds(org.Hs.eg.db, keys = sig_genes, column = "ENTREZID", 
                         keytype = "SYMBOL", multiVals = "first")
    all_entrez <- mapIds(org.Hs.eg.db, keys = all_genes, column = "ENTREZID", 
                         keytype = "SYMBOL", multiVals = "first")
    
    sig_entrez <- sig_entrez[!is.na(sig_entrez)]
    all_entrez <- all_entrez[!is.na(all_entrez)]
    
    cat("  ", direction_name, "genes:", length(sig_entrez), "\n")
    
    # GO analysis
    go_results <- tryCatch({
      go <- goana(de = sig_entrez, universe = all_entrez, species = "Hs")
      go <- go[go$Ont == "BP", ]
      go$FDR <- p.adjust(go$P.DE, method = "BH")
      go <- go[order(go$P.DE), ]
      go
    }, error = function(e) NULL)
    
    # KEGG analysis
    kegg_results <- tryCatch({
      kegg <- kegga(de = sig_entrez, universe = all_entrez, species = "Hs")
      kegg$FDR <- p.adjust(kegg$P.DE, method = "BH")
      kegg <- kegg[order(kegg$P.DE), ]
      kegg
    }, error = function(e) NULL)
    
    if (!is.null(go_results)) {
      cat("  ", direction_name, "GO BP terms (P<0.05):", sum(go_results$P.DE < 0.05), "\n")
    }
    if (!is.null(kegg_results)) {
      cat("  ", direction_name, "KEGG pathways (P<0.05):", sum(kegg_results$P.DE < 0.05), "\n")
    }
    
    return(list(GO = go_results, KEGG = kegg_results))
  }
  
  # Process both directions
  hyper_results <- process_direction(hyper_cpgs, "Hypermethylated")
  hypo_results <- process_direction(hypo_cpgs, "Hypomethylated")
  
  return(list(
    hyper = hyper_results,
    hypo = hypo_results
  ))
}

############################################
# Run directional analysis
############################################

dir_gsa_IM_vs_Gas <- run_directional_gsa(dmp_IM_vs_Gas, "IM vs Gastric")
dir_gsa_Duo_vs_Gas <- run_directional_gsa(dmp_Duo_vs_Gas, "Duodenal vs Gastric")
dir_gsa_IM_vs_Duo <- run_directional_gsa(dmp_IM_vs_Duo, "IM vs Duodenal")

############################################
# Save directional results
############################################

save_directional_results <- function(dir_gsa, comparison_name) {
  
  # Hypermethylated GO
  if (!is.null(dir_gsa$hyper$GO)) {
    write.csv(dir_gsa$hyper$GO, 
              paste0("GSA_HYPER_GO_", gsub(" ", "_", comparison_name), ".csv"))
  }
  
  # Hypermethylated KEGG
  if (!is.null(dir_gsa$hyper$KEGG)) {
    write.csv(dir_gsa$hyper$KEGG, 
              paste0("GSA_HYPER_KEGG_", gsub(" ", "_", comparison_name), ".csv"))
  }
  
  # Hypomethylated GO
  if (!is.null(dir_gsa$hypo$GO)) {
    write.csv(dir_gsa$hypo$GO, 
              paste0("GSA_HYPO_GO_", gsub(" ", "_", comparison_name), ".csv"))
  }
  
  # Hypomethylated KEGG
  if (!is.null(dir_gsa$hypo$KEGG)) {
    write.csv(dir_gsa$hypo$KEGG, 
              paste0("GSA_HYPO_KEGG_", gsub(" ", "_", comparison_name), ".csv"))
  }
}

save_directional_results(dir_gsa_IM_vs_Gas, "IM_vs_Gastric")
save_directional_results(dir_gsa_Duo_vs_Gas, "Duodenal_vs_Gastric")
save_directional_results(dir_gsa_IM_vs_Duo, "IM_vs_Duodenal")

############################################
# Directional dot plot
############################################

create_directional_dotplot <- function(dir_gsa, comparison_name, type = "GO", n_terms = 10) {
  
  # Get results for both directions
  hyper <- if (type == "GO") dir_gsa$hyper$GO else dir_gsa$hyper$KEGG
  hypo <- if (type == "GO") dir_gsa$hypo$GO else dir_gsa$hypo$KEGG
  
  if (is.null(hyper) && is.null(hypo)) {
    cat("No directional", type, "results for", comparison_name, "\n")
    return(NULL)
  }
  
  # Prepare data
  plot_data <- data.frame()
  
  if (!is.null(hyper) && nrow(hyper) > 0) {
    top_hyper <- head(hyper, n_terms)
    top_hyper$Direction <- "Hypermethylated"
    top_hyper$Term_display <- if (type == "GO") top_hyper$Term else top_hyper$Pathway
    top_hyper$GeneRatio <- top_hyper$DE / top_hyper$N
    top_hyper$negLogP <- -log10(top_hyper$P.DE)
    plot_data <- rbind(plot_data, top_hyper[, c("Term_display", "Direction", "DE", "GeneRatio", "negLogP", "P.DE")])
  }
  
  if (!is.null(hypo) && nrow(hypo) > 0) {
    top_hypo <- head(hypo, n_terms)
    top_hypo$Direction <- "Hypomethylated"
    top_hypo$Term_display <- if (type == "GO") top_hypo$Term else top_hypo$Pathway
    top_hypo$GeneRatio <- top_hypo$DE / top_hypo$N
    top_hypo$negLogP <- -log10(top_hypo$P.DE)
    plot_data <- rbind(plot_data, top_hypo[, c("Term_display", "Direction", "DE", "GeneRatio", "negLogP", "P.DE")])
  }
  
  if (nrow(plot_data) == 0) return(NULL)
  
  # Shorten term names
  plot_data$Term_display <- ifelse(
    nchar(plot_data$Term_display) > 55,
    paste0(substr(plot_data$Term_display, 1, 52), "..."),
    plot_data$Term_display
  )
  
  # Order by direction and p-value
  plot_data <- plot_data[order(plot_data$Direction, -plot_data$negLogP), ]
  plot_data$Term_display <- factor(plot_data$Term_display, levels = unique(plot_data$Term_display))
  
  # Create plot
  p <- ggplot(plot_data, aes(x = GeneRatio, y = Term_display)) +
    geom_point(aes(size = DE, color = negLogP)) +
    scale_color_gradient(low = "blue", high = "red", name = "-log10(P)") +
    scale_size_continuous(name = "# Genes", range = c(3, 10)) +
    facet_wrap(~Direction, scales = "free_y", ncol = 1) +
    labs(
      title = paste0(type, " Enrichment (Directional): ", comparison_name),
      subtitle = paste0("Top ", n_terms, " terms for each direction"),
      x = "Gene Ratio",
      y = NULL
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 9),
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "right"
    )
  
  filename <- paste0("GSA_DIRECTIONAL_", type, "_", gsub(" ", "_", comparison_name), ".png")
  ggsave(filename, p, width = 12, height = 12, dpi = 300)
  cat("✅ Saved directional plot:", filename, "\n")
  
  return(p)
}

############################################
# Create comparison heatmap
############################################

create_comparison_heatmap <- function(gsa_list, type = "GO", top_n = 30) {
  
  cat("\n--- Creating comparison heatmap for", type, "---\n")
  
  # Extract results from all comparisons
  all_terms <- list()
  
  for (comp_name in names(gsa_list)) {
    results <- if (type == "GO") gsa_list[[comp_name]]$GO else gsa_list[[comp_name]]$KEGG
    
    if (!is.null(results) && nrow(results) > 0) {
      # Get term names
      results$Term_name <- if (type == "GO") results$Term else results$Pathway
      
      # Get top terms by p-value
      top <- head(results[order(results$P.DE), ], top_n)
      all_terms[[comp_name]] <- setNames(top$P.DE, top$Term_name)
    }
  }
  
  if (length(all_terms) == 0) {
    cat("No results to create heatmap\n")
    return(NULL)
  }
  
  # Get union of all terms
  unique_terms <- unique(unlist(lapply(all_terms, names)))
  
  # Create matrix
  heatmap_matrix <- matrix(NA, nrow = length(unique_terms), ncol = length(all_terms))
  rownames(heatmap_matrix) <- unique_terms
  colnames(heatmap_matrix) <- names(all_terms)
  
  # Fill matrix with -log10(P-values)
  for (comp in names(all_terms)) {
    terms <- names(all_terms[[comp]])
    pvals <- all_terms[[comp]]
    heatmap_matrix[terms, comp] <- -log10(pvals)
  }
  
  # Keep only terms that appear in at least one comparison with P < 0.05
  sig_rows <- apply(heatmap_matrix, 1, function(x) any(x > -log10(0.05), na.rm = TRUE))
  heatmap_matrix <- heatmap_matrix[sig_rows, , drop = FALSE]
  
  # Take top 40 by max significance across comparisons
  row_max <- apply(heatmap_matrix, 1, max, na.rm = TRUE)
  top_rows <- head(order(row_max, decreasing = TRUE), 40)
  heatmap_matrix <- heatmap_matrix[top_rows, , drop = FALSE]
  
  # Shorten row names
  rownames(heatmap_matrix) <- ifelse(
    nchar(rownames(heatmap_matrix)) > 60,
    paste0(substr(rownames(heatmap_matrix), 1, 57), "..."),
    rownames(heatmap_matrix)
  )
  
  # Replace NA with 0 for visualization
  heatmap_matrix[is.na(heatmap_matrix)] <- 0
  
  # Create heatmap
  filename <- paste0("GSA_Comparison_Heatmap_", type, ".png")
  
  png(filename, width = 12, height = 14, units = "in", res = 300)
  
  pheatmap(
    heatmap_matrix,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100),
    breaks = seq(0, max(heatmap_matrix, na.rm = TRUE), length.out = 101),
    main = paste0(type, " Enrichment Across Comparisons\n(-log10 P-value)"),
    fontsize_row = 8,
    fontsize_col = 11,
    angle_col = 45,
    border_color = "gray80",
    na_col = "gray95"
  )
  
  dev.off()
  cat("✅ Saved heatmap:", filename, "\n")
}

############################################
# Generate all directional plots and heatmaps
############################################

cat("\n=== CREATING DIRECTIONAL VISUALIZATIONS ===\n")

# Directional plots for each comparison
create_directional_dotplot(dir_gsa_IM_vs_Gas, "IM vs Gastric", "GO", 10)
create_directional_dotplot(dir_gsa_IM_vs_Gas, "IM vs Gastric", "KEGG", 10)

create_directional_dotplot(dir_gsa_Duo_vs_Gas, "Duodenal vs Gastric", "GO", 10)
create_directional_dotplot(dir_gsa_Duo_vs_Gas, "Duodenal vs Gastric", "KEGG", 10)

create_directional_dotplot(dir_gsa_IM_vs_Duo, "IM vs Duodenal", "GO", 10)
create_directional_dotplot(dir_gsa_IM_vs_Duo, "IM vs Duodenal", "KEGG", 10)

# Comparison heatmaps
gsa_list <- list(
  "IM vs Gastric" = gsa_IM_vs_Gas,
  "Duodenal vs Gastric" = gsa_Duo_vs_Gas,
  "IM vs Duodenal" = gsa_IM_vs_Duo
)

create_comparison_heatmap(gsa_list, "GO", top_n = 50)
create_comparison_heatmap(gsa_list, "KEGG", top_n = 50)

cat("\n✅ ALL DIRECTIONAL ANALYSES COMPLETE!\n")





############################################
# Summary Statistics and Overview Figures
############################################

library(ggplot2)
library(gridExtra)
library(scales)

cat("\n=== CREATING SUMMARY STATISTICS AND FIGURES ===\n")

############################################
# Function to extract summary statistics
############################################

extract_summary_stats <- function(gsa_results, dir_gsa_results, comparison_name) {
  
  # Overall enrichment
  go_total <- if (!is.null(gsa_results$GO)) nrow(gsa_results$GO) else 0
  go_sig <- if (!is.null(gsa_results$GO)) sum(gsa_results$GO$P.DE < 0.05) else 0
  
  kegg_total <- if (!is.null(gsa_results$KEGG)) nrow(gsa_results$KEGG) else 0
  kegg_sig <- if (!is.null(gsa_results$KEGG)) sum(gsa_results$KEGG$P.DE < 0.05) else 0
  
  # Directional enrichment
  go_hyper_sig <- if (!is.null(dir_gsa_results$hyper$GO)) sum(dir_gsa_results$hyper$GO$P.DE < 0.05) else 0
  go_hypo_sig <- if (!is.null(dir_gsa_results$hypo$GO)) sum(dir_gsa_results$hypo$GO$P.DE < 0.05) else 0
  
  kegg_hyper_sig <- if (!is.null(dir_gsa_results$hyper$KEGG)) sum(dir_gsa_results$hyper$KEGG$P.DE < 0.05) else 0
  kegg_hypo_sig <- if (!is.null(dir_gsa_results$hypo$KEGG)) sum(dir_gsa_results$hypo$KEGG$P.DE < 0.05) else 0
  
  # Top terms/pathways
  go_top5 <- if (!is.null(gsa_results$GO) && nrow(gsa_results$GO) > 0) {
    top <- head(gsa_results$GO[order(gsa_results$GO$P.DE), ], 5)
    data.frame(
      Term = top$Term,
      P.DE = top$P.DE,
      FDR = top$FDR,
      N_genes = top$N,
      DE_genes = top$DE
    )
  } else NULL
  
  kegg_top5 <- if (!is.null(gsa_results$KEGG) && nrow(gsa_results$KEGG) > 0) {
    top <- head(gsa_results$KEGG[order(gsa_results$KEGG$P.DE), ], 5)
    data.frame(
      Pathway = top$Pathway,
      P.DE = top$P.DE,
      FDR = top$FDR,
      N_genes = top$N,
      DE_genes = top$DE
    )
  } else NULL
  
  return(list(
    comparison = comparison_name,
    go_total = go_total,
    go_sig = go_sig,
    kegg_total = kegg_total,
    kegg_sig = kegg_sig,
    go_hyper_sig = go_hyper_sig,
    go_hypo_sig = go_hypo_sig,
    kegg_hyper_sig = kegg_hyper_sig,
    kegg_hypo_sig = kegg_hypo_sig,
    go_top5 = go_top5,
    kegg_top5 = kegg_top5
  ))
}

############################################
# Extract stats for all comparisons
############################################

stats_IM_Gas <- extract_summary_stats(gsa_IM_vs_Gas, dir_gsa_IM_vs_Gas, "IM vs Gastric")
stats_Duo_Gas <- extract_summary_stats(gsa_Duo_vs_Gas, dir_gsa_Duo_vs_Gas, "Duodenal vs Gastric")
stats_IM_Duo <- extract_summary_stats(gsa_IM_vs_Duo, dir_gsa_IM_vs_Duo, "IM vs Duodenal")

all_stats <- list(stats_IM_Gas, stats_Duo_Gas, stats_IM_Duo)

############################################
# Print summary statistics to console
############################################

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║          GENE SET ENRICHMENT ANALYSIS - SUMMARY               ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

for (stats in all_stats) {
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat(stats$comparison, "\n")
  cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")
  
  cat("GO BIOLOGICAL PROCESS:\n")
  cat("  Total terms tested:", stats$go_total, "\n")
  cat("  Significantly enriched (P < 0.05):", stats$go_sig, "\n")
  cat("    └─ Hypermethylated:", stats$go_hyper_sig, "\n")
  cat("    └─ Hypomethylated:", stats$go_hypo_sig, "\n\n")
  
  cat("KEGG PATHWAYS:\n")
  cat("  Total pathways tested:", stats$kegg_total, "\n")
  cat("  Significantly enriched (P < 0.05):", stats$kegg_sig, "\n")
  cat("    └─ Hypermethylated:", stats$kegg_hyper_sig, "\n")
  cat("    └─ Hypomethylated:", stats$kegg_hypo_sig, "\n\n")
  
  if (!is.null(stats$go_top5)) {
    cat("TOP 5 GO TERMS:\n")
    for (i in 1:min(5, nrow(stats$go_top5))) {
      term <- stats$go_top5[i, ]
      cat(sprintf("  %d. %s\n", i, substr(term$Term, 1, 70)))
      cat(sprintf("     P = %.2e | FDR = %.2e | %d/%d genes\n", 
                  term$P.DE, term$FDR, term$DE_genes, term$N_genes))
    }
    cat("\n")
  }
  
  if (!is.null(stats$kegg_top5)) {
    cat("TOP 5 KEGG PATHWAYS:\n")
    for (i in 1:min(5, nrow(stats$kegg_top5))) {
      pathway <- stats$kegg_top5[i, ]
      cat(sprintf("  %d. %s\n", i, substr(pathway$Pathway, 1, 70)))
      cat(sprintf("     P = %.2e | FDR = %.2e | %d/%d genes\n", 
                  pathway$P.DE, pathway$FDR, pathway$DE_genes, pathway$N_genes))
    }
    cat("\n")
  }
  
  cat("\n")
}

############################################
# Create summary bar chart
############################################

create_summary_barplot <- function(all_stats) {
  
  # Prepare data
  summary_data <- data.frame(
    Comparison = rep(sapply(all_stats, function(x) x$comparison), each = 2),
    Type = rep(c("GO BP", "KEGG"), times = length(all_stats)),
    Significant = c(
      sapply(all_stats, function(x) c(x$go_sig, x$kegg_sig))
    )
  )
  
  summary_data$Comparison <- factor(summary_data$Comparison, 
                                    levels = c("IM vs Gastric", "Duodenal vs Gastric", "IM vs Duodenal"))
  
  # Create plot
  p <- ggplot(summary_data, aes(x = Comparison, y = Significant, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_text(aes(label = Significant), 
              position = position_dodge(width = 0.7), 
              vjust = -0.5, size = 4) +
    scale_fill_manual(values = c("GO BP" = "#3498db", "KEGG" = "#e74c3c")) +
    labs(
      title = "Significantly Enriched Pathways Across Comparisons",
      subtitle = "Number of GO Biological Process terms and KEGG pathways (P < 0.05)",
      x = NULL,
      y = "Number of Significantly Enriched Terms/Pathways",
      fill = "Database"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "gray40"),
      axis.text.x = element_text(size = 11, angle = 0),
      axis.text.y = element_text(size = 10),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      panel.grid.major.x = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
  ggsave("GSA_Summary_BarChart.png", p, width = 10, height = 6, dpi = 300)
  cat("✅ Saved: GSA_Summary_BarChart.png\n")
  
  return(p)
}

############################################
# Create directional stacked bar chart
############################################

create_directional_summary <- function(all_stats) {
  
  # Prepare data for GO
  go_data <- data.frame(
    Comparison = rep(sapply(all_stats, function(x) x$comparison), each = 2),
    Direction = rep(c("Hypermethylated", "Hypomethylated"), times = length(all_stats)),
    Count = c(sapply(all_stats, function(x) c(x$go_hyper_sig, x$go_hypo_sig)))
  )
  
  go_data$Comparison <- factor(go_data$Comparison, 
                               levels = c("IM vs Gastric", "Duodenal vs Gastric", "IM vs Duodenal"))
  
  # Prepare data for KEGG
  kegg_data <- data.frame(
    Comparison = rep(sapply(all_stats, function(x) x$comparison), each = 2),
    Direction = rep(c("Hypermethylated", "Hypomethylated"), times = length(all_stats)),
    Count = c(sapply(all_stats, function(x) c(x$kegg_hyper_sig, x$kegg_hypo_sig)))
  )
  
  kegg_data$Comparison <- factor(kegg_data$Comparison, 
                                 levels = c("IM vs Gastric", "Duodenal vs Gastric", "IM vs Duodenal"))
  
  # Create GO plot
  p_go <- ggplot(go_data, aes(x = Comparison, y = Count, fill = Direction)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_text(aes(label = Count), 
              position = position_dodge(width = 0.7), 
              vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = c("Hypermethylated" = "#d62728", "Hypomethylated" = "#2ca02c")) +
    labs(
      title = "GO BP Enrichment by Direction",
      x = NULL,
      y = "Number of Enriched Terms (P < 0.05)",
      fill = "Direction"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(size = 10, angle = 0),
      legend.position = "top"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
  # Create KEGG plot
  p_kegg <- ggplot(kegg_data, aes(x = Comparison, y = Count, fill = Direction)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_text(aes(label = Count), 
              position = position_dodge(width = 0.7), 
              vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = c("Hypermethylated" = "#d62728", "Hypomethylated" = "#2ca02c")) +
    labs(
      title = "KEGG Pathway Enrichment by Direction",
      x = NULL,
      y = "Number of Enriched Pathways (P < 0.05)",
      fill = "Direction"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(size = 10, angle = 0),
      legend.position = "top"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
  # Combine plots
  combined <- grid.arrange(p_go, p_kegg, ncol = 1)
  
  ggsave("GSA_Directional_Summary.png", combined, width = 10, height = 10, dpi = 300)
  cat("✅ Saved: GSA_Directional_Summary.png\n")
  
  return(combined)
}

############################################
# Create comprehensive summary table
############################################

create_summary_table_plot <- function(all_stats) {
  
  # Create table data
  table_data <- data.frame(
    Comparison = character(),
    Metric = character(),
    Value = character(),
    stringsAsFactors = FALSE
  )
  
  for (stats in all_stats) {
    table_data <- rbind(table_data,
                        data.frame(
                          Comparison = stats$comparison,
                          Metric = c(
                            "GO BP (significant)",
                            "GO BP (hypermeth)",
                            "GO BP (hypometh)",
                            "KEGG (significant)",
                            "KEGG (hypermeth)",
                            "KEGG (hypometh)"
                          ),
                          Value = c(
                            stats$go_sig,
                            stats$go_hyper_sig,
                            stats$go_hypo_sig,
                            stats$kegg_sig,
                            stats$kegg_hyper_sig,
                            stats$kegg_hypo_sig
                          )
                        )
    )
  }
  
  # Reshape for table
  table_wide <- reshape(table_data, 
                        idvar = "Metric", 
                        timevar = "Comparison", 
                        direction = "wide")
  
  colnames(table_wide) <- gsub("Value.", "", colnames(table_wide))
  
  # Create visual table
  library(gridExtra)
  
  # Format the table
  tt <- ttheme_default(
    core = list(
      fg_params = list(fontsize = 10),
      bg_params = list(fill = c(rep(c("white", "gray95"), length.out = 6)))
    ),
    colhead = list(
      fg_params = list(fontsize = 11, fontface = "bold"),
      bg_params = list(fill = "lightblue")
    ),
    rowhead = list(
      fg_params = list(fontsize = 10, fontface = "bold")
    )
  )
  
  table_grob <- tableGrob(table_wide, rows = NULL, theme = tt)
  
  # Add title
  title <- textGrob("Gene Set Enrichment Summary Statistics", 
                    gp = gpar(fontsize = 14, fontface = "bold"))
  
  combined <- grid.arrange(title, table_grob, 
                           heights = c(0.1, 0.9),
                           ncol = 1)
  
  ggsave("GSA_Summary_Table.png", combined, width = 12, height = 5, dpi = 300)
  cat("✅ Saved: GSA_Summary_Table.png\n")
  
  # Also save as CSV
  write.csv(table_wide, "GSA_Summary_Statistics.csv", row.names = FALSE)
  cat("✅ Saved: GSA_Summary_Statistics.csv\n")
  
  return(combined)
}

############################################
# Create pie charts showing proportions
############################################

create_proportion_pies <- function(all_stats) {
  
  par(mfrow = c(2, 3), mar = c(2, 2, 3, 2))
  
  for (i in 1:length(all_stats)) {
    stats <- all_stats[[i]]
    
    # GO pie chart
    if (stats$go_sig > 0) {
      go_values <- c(stats$go_hyper_sig, stats$go_hypo_sig)
      go_labels <- c(
        paste0("Hyper: ", stats$go_hyper_sig),
        paste0("Hypo: ", stats$go_hypo_sig)
      )
      go_pct <- round(100 * go_values / sum(go_values), 1)
      go_labels <- paste(go_labels, "\n(", go_pct, "%)", sep = "")
      
      pie(go_values, 
          labels = go_labels, 
          col = c("#d62728", "#2ca02c"),
          main = paste0("GO BP: ", stats$comparison))
    }
    
    # KEGG pie chart
    if (stats$kegg_sig > 0) {
      kegg_values <- c(stats$kegg_hyper_sig, stats$kegg_hypo_sig)
      kegg_labels <- c(
        paste0("Hyper: ", stats$kegg_hyper_sig),
        paste0("Hypo: ", stats$kegg_hypo_sig)
      )
      kegg_pct <- round(100 * kegg_values / sum(kegg_values), 1)
      kegg_labels <- paste(kegg_labels, "\n(", kegg_pct, "%)", sep = "")
      
      pie(kegg_values, 
          labels = kegg_labels, 
          col = c("#d62728", "#2ca02c"),
          main = paste0("KEGG: ", stats$comparison))
    }
  }
  
  dev.print(png, "GSA_Proportion_Pies.png", width = 1200, height = 800, res = 100)
  cat("✅ Saved: GSA_Proportion_Pies.png\n")
  
  par(mfrow = c(1, 1))
}

############################################
# Generate all summary visualizations
############################################

cat("\n--- Generating summary visualizations ---\n")

create_summary_barplot(all_stats)
create_directional_summary(all_stats)
create_summary_table_plot(all_stats)
create_proportion_pies(all_stats)

cat("\n✅ ALL SUMMARY STATISTICS AND FIGURES COMPLETE!\n")
cat("\nGenerated files:\n")
cat("  1. GSA_Summary_BarChart.png - Overall enrichment comparison\n")
cat("  2. GSA_Directional_Summary.png - Hyper vs hypo breakdown\n")
cat("  3. GSA_Summary_Table.png - Statistical summary table\n")
cat("  4. GSA_Summary_Statistics.csv - Raw numbers in CSV format\n")
cat("  5. GSA_Proportion_Pies.png - Proportion pie charts\n")



############################################
# Create pie charts showing proportions
############################################

create_proportion_pies <- function(all_stats) {
  
  # Open PNG device directly
  png("GSA_Proportion_Pies.png", width = 1200, height = 800, res = 100)
  
  par(mfrow = c(2, 3), mar = c(2, 2, 3, 2))
  
  for (i in 1:length(all_stats)) {
    stats <- all_stats[[i]]
    
    # GO pie chart
    if (stats$go_sig > 0) {
      go_values <- c(stats$go_hyper_sig, stats$go_hypo_sig)
      go_labels <- c(
        paste0("Hyper: ", stats$go_hyper_sig),
        paste0("Hypo: ", stats$go_hypo_sig)
      )
      go_pct <- round(100 * go_values / sum(go_values), 1)
      go_labels <- paste(go_labels, "\n(", go_pct, "%)", sep = "")
      
      pie(go_values, 
          labels = go_labels, 
          col = c("#d62728", "#2ca02c"),
          main = paste0("GO BP: ", stats$comparison))
    }
  }
  
  for (i in 1:length(all_stats)) {
    stats <- all_stats[[i]]
    
    # KEGG pie chart
    if (stats$kegg_sig > 0) {
      kegg_values <- c(stats$kegg_hyper_sig, stats$kegg_hypo_sig)
      kegg_labels <- c(
        paste0("Hyper: ", stats$kegg_hyper_sig),
        paste0("Hypo: ", stats$kegg_hypo_sig)
      )
      kegg_pct <- round(100 * kegg_values / sum(kegg_values), 1)
      kegg_labels <- paste(kegg_labels, "\n(", kegg_pct, "%)", sep = "")
      
      pie(kegg_values, 
          labels = kegg_labels, 
          col = c("#d62728", "#2ca02c"),
          main = paste0("KEGG: ", stats$comparison))
    }
  }
  
  dev.off()
  
  par(mfrow = c(1, 1))
  
  cat("✅ Saved: GSA_Proportion_Pies.png\n")
}

# Run it again
create_proportion_pies(all_stats)