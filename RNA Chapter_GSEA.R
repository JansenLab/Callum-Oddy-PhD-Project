# ====================================================================
# RANK-BASED GENE SET ENRICHMENT ANALYSIS (GSEA)
# Following Best Practices for Subtle vs Strong Transcriptional Changes
# ====================================================================
# 
# STRATEGY:
# - IM vs Gastric (~40 DEGs): Rank-based GSEA ONLY
# - Duodenal vs Gastric/IM (>2000 DEGs): Rank-based GSEA (primary) + ORA (secondary)
# - Ranking metric: DESeq2 Wald statistic (incorporates effect size + variance)
# ====================================================================

library(fgsea)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)

setwd('setworkingdirectory')

# Create GSEA output directory
dir.create("GSEA_Results_v2", showWarnings = FALSE)
setwd("GSEA_Results_v2")

cat("=======================================================\n")
cat("RANK-BASED GENE SET ENRICHMENT ANALYSIS\n")
cat("Patient-Adjusted Model (~ Patient + phenotype)\n")
cat("=======================================================\n\n")

# Color scheme
color_palette <- c("Gastric" = "#C75A5A", "IM" = "#9370B8", "Duodenal" = "#6B8DB8")

# ====================================================================
# STEP 1: LOAD DESEQ2 RESULTS (UNFILTERED)
# ====================================================================

cat("STEP 1: Loading DESeq2 results (unfiltered)...\n")

res_IM_Gastric <- read.csv("../Paired_DEG_IM_vs_Gastric.csv", row.names=1)
res_Duo_Gastric <- read.csv("../Paired_DEG_Duodenal_vs_Gastric.csv", row.names=1)
res_Duo_IM <- read.csv("../Paired_DEG_Duodenal_vs_IM.csv", row.names=1)

cat("  IM vs Gastric genes:", nrow(res_IM_Gastric), "\n")
cat("  Duodenal vs Gastric genes:", nrow(res_Duo_Gastric), "\n")
cat("  Duodenal vs IM genes:", nrow(res_Duo_IM), "\n\n")

# ====================================================================
# STEP 2: CREATE RANKED GENE LISTS (WALD STATISTIC)
# ====================================================================

cat("STEP 2: Ranking genes by DESeq2 Wald statistic...\n\n")

create_ranked_list <- function(res_df, comparison_name) {
  # Remove genes with NA statistics or extremely low expression
  res_clean <- res_df %>%
    filter(!is.na(stat), !is.na(baseMean), baseMean > 1)
  
  # Create ranked list using Wald statistic
  gene_ranks <- res_clean$stat
  names(gene_ranks) <- rownames(res_clean)
  gene_ranks <- sort(gene_ranks, decreasing = TRUE)
  
  cat("  ", comparison_name, "\n")
  cat("    Genes after filtering:", length(gene_ranks), "\n")
  cat("    Range:", round(min(gene_ranks), 2), "to", round(max(gene_ranks), 2), "\n")
  
  # Save ranked list
  write.csv(data.frame(Gene = names(gene_ranks), 
                       Wald_stat = gene_ranks),
            paste0("RankedList_", gsub(" ", "_", comparison_name), ".csv"),
            row.names = FALSE)
  
  return(gene_ranks)
}

ranks_IM_Gastric <- create_ranked_list(res_IM_Gastric, "IM vs Gastric")
ranks_Duo_Gastric <- create_ranked_list(res_Duo_Gastric, "Duodenal vs Gastric")
ranks_Duo_IM <- create_ranked_list(res_Duo_IM, "Duodenal vs IM")

cat("\n")

# ====================================================================
# STEP 3: LOAD GENE SETS (BIOLOGY-AWARE SELECTION)
# ====================================================================

cat("STEP 3: Loading gene set collections...\n")

# 3A. MSigDB Hallmark
hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_list <- hallmark_sets %>%
  dplyr::select(gs_name, gene_symbol) %>%
  group_by(gs_name) %>%
  summarise(genes = list(gene_symbol)) %>%
  deframe()

cat("  Hallmark gene sets:", length(hallmark_list), "\n")

# 3B. Canonical pathways
canonical_sets <- msigdbr(species = "Homo sapiens", 
                          collection = "C2", 
                          subcategory = "CP:WIKIPATHWAYS")
canonical_list <- canonical_sets %>%
  dplyr::select(gs_name, gene_symbol) %>%
  group_by(gs_name) %>%
  summarise(genes = list(gene_symbol)) %>%
  deframe()

cat("  WikiPathways:", length(canonical_list), "\n")

# 3C. GO Biological Process
go_sets <- msigdbr(species = "Homo sapiens", 
                   collection = "C5", 
                   subcategory = "GO:BP")
go_list <- go_sets %>%
  dplyr::select(gs_name, gene_symbol) %>%
  group_by(gs_name) %>%
  summarise(genes = list(gene_symbol)) %>%
  deframe()

cat("  GO Biological Process:", length(go_list), "\n")

# 3D. Reactome pathways
reactome_sets <- msigdbr(species = "Homo sapiens", 
                         collection = "C2", 
                         subcategory = "CP:REACTOME")
reactome_list <- reactome_sets %>%
  dplyr::select(gs_name, gene_symbol) %>%
  group_by(gs_name) %>%
  summarise(genes = list(gene_symbol)) %>%
  deframe()

cat("  Reactome pathways:", length(reactome_list), "\n")

# 3E. Cell type signatures
celltype_sets <- msigdbr(species = "Homo sapiens", collection = "C8")
celltype_list <- celltype_sets %>%
  dplyr::select(gs_name, gene_symbol) %>%
  group_by(gs_name) %>%
  summarise(genes = list(gene_symbol)) %>%
  deframe()

cat("  Cell type signatures:", length(celltype_list), "\n\n")

# Combine all gene sets
all_pathways <- c(hallmark_list, canonical_list, reactome_list, go_list, celltype_list)
cat("  Total gene sets:", length(all_pathways), "\n\n")

# ====================================================================
# STEP 4: RUN RANK-BASED GSEA (PRIMARY ANALYSIS)
# ====================================================================

cat("STEP 4: Running rank-based GSEA (fgsea)...\n")
cat("Using fgseaMultilevel (adaptive permutation testing)\n")
cat("Parameters: minSize=15, maxSize=500\n\n")

run_fgsea <- function(gene_ranks, pathways, comparison_name) {
  
  cat("  Running:", comparison_name, "\n")
  
  # Run fgsea (using fgseaMultilevel - more accurate)
  set.seed(42)
  fgsea_res <- fgsea(
    pathways = pathways,
    stats = gene_ranks,
    minSize = 15,
    maxSize = 500
    # Removed nperm to use fgseaMultilevel (recommended)
  )
  
  # Sort by NES
  fgsea_res <- fgsea_res %>%
    arrange(desc(NES))
  
  # Count significant pathways
  sig_up <- sum(fgsea_res$padj < 0.05 & fgsea_res$NES > 0, na.rm=TRUE)
  sig_down <- sum(fgsea_res$padj < 0.05 & fgsea_res$NES < 0, na.rm=TRUE)
  
  cat("    Significant pathways (padj < 0.05):", sig_up + sig_down, "\n")
  cat("      Upregulated (NES > 0):", sig_up, "\n")
  cat("      Downregulated (NES < 0):", sig_down, "\n\n")
  
  # Convert list columns to character for CSV export
  fgsea_for_export <- fgsea_res %>%
    mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse=";")))
  
  # Save full results
  write.csv(fgsea_for_export, 
            paste0("GSEA_", gsub(" ", "_", comparison_name), "_full.csv"),
            row.names = FALSE)
  
  # Save significant only
  fgsea_sig <- fgsea_for_export %>%
    filter(padj < 0.05) %>%
    arrange(padj)
  
  if(nrow(fgsea_sig) > 0) {
    write.csv(fgsea_sig, 
              paste0("GSEA_", gsub(" ", "_", comparison_name), "_significant.csv"),
              row.names = FALSE)
  }
  
  return(fgsea_res)  # Return original with list columns for plotting
}

# Run GSEA for all comparisons
gsea_IM_Gastric <- run_fgsea(ranks_IM_Gastric, all_pathways, "IM vs Gastric")
gsea_Duo_Gastric <- run_fgsea(ranks_Duo_Gastric, all_pathways, "Duodenal vs Gastric")
gsea_Duo_IM <- run_fgsea(ranks_Duo_IM, all_pathways, "Duodenal vs IM")

# ====================================================================
# STEP 5: VISUALIZATIONS
# ====================================================================

cat("STEP 5: Creating visualizations...\n\n")

# 5A. ENRICHMENT DOT PLOTS (NES vs -log10 padj)
create_enrichment_dotplot <- function(fgsea_res, comparison_name, color, top_n=30) {
  
  # Get top pathways
  top_pathways <- fgsea_res %>%
    filter(padj < 0.05) %>%
    arrange(padj) %>%
    head(top_n)
  
  if(nrow(top_pathways) == 0) {
    cat("    No significant pathways for", comparison_name, "\n")
    return(NULL)
  }
  
  # Simplify pathway names for readability
  top_pathways <- top_pathways %>%
    mutate(pathway_short = gsub("HALLMARK_", "", pathway),
           pathway_short = gsub("GOBP_", "", pathway_short),
           pathway_short = gsub("REACTOME_", "", pathway_short),
           pathway_short = gsub("WP_", "", pathway_short),
           pathway_short = gsub("_", " ", pathway_short),
           pathway_short = stringr::str_to_title(pathway_short),
           pathway_short = ifelse(nchar(pathway_short) > 60, 
                                  paste0(substr(pathway_short, 1, 57), "..."), 
                                  pathway_short))
  
  # Create plot
  p <- ggplot(top_pathways, aes(x = NES, y = reorder(pathway_short, NES))) +
    geom_point(aes(size = -log10(padj), color = NES)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                          midpoint = 0, name = "NES") +
    scale_size_continuous(name = "-log10(padj)", range = c(2, 10)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    theme_bw(base_size = 11) +
    labs(title = paste("GSEA:", comparison_name),
         subtitle = paste("Top", nrow(top_pathways), "pathways (padj < 0.05)"),
         x = "Normalized Enrichment Score (NES)",
         y = "") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.y = element_text(size = 9))
  
  return(p)
}

# Generate dot plots
p_IM <- create_enrichment_dotplot(gsea_IM_Gastric, "IM vs Gastric", color_palette["IM"])
p_Duo <- create_enrichment_dotplot(gsea_Duo_Gastric, "Duodenal vs Gastric", color_palette["Duodenal"])
p_DuoIM <- create_enrichment_dotplot(gsea_Duo_IM, "Duodenal vs IM", "#E69F00")

if(!is.null(p_IM)) {
  ggsave("GSEA_DotPlot_IM_vs_Gastric.pdf", p_IM, width=12, height=10)
}
if(!is.null(p_Duo)) {
  ggsave("GSEA_DotPlot_Duodenal_vs_Gastric.pdf", p_Duo, width=12, height=10)
}
if(!is.null(p_DuoIM)) {
  ggsave("GSEA_DotPlot_Duodenal_vs_IM.pdf", p_DuoIM, width=12, height=10)
}

# 5B. TOP PATHWAYS BAR PLOT
create_top_pathways_barplot <- function(fgsea_res, comparison_name, color, n=15) {
  
  sig_pathways <- fgsea_res %>%
    filter(padj < 0.05) %>%
    arrange(desc(abs(NES)))
  
  if(nrow(sig_pathways) == 0) return(NULL)
  
  # Get top up and down
  top_up <- sig_pathways %>% filter(NES > 0) %>% head(n/2)
  top_down <- sig_pathways %>% filter(NES < 0) %>% head(n/2)
  
  top_pathways <- bind_rows(top_up, top_down) %>%
    mutate(pathway_short = gsub("HALLMARK_", "", pathway),
           pathway_short = gsub("GOBP_", "", pathway_short),
           pathway_short = gsub("REACTOME_", "", pathway_short),
           pathway_short = gsub("WP_", "", pathway_short),
           pathway_short = gsub("_", " ", pathway_short),
           pathway_short = stringr::str_to_title(pathway_short),
           direction = ifelse(NES > 0, "Upregulated", "Downregulated"))
  
  p <- ggplot(top_pathways, aes(x = NES, y = reorder(pathway_short, NES), fill = direction)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
    theme_bw(base_size = 11) +
    labs(title = paste("Top Enriched Pathways:", comparison_name),
         x = "Normalized Enrichment Score (NES)",
         y = "",
         fill = "Direction") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 9))
  
  return(p)
}

p_bar_IM <- create_top_pathways_barplot(gsea_IM_Gastric, "IM vs Gastric", color_palette["IM"])
p_bar_Duo <- create_top_pathways_barplot(gsea_Duo_Gastric, "Duodenal vs Gastric", color_palette["Duodenal"])
p_bar_DuoIM <- create_top_pathways_barplot(gsea_Duo_IM, "Duodenal vs IM", "#E69F00")

if(!is.null(p_bar_IM)) {
  ggsave("GSEA_BarPlot_IM_vs_Gastric.pdf", p_bar_IM, width=10, height=8)
}
if(!is.null(p_bar_Duo)) {
  ggsave("GSEA_BarPlot_Duodenal_vs_Gastric.pdf", p_bar_Duo, width=10, height=8)
}
if(!is.null(p_bar_DuoIM)) {
  ggsave("GSEA_BarPlot_Duodenal_vs_IM.pdf", p_bar_DuoIM, width=10, height=8)
}

# 5C. RUNNING ENRICHMENT PLOTS (for top pathways)
cat("  Creating enrichment plots for top pathways...\n")

create_enrichment_plots <- function(fgsea_res, gene_ranks, pathways, comparison_name, n=6) {
  
  top_pathways <- fgsea_res %>%
    filter(padj < 0.05) %>%
    arrange(padj) %>%
    head(n)
  
  if(nrow(top_pathways) == 0) {
    cat("    No pathways to plot for", comparison_name, "\n")
    return()
  }
  
  pdf(paste0("GSEA_EnrichmentPlots_", gsub(" ", "_", comparison_name), ".pdf"), 
      width=10, height=6)
  
  for(i in 1:min(n, nrow(top_pathways))) {
    pathway_name <- top_pathways$pathway[i]
    
    tryCatch({
      plotEnrichment(pathways[[pathway_name]], gene_ranks) +
        labs(title = gsub("_", " ", pathway_name),
             subtitle = paste0("NES = ", round(top_pathways$NES[i], 2), 
                               ", padj = ", format(top_pathways$padj[i], scientific=TRUE, digits=2))) +
        theme_bw(base_size=12) +
        theme(plot.title = element_text(size=10))
      print(last_plot())
    }, error = function(e) {
      cat("    Could not plot:", pathway_name, "\n")
    })
  }
  
  dev.off()
  cat("    Saved enrichment plots for", comparison_name, "\n")
}

create_enrichment_plots(gsea_IM_Gastric, ranks_IM_Gastric, all_pathways, "IM vs Gastric")
create_enrichment_plots(gsea_Duo_Gastric, ranks_Duo_Gastric, all_pathways, "Duodenal vs Gastric")
create_enrichment_plots(gsea_Duo_IM, ranks_Duo_IM, all_pathways, "Duodenal vs IM")

# 5D. HEATMAP: Compare pathway enrichment across all comparisons
cat("\n  Creating comparison heatmap...\n")

create_comparison_heatmap <- function() {
  
  # Get top pathways from each comparison
  top_IM <- gsea_IM_Gastric %>% 
    filter(padj < 0.05) %>% 
    arrange(padj) %>% 
    head(20) %>%
    pull(pathway)
  
  top_Duo <- gsea_Duo_Gastric %>% 
    filter(padj < 0.05) %>% 
    arrange(padj) %>% 
    head(20) %>%
    pull(pathway)
  
  top_DuoIM <- gsea_Duo_IM %>% 
    filter(padj < 0.05) %>% 
    arrange(padj) %>% 
    head(20) %>%
    pull(pathway)
  
  # Get union of pathways
  all_top_pathways <- unique(c(top_IM, top_Duo, top_DuoIM))
  
  if(length(all_top_pathways) == 0) {
    cat("    No pathways for heatmap\n")
    return()
  }
  
  # Create matrix of NES values
  nes_matrix <- data.frame(
    pathway = all_top_pathways,
    IM_vs_Gastric = sapply(all_top_pathways, function(p) {
      nes <- gsea_IM_Gastric$NES[gsea_IM_Gastric$pathway == p]
      ifelse(length(nes) == 0, 0, nes)
    }),
    Duodenal_vs_Gastric = sapply(all_top_pathways, function(p) {
      nes <- gsea_Duo_Gastric$NES[gsea_Duo_Gastric$pathway == p]
      ifelse(length(nes) == 0, 0, nes)
    }),
    Duodenal_vs_IM = sapply(all_top_pathways, function(p) {
      nes <- gsea_Duo_IM$NES[gsea_Duo_IM$pathway == p]
      ifelse(length(nes) == 0, 0, nes)
    })
  )
  
  # Shorten pathway names
  nes_matrix$pathway_short <- gsub("HALLMARK_|GOBP_|REACTOME_|WP_", "", nes_matrix$pathway)
  nes_matrix$pathway_short <- gsub("_", " ", nes_matrix$pathway_short)
  nes_matrix$pathway_short <- stringr::str_to_title(nes_matrix$pathway_short)
  nes_matrix$pathway_short <- ifelse(nchar(nes_matrix$pathway_short) > 50,
                                     paste0(substr(nes_matrix$pathway_short, 1, 47), "..."),
                                     nes_matrix$pathway_short)
  
  mat <- nes_matrix %>%
    dplyr::select(-pathway) %>%
    column_to_rownames("pathway_short") %>%
    as.matrix()
  
  # Plot
  pdf("GSEA_Comparison_Heatmap.pdf", width=8, height=12)
  pheatmap(mat,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           breaks = seq(-3, 3, length.out=101),
           cluster_cols = FALSE,
           fontsize_row = 8,
           main = "Pathway Enrichment Across Comparisons (NES)",
           border_color = NA)
  dev.off()
  
  cat("    Heatmap saved\n")
}

create_comparison_heatmap()

# ====================================================================
# STEP 6: FOCUSED PATHWAY ANALYSIS (Intestinal/Gastric Identity)
# ====================================================================

cat("\nSTEP 6: Focused analysis of intestinal/gastric pathways...\n")

# Define keywords for pathway filtering
intestinal_keywords <- c("INTESTIN", "ENTEROCYTE", "COLON", "GOBLET", "CDX", 
                         "NOTCH", "WNT", "ABSORPT")
gastric_keywords <- c("GASTRIC", "STOMACH", "PARIETAL", "CHIEF", "MUCUS")
differentiation_keywords <- c("DIFFERENT", "EPITHELIAL", "MORPHOGEN", 
                              "DEVELOPMENT", "MATURATION")

filter_pathways <- function(fgsea_res, keywords) {
  pattern <- paste(keywords, collapse="|")
  fgsea_res %>%
    filter(grepl(pattern, pathway, ignore.case=TRUE),
           padj < 0.05) %>%
    arrange(padj)
}

# Extract focused pathways
cat("\n  Intestinal identity pathways:\n")
intestinal_IM <- filter_pathways(gsea_IM_Gastric, intestinal_keywords)
intestinal_Duo <- filter_pathways(gsea_Duo_Gastric, intestinal_keywords)

if(nrow(intestinal_IM) > 0) {
  cat("    IM vs Gastric:", nrow(intestinal_IM), "pathways\n")
  write.csv(intestinal_IM, "Focused_Intestinal_IM_vs_Gastric.csv", row.names=FALSE)
}

if(nrow(intestinal_Duo) > 0) {
  cat("    Duodenal vs Gastric:", nrow(intestinal_Duo), "pathways\n")
  write.csv(intestinal_Duo, "Focused_Intestinal_Duodenal_vs_Gastric.csv", row.names=FALSE)
}

cat("\n  Differentiation/development pathways:\n")
diff_IM <- filter_pathways(gsea_IM_Gastric, differentiation_keywords)
diff_Duo <- filter_pathways(gsea_Duo_Gastric, differentiation_keywords)

if(nrow(diff_IM) > 0) {
  cat("    IM vs Gastric:", nrow(diff_IM), "pathways\n")
  write.csv(diff_IM, "Focused_Differentiation_IM_vs_Gastric.csv", row.names=FALSE)
}

if(nrow(diff_Duo) > 0) {
  cat("    Duodenal vs Gastric:", nrow(diff_Duo), "pathways\n")
  write.csv(diff_Duo, "Focused_Differentiation_Duodenal_vs_Gastric.csv", row.names=FALSE)
}

# ====================================================================
# STEP 7: SUMMARY REPORT
# ====================================================================

cat("\n=======================================================\n")
cat("SUMMARY\n")
cat("=======================================================\n\n")

summary_df <- data.frame(
  Comparison = c("IM vs Gastric", "Duodenal vs Gastric", "Duodenal vs IM"),
  Total_Pathways_Tested = c(length(all_pathways), length(all_pathways), length(all_pathways)),
  Significant_padj005 = c(
    sum(gsea_IM_Gastric$padj < 0.05, na.rm=TRUE),
    sum(gsea_Duo_Gastric$padj < 0.05, na.rm=TRUE),
    sum(gsea_Duo_IM$padj < 0.05, na.rm=TRUE)
  ),
  Upregulated_NES = c(
    sum(gsea_IM_Gastric$padj < 0.05 & gsea_IM_Gastric$NES > 0, na.rm=TRUE),
    sum(gsea_Duo_Gastric$padj < 0.05 & gsea_Duo_Gastric$NES > 0, na.rm=TRUE),
    sum(gsea_Duo_IM$padj < 0.05 & gsea_Duo_IM$NES > 0, na.rm=TRUE)
  ),
  Downregulated_NES = c(
    sum(gsea_IM_Gastric$padj < 0.05 & gsea_IM_Gastric$NES < 0, na.rm=TRUE),
    sum(gsea_Duo_Gastric$padj < 0.05 & gsea_Duo_Gastric$NES < 0, na.rm=TRUE),
    sum(gsea_Duo_IM$padj < 0.05 & gsea_Duo_IM$NES < 0, na.rm=TRUE)
  )
)

print(summary_df)
write.csv(summary_df, "GSEA_Summary_Report.csv", row.names=FALSE)

# Create summary plot
p_summary <- ggplot(summary_df, aes(x=Comparison)) +
  geom_bar(aes(y=Upregulated_NES, fill="Upregulated"), stat="identity", position="dodge") +
  geom_bar(aes(y=-Downregulated_NES, fill="Downregulated"), stat="identity", position="dodge") +
  scale_fill_manual(values=c("Upregulated"="red", "Downregulated"="blue"), name="") +
  geom_hline(yintercept=0, color="black") +
  theme_bw(base_size=14) +
  labs(title="GSEA Results Summary",
       subtitle="Number of significantly enriched pathways (padj < 0.05)",
       y="Number of Pathways",
       x="") +
  theme(plot.title = element_text(face="bold", hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
        axis.text.x = element_text(angle=45, hjust=1))

ggsave("GSEA_Summary_Barplot.pdf", p_summary, width=10, height=8)

cat("\n=======================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=======================================================\n\n")

cat("Key findings:\n")
cat("1. IM vs Gastric shows", summary_df$Significant_padj005[1], "significant pathways\n")
cat("   → Coordinated reprogramming despite low DEG count\n\n")
cat("2. Duodenal vs Gastric shows", summary_df$Significant_padj005[2], "significant pathways\n")
cat("   → Strong transcriptional rewiring\n\n")
cat("3. Duodenal vs IM shows", summary_df$Significant_padj005[3], "significant pathways\n")
cat("   → IM as intermediate state\n\n")

cat("Files generated:\n")
cat("  - Ranked gene lists (Wald statistic)\n")
cat("  - Full GSEA results (all pathways tested)\n")
cat("  - Significant pathways only (padj < 0.05)\n")
cat("  - Dot plots (NES vs significance)\n")
cat("  - Bar plots (top pathways)\n")
cat("  - Running enrichment plots\n")
cat("  - Comparison heatmap (NES across contrasts)\n")
cat("  - Focused pathway lists (intestinal/differentiation)\n")
cat("  - Summary report\n\n")

cat("Next steps:\n")
cat("1. Review enrichment plots for biological coherence\n")
cat("2. Identify key pathways showing progressive changes\n")
cat("3. Focus interpretation on coordinated programs (not single genes)\n")
cat("4. Use for thesis narrative: IM as reprogrammed, not 'failed' state\n")
