library(Seurat)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(harmony)

setwd("/Users/callu/OneDrive - University College London/Protocol/scRNA/")

# 1. Load Data
message("Loading Seurat object...")
seu <- readRDS("seu_cytotrace_complete.rds")

head(seu@meta.data)

table(seu$final_annotation)

seu <- RunHarmony(seu, "dataset")
seu <- RunUMAP(seu, reduction = "harmony")

# FORCE consistency: Replace all underscores in sample IDs with dashes
seu$orig.ident <- gsub("_", "-", seu$orig.ident)

# Create a clean mapping table from the cell-level metadata
meta_map <- seu@meta.data %>%
  select(orig.ident, patient, tissue_group) %>%
  distinct() %>%
  filter(!is.na(patient)) 

# 2. Aggregation
message("Aggregating counts...")
# We use a unique separator that won't be in the lineage names
seu$pseudo_group <- paste0(seu$final_annotation, "---", seu$orig.ident)

pseudo <- AggregateExpression(seu, 
                              group.by = "pseudo_group", 
                              assays = "RNA", 
                              return.seurat = FALSE)
cts <- pseudo$RNA

# 3. Build Metadata and CHECK for NAs
message("Building metadata...")
col_names <- colnames(cts)

coldata <- data.frame(id = col_names) %>%
  mutate(
    # Everything before '---' is lineage, everything after is sample
    final_annotation = sub("---.*", "", id),
    orig.ident = sub(".*---", "", id)
  )

# Join with our Source of Truth
coldata <- coldata %>%
  left_join(meta_map, by = "orig.ident")

# --- DEBUGGING CHECK ---
if (any(is.na(coldata$patient))) {
  message("CRITICAL ERROR: Patient IDs still contain NAs!")
  print(head(coldata[is.na(coldata$patient), ]))
  stop("Fix the mapping above before continuing.")
}

# Final factor setup
coldata$tissue_group <- factor(coldata$tissue_group, levels = c("GAS", "DUO", "IM"))
coldata$patient <- factor(coldata$patient)
rownames(coldata) <- coldata$id

# 4. The Loop (Same logic as before, now with clean data)
dir.create("DE_results", showWarnings = FALSE)
lineages <- unique(coldata$final_annotation)

for (l in lineages) {
  message(paste("Processing:", l))
  l_info <- coldata[coldata$final_annotation == l, ]
  l_cts <- cts[, rownames(l_info), drop = FALSE]
  
  if (ncol(l_cts) < 4) {
    message(paste("   Skipping", l, "- insufficient samples."))
    next
  }
  
  try({
    dds <- DESeqDataSetFromMatrix(countData = l_cts, colData = l_info, design = ~ patient + tissue_group)
    dds <- DESeq(dds)
    
    available_groups <- unique(as.character(l_info$tissue_group))
    
    # Run the Comparisons
    if (all(c("IM", "GAS") %in% available_groups)) {
      res <- results(dds, contrast=c("tissue_group", "IM", "GAS"))
      write.csv(as.data.frame(res), paste0("DE_results/", l, "_IM_vs_GAS.csv"))
    }
    
    # ... (Repeat for other contrasts as in the previous block) ...
    
    message(paste("   Successfully processed:", l))
  })
}





#### Pseudobulk DEG 
# =============================================================================
# scRNA-seq Analysis: Stem/Progenitor vs Central Mixed
# DESeq2 pseudobulk + tissue_group heatmap + GSEA
# =============================================================================

# --- 0. LIBRARIES -------------------------------------------------------------
library(Seurat)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(fgsea)
library(msigdbr)
library(tibble)
library(tidyr)
library(RColorBrewer)

# =============================================================================
# SECTION 1: CELL SELECTION VIA UMAP
# =============================================================================

# Plot UMAP coloured by tissue_group — use this to orient yourself
p <- DimPlot(
  seu,
  reduction  = "umap",
  group.by   = "tissue_group",
  cols       = c("GAS" = "#D6604D", "IM" = "#9970AB", "DUO" = "#4393C3"),
  pt.size    = 0.5
) + ggtitle("UMAP — select populations using CellSelector()")

print(p)

# ------ Manual selection -------
# A window will open. Draw a polygon around each population, then click 'Done'.
message(">>> Select the CENTRAL MIXED population on the UMAP, then click Done.")
central_cells <- CellSelector(p)

message(">>> Now select the STEM/PROGENITOR population, then click Done.")
stem_cells <- CellSelector(p)

# ------ Tag cells in metadata -------
seu$umap_region <- "Other"
seu$umap_region[central_cells] <- "Central_Mixed"
seu$umap_region[stem_cells]    <- "Stem_Progenitor"

# Quick sanity check
message(sprintf(
  "Cells selected — Central Mixed: %d | Stem/Progenitor: %d | Other: %d",
  sum(seu$umap_region == "Central_Mixed"),
  sum(seu$umap_region == "Stem_Progenitor"),
  sum(seu$umap_region == "Other")
))

# =============================================================================
# SECTION 2: PSEUDOBULK DESeq2  (Stem/Progenitor vs Central Mixed)
# =============================================================================

# Subset to the two populations of interest
seu_sub <- subset(seu, umap_region %in% c("Central_Mixed", "Stem_Progenitor"))

# Pseudobulk: sum counts per patient × umap_region
bulk <- AggregateExpression(
  seu_sub,
  group.by      = c("patient", "umap_region"),
  assays        = "RNA",
  return.seurat = TRUE
)

# Clean up metadata factors
bulk$umap_region <- factor(bulk$umap_region, levels = c("Central_Mixed", "Stem_Progenitor"))
bulk$patient     <- as.factor(bulk$patient)

# Propagate tissue_group onto pseudobulk samples
# (majority tissue_group per patient × region combination)
tissue_map <- seu_sub@meta.data %>%
  group_by(patient, umap_region) %>%
  count(tissue_group) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(sample_id = paste(patient, umap_region, sep = "_"))

bulk$tissue_group <- tissue_map$tissue_group[
  match(paste(bulk$patient, bulk$umap_region, sep = "_"), tissue_map$sample_id)
]
bulk$tissue_group <- factor(bulk$tissue_group, levels = c("GAS", "IM", "DUO"))

# Verify no NAs
stopifnot(sum(is.na(bulk$umap_region)) == 0)
message("Pseudobulk sample table:")
print(bulk@meta.data[, c("patient", "umap_region", "tissue_group")])

# Build DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = LayerData(bulk, assay = "RNA", layer = "counts"),
  colData   = bulk@meta.data,
  design    = ~ patient + umap_region   
)

# Filter lowly expressed genes (keep genes with ≥10 counts in ≥2 samples)
keep <- rowSums(counts(dds) >= 10) >= 2
dds  <- dds[keep, ]
message(sprintf("Genes retained after low-count filter: %d", sum(keep)))

# Run DESeq2
dds <- DESeq(dds)

# Extract results: positive LFC = higher in Stem_Progenitor vs Central_Mixed
res <- results(
  dds,
  contrast = c("umap_region", "Stem_Progenitor", "Central_Mixed"),
  alpha    = 0.05
)

res_df <- as.data.frame(res) %>%
  filter(!is.na(padj)) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

message(sprintf(
  "DEGs (padj < 0.05): %d up in Stem, %d up in Central",
  sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE),
  sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE)
))

# Save DEG table
write.csv(res_df, "DEGs_Stem_vs_Central.csv", row.names = FALSE)

# Volcano plot
res_df <- res_df %>%
  mutate(
    sig       = padj < 0.05 & abs(log2FoldChange) > 1,
    direction = case_when(
      sig & log2FoldChange > 0 ~ "Up in Stem",
      sig & log2FoldChange < 0 ~ "Up in Central",
      TRUE                     ~ "NS"
    )
  )

top_labels <- res_df %>%
  filter(sig) %>%
  slice_max(abs(log2FoldChange), n = 20)

volcano_plot <- ggplot(res_df, aes(log2FoldChange, -log10(padj), colour = direction)) +
  geom_point(alpha = 0.6, size = 1.2) +
  ggrepel::geom_text_repel(
    data        = top_labels,
    aes(label   = gene),
    size        = 3,
    max.overlaps = 20
  ) +
  scale_colour_manual(values = c(
    "Up in Stem"    = "#E64B35",
    "Up in Central" = "#4DBBD5",
    "NS"            = "grey70"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  labs(
    title    = "Stem/Progenitor vs Central Mixed (DESeq2 pseudobulk)",
    subtitle = "Positive LFC = higher in Stem/Progenitor",
    x        = "log2 Fold Change",
    y        = "-log10(padj)",
    colour   = NULL
  ) +
  theme_classic(base_size = 13)

print(volcano_plot)
ggsave("volcano_Stem_vs_Central.pdf", volcano_plot, width = 8, height = 6)

# =============================================================================
# SECTION 3: HEATMAP  (top DEGs × tissue_group × population)
# =============================================================================

# Select top DEGs by significance + effect size
n_genes <- 50   

top_genes <- res_df %>%
  filter(padj < 0.05) %>%
  slice_max(abs(log2FoldChange), n = n_genes) %>%
  pull(gene)

# Get VST-normalised counts for visualisation
vst_mat <- assay(vst(dds, blind = FALSE))
hm_mat  <- vst_mat[top_genes, , drop = FALSE]

# Build column annotation: tissue_group and umap_region per pseudobulk sample
col_anno_df <- bulk@meta.data[colnames(hm_mat), c("umap_region", "tissue_group")] %>%
  as.data.frame()

anno_colours <- list(
  umap_region  = c("Central_Mixed"   = "#8DA0CB", "Stem_Progenitor" = "#FC8D62"),
  tissue_group = c("GAS" = "#D6604D", "IM"  = "#9970AB",  "DUO" = "#4393C3")
)

col_annotation <- HeatmapAnnotation(
  df             = col_anno_df,
  col            = anno_colours,
  annotation_name_side = "left"
)

# Row annotation: direction of change
row_anno_df <- data.frame(
  Direction = ifelse(
    res_df$log2FoldChange[match(top_genes, res_df$gene)] > 0,
    "Up in Stem", "Up in Central"
  ),
  row.names = top_genes
)

row_annotation <- rowAnnotation(
  df  = row_anno_df,
  col = list(Direction = c("Up in Stem" = "#E64B35", "Up in Central" = "#4DBBD5"))
)

# Scale rows (z-score per gene)
hm_scaled <- t(scale(t(hm_mat)))

# Colour scale
col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#4575B4", "white", "#D73027")
)

hm <- Heatmap(
  hm_scaled,
  name                  = "z-score",
  col                   = col_fun,
  top_annotation        = col_annotation,
  right_annotation      = row_annotation,
  show_column_names     = TRUE,
  show_row_names        = TRUE,
  row_names_gp          = gpar(fontsize = 8),
  column_names_gp       = gpar(fontsize = 9),
  cluster_columns       = TRUE,
  cluster_rows          = TRUE,
  column_split          = col_anno_df$umap_region,   
  column_title_gp       = gpar(fontsize = 11, fontface = "bold"),
  heatmap_legend_param  = list(title = "z-score"),
  border                = TRUE
)

pdf("heatmap_DEGs_tissue_x_population.pdf", width = 12, height = 10)
draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
message("Heatmap saved.")

# =============================================================================
# SECTION 4: GENE SET ENRICHMENT ANALYSIS (GSEA)
# =============================================================================
message("Running GSEA...")

# 1. Prepare the ranked gene list
# The DESeq2 Wald statistic ('stat') is used for ranking.
# Positive stat = higher in Stem-Progenitor; Negative stat = higher in Central-Mixed
ranked_genes <- res_df %>%
  filter(!is.na(stat)) %>%           
  arrange(desc(stat)) %>%           
  select(gene, stat) %>%
  deframe()                          

# 2. Fetch Human Pathways from MSigDB
# Starting with the 'H' (Hallmark) gene sets for broad biological themes
m_df <- msigdbr(species = "Homo sapiens", category = "H") 
pathways <- split(x = m_df$gene_symbol, f = m_df$gs_name)

# 3. Run fgseaMultilevel
set.seed(42) 
fgsea_res <- fgseaMultilevel(
  pathways = pathways,
  stats    = ranked_genes,
  minSize  = 15,
  maxSize  = 500
)

# 4. Process results for plotting
# Grab the top enriched pathways for both populations
top_pathways <- fgsea_res %>%
  filter(padj < 0.05) %>%
  arrange(desc(NES))

# If there are many pathways, grab the top 10 from each side to keep the plot clean
if(nrow(top_pathways) > 20) {
  top_pathways <- top_pathways %>% 
    slice(c(1:10, (n()-9):n()))
}

top_pathways <- top_pathways %>%
  mutate(
    Direction = ifelse(NES > 0, "Enriched in Stem", "Enriched in Central Mixed"),
    # Clean up the pathway names for aesthetics (remove "HALLMARK_" and replace underscores)
    Pathway = gsub("^HALLMARK_", "", pathway),
    Pathway = gsub("_", " ", Pathway)
  )

# 5. Plot the Normalized Enrichment Scores (NES)
gsea_plot <- ggplot(top_pathways, aes(x = reorder(Pathway, NES), y = NES, fill = Direction)) +
  geom_col(color = "black", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Enriched in Central Mixed" = "#4DBBD5", 
    "Enriched in Stem"          = "#E64B35"
  )) +
  theme_minimal(base_size = 13) +
  labs(
    title    = "GSEA: Stem/Progenitor vs Central Mixed",
    subtitle = "MSigDB Hallmark Pathways (padj < 0.05)",
    x        = NULL,
    y        = "Normalized Enrichment Score (NES)"
  ) +
  theme(
    legend.position = "bottom",
    legend.title    = element_blank(),
    panel.grid.major.y = element_blank() 
  )

print(gsea_plot)
ggsave("GSEA_Stem_vs_Central_Hallmarks.pdf", gsea_plot, width = 9, height = 6)
message("GSEA plot saved.")

# Save the full GSEA results table
fgsea_res_out <- fgsea_res %>% 
  arrange(desc(NES)) %>% 
  select(-leadingEdge) 

write.csv(fgsea_res_out, "GSEA_Results_Stem_vs_Central.csv", row.names = FALSE)









# ============================================================
# Progenitor DE Analysis
# Definition: CytoTRACE2_Potency != "Differentiated"
# Groups: GAS / IM / DUO (via tissue_group, collapsing IM1+IM2)
# Method: pseudo-bulk DESeq2, exploratory (n=2 donors)
# ============================================================

library(Seurat)
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)

# ── Unified colour palette ─────────────────────────────────────
tissue_cols <- c(
  GAS = "#B5546A",   # muted red
  IM  = "#7B6BA8",   # muted purple
  DUO = "#4F86A8"    # muted blue
)

donor_cols <- c(
  "CO-06" = "#888780",
  "CO-07" = "#2C2C2A"
)

# ── 1. Subset to progenitors ──────────────────────────────────
progenitors <- subset(seu, subset = CytoTRACE2_Potency != "Differentiated")

cat("Total progenitor cells:", ncol(progenitors), "\n")
print(table(progenitors$tissue_group))
print(table(progenitors$orig.ident))

# ── 2. Build pseudo-bulk count matrix ─────────────────────────
DefaultAssay(progenitors) <- "RNA"

pb <- AggregateExpression(
  progenitors,
  group.by      = "orig.ident",
  assays        = "RNA",
  return.seurat = FALSE
)$RNA


cat("\nPseudo-bulk matrix dimensions:", nrow(pb), "genes x", ncol(pb), "samples\n")
cat("Column names:", colnames(pb), "\n")

# ── 3. Build matching sample metadata ─────────────────────────
# Parse directly from pb column names (CO-06-GAS etc.)
# to avoid underscore/hyphen mismatch from AggregateExpression
pb_meta <- data.frame(
  orig.ident = colnames(pb)
) %>%
  mutate(
    patient     = sub("^(CO-\\d+)-.*$", "\\1", orig.ident),
    ident_label = sub("^CO-\\d+-(.+)$", "\\1", orig.ident),
    tissue_group = case_when(
      ident_label %in% c("IM1", "IM2") ~ "IM",
      TRUE ~ ident_label
    ),
    tissue_group = factor(tissue_group, levels = c("GAS", "IM", "DUO")),
    patient      = factor(patient)
  )
rownames(pb_meta) <- pb_meta$orig.ident

cat("\nSample metadata:\n")
print(pb_meta[, c("orig.ident", "patient", "tissue_group")])

# ── 4. DE helper function ──────────────────────────────────────
# group_A vs group_B: positive log2FC = higher in group_A
# Tries ~ patient + tissue_group first; falls back to ~ tissue_group
# if patient is confounded with tissue_group (expected with n=2 donors)

run_pb_de <- function(counts_mat, sample_meta,
                      group_A, group_B,
                      lfc_cutoff  = 1.5,
                      padj_cutoff = 0.05) {
  
  keep <- sample_meta$tissue_group %in% c(group_A, group_B)
  mat  <- counts_mat[, keep]
  meta <- sample_meta[keep, , drop = FALSE]
  meta$tissue_group <- factor(meta$tissue_group,
                              levels = c(group_B, group_A))
  
  # Filter lowly expressed genes
  mat <- mat[rowSums(mat) >= 10, ]
  
  # Try model with patient covariate; fall back if rank-deficient
  dds <- tryCatch({
    d <- DESeqDataSetFromMatrix(mat, meta, design = ~ patient + tissue_group)
    DESeq(d, test = "Wald", fitType = "parametric", quiet = TRUE)
  }, error = function(e) {
    message("Patient covariate caused rank deficiency — fitting ~ tissue_group only")
    d <- DESeqDataSetFromMatrix(mat, meta, design = ~ tissue_group)
    DESeq(d, test = "Wald", fitType = "parametric", quiet = TRUE)
  })
  
  # Extract results
  res <- results(dds,
                 contrast = c("tissue_group", group_A, group_B),
                 alpha    = padj_cutoff)
  
  # LFC shrinkage
  res_shrunk <- tryCatch(
    lfcShrink(dds,
              contrast = c("tissue_group", group_A, group_B),
              res      = res,
              type     = "ashr"),
    error = function(e) {
      message("ashr shrinkage failed — returning unshrunk estimates")
      res
    }
  )
  
  as.data.frame(res_shrunk) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(padj) %>%
    mutate(
      comparison = paste0(group_A, "_vs_", group_B),
      sig        = !is.na(padj) &
        padj < padj_cutoff &
        abs(log2FoldChange) > lfc_cutoff,
      direction  = case_when(
        sig & log2FoldChange > 0 ~ paste0("Up in ", group_A),
        sig & log2FoldChange < 0 ~ paste0("Up in ", group_B),
        TRUE ~ "NS"
      )
    )
}

# ── 5. Run all three comparisons ──────────────────────────────
cat("\nRunning IM vs GAS...\n")
res_IM_GAS  <- run_pb_de(pb, pb_meta, group_A = "IM",  group_B = "GAS")

cat("Running IM vs DUO...\n")
res_IM_DUO  <- run_pb_de(pb, pb_meta, group_A = "IM",  group_B = "DUO")

cat("Running DUO vs GAS...\n")
res_DUO_GAS <- run_pb_de(pb, pb_meta, group_A = "DUO", group_B = "GAS")

# ── 6. Summarise results ──────────────────────────────────────
cat("\n── DE Summary ──────────────────────────────────────────\n")
for (res in list(res_IM_GAS, res_IM_DUO, res_DUO_GAS)) {
  comp   <- res$comparison[1]
  groupA <- strsplit(comp, "_vs_")[[1]][1]
  groupB <- strsplit(comp, "_vs_")[[1]][2]
  sig    <- filter(res, sig)
  cat(comp, "→", nrow(sig), "DEGs |",
      sum(sig$direction == paste0("Up in ", groupA)), "up in", groupA, "|",
      sum(sig$direction == paste0("Up in ", groupB)), "up in", groupB, "\n")
}

# ── 7. Save CSVs ──────────────────────────────────────────────
write.csv(res_IM_GAS,  "DE_progenitors_IM_vs_GAS.csv",  row.names = FALSE)
write.csv(res_IM_DUO,  "DE_progenitors_IM_vs_DUO.csv",  row.names = FALSE)
write.csv(res_DUO_GAS, "DE_progenitors_DUO_vs_GAS.csv", row.names = FALSE)

# ── 8. Volcano plot function ───────────────────────────────────
plot_volcano <- function(res_df, top_n = 25) {
  
  comp   <- res_df$comparison[1]
  groupA <- strsplit(comp, "_vs_")[[1]][1]
  groupB <- strsplit(comp, "_vs_")[[1]][2]
  
  cols <- c("NS" = "grey75")
  cols[paste0("Up in ", groupA)] <- tissue_cols[groupA]
  cols[paste0("Up in ", groupB)] <- tissue_cols[groupB]
  
  top_genes <- res_df %>%
    filter(sig) %>%
    slice_min(padj, n = top_n, with_ties = FALSE)
  
  ggplot(res_df, aes(log2FoldChange, -log10(padj), colour = direction)) +
    geom_point(alpha = 0.5, size = 1.2) +
    geom_text_repel(data = top_genes, aes(label = gene),
                    size = 3, max.overlaps = 25, box.padding = 0.4) +
    scale_colour_manual(values = cols) +
    geom_vline(xintercept = c(-1.5, 1.5),
               linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    labs(title    = paste("Progenitors:", gsub("_", " ", comp)),
         subtitle = "Pseudo-bulk DESeq2 | |log2FC| > 1.5 | padj < 0.05 | exploratory",
         x        = paste0("log2FC  (positive = higher in ", groupA, ")"),
         y        = "-log10(adjusted p-value)",
         colour   = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top")
}

p1 <- plot_volcano(res_IM_GAS)
p2 <- plot_volcano(res_IM_DUO)
p3 <- plot_volcano(res_DUO_GAS)

ggsave("volcano_progenitors_IM_vs_GAS.pdf",  p1, width = 8, height = 7)
ggsave("volcano_progenitors_IM_vs_DUO.pdf",  p2, width = 8, height = 7)
ggsave("volcano_progenitors_DUO_vs_GAS.pdf", p3, width = 8, height = 7)

# ── 9. Heatmap of top DEGs across all comparisons ─────────────
# Fit full 3-group DESeq2 model for VST normalisation
dds_full <- DESeqDataSetFromMatrix(
  countData = pb[rowSums(pb) >= 10, ],
  colData   = pb_meta,
  design    = ~ tissue_group
)
dds_full <- DESeq(dds_full, quiet = TRUE)
vst_mat  <- assay(vst(dds_full, blind = FALSE))

# Top 30 sig genes per comparison by padj
top_genes <- bind_rows(res_IM_GAS, res_IM_DUO, res_DUO_GAS) %>%
  filter(sig) %>%
  group_by(comparison) %>%
  slice_min(padj, n = 30, with_ties = FALSE) %>%
  pull(gene) %>%
  unique()

heat_mat    <- vst_mat[intersect(top_genes, rownames(vst_mat)), ]
heat_scaled <- t(scale(t(heat_mat)))

col_anno <- data.frame(
  Tissue = pb_meta$tissue_group,
  Donor  = pb_meta$patient,
  row.names = pb_meta$orig.ident
)

anno_cols <- list(
  Tissue = tissue_cols,
  Donor  = donor_cols
)

pdf("heatmap_progenitors_DEGs.pdf", width = 10, height = 14)
pheatmap(
  heat_scaled,
  annotation_col    = col_anno,
  annotation_colors = anno_cols,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  clustering_method = "ward.D2",
  color             = colorRampPalette(c("#4F86A8", "white", "#B5546A"))(100),
  fontsize_row      = 7,
  fontsize_col      = 9,
  main              = "Top DEGs — progenitor cells (z-scored VST)\nExploratory: n=2 donors"
)
dev.off()








