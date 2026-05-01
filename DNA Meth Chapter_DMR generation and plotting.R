workdir= "setworkingdirectory"
setwd(workdir)

library(ChAMP)
library(knitr)
library(dplyr)
library(tidyr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(minfi)
library(shinyMethyl)
library(ggplot2)
library(DNAmArray)
library(Gviz)
library(DMRcate)
library(stringr)
library(IlluminaHumanMethylationEPICv2manifest)
library("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
library("ggsignif")
library("gplots")
library("DNAmArray")
library("sva")
library("DNAcopy")
library("impute")
library("wateRmelon")
library("ggfortify")
library("irlba")
library("devtools")
library("fastcluster")
library("ggpubr")
library("qqman")

####ChAMP Load in of data, this section will provide QC plots for both pre and post normalization
Sample_Data <- read.csv("setworkingdirectory/Sample_Data.csv")
myImport <- champ.import("setworkingdirectory/IDATs/All Data", arraytype = "EPICv2")
myLoad <-  champ.filter(beta = myImport$beta,
                        pd=myImport$pd,
                        detP=myImport$detP,
                        beadcount=myImport$beadcount,
                        ProbeCutoff=0.1,
                        arraytype = "EPICv2")
myNorm <- champ.norm(myLoad$beta, method = "BMIQ", arraytype = "EPICv2")


############################################
# CORRECTED: Differential Methylation Analysis
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
# STEP 9: Save results (without Delta_Beta)
############################################

cat("\nSaving initial results...\n")
write.csv(dmp_IM_vs_Gas, "DMPs_IM_vs_Gastric_initial.csv", row.names = FALSE)
write.csv(dmp_Duo_vs_Gas, "DMPs_Duodenal_vs_Gastric_initial.csv", row.names = FALSE)
write.csv(dmp_IM_vs_Duo, "DMPs_IM_vs_Duodenal_initial.csv", row.names = FALSE)


############################################
# STEP 7b: Calculate Beta Differences (Delta Beta)
############################################

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
dmp_IM_vs_Gas$Mean_Gastric <- mean_Gastric[rownames(dmp_IM_vs_Gas)]
dmp_IM_vs_Gas$Mean_IM      <- mean_IM[rownames(dmp_IM_vs_Gas)]
dmp_IM_vs_Gas$Delta_Beta   <- dmp_IM_vs_Gas$Mean_IM - dmp_IM_vs_Gas$Mean_Gastric

# --- Process Duodenal vs Gastric ---
dmp_Duo_vs_Gas$Mean_Gastric <- mean_Gastric[rownames(dmp_Duo_vs_Gas)]
dmp_Duo_vs_Gas$Mean_Duodenal<- mean_Duo[rownames(dmp_Duo_vs_Gas)]
dmp_Duo_vs_Gas$Delta_Beta   <- dmp_Duo_vs_Gas$Mean_Duodenal - dmp_Duo_vs_Gas$Mean_Gastric

# --- Process IM vs Duodenal ---
dmp_IM_vs_Duo$Mean_IM       <- mean_IM[rownames(dmp_IM_vs_Duo)]
dmp_IM_vs_Duo$Mean_Duodenal <- mean_Duo[rownames(dmp_IM_vs_Duo)]
dmp_IM_vs_Duo$Delta_Beta    <- dmp_IM_vs_Duo$Mean_IM - dmp_IM_vs_Duo$Mean_Duodenal

# Save updated files WITH Delta_Beta
cat("\nSaving final results WITH Delta_Beta...\n")
write.csv(dmp_IM_vs_Gas, "DMPs_IM_vs_Gastric.csv", row.names = FALSE)
write.csv(dmp_Duo_vs_Gas, "DMPs_Duodenal_vs_Gastric.csv", row.names = FALSE)
write.csv(dmp_IM_vs_Duo, "DMPs_IM_vs_Duodenal.csv", row.names = FALSE)

cat("\n✅ COMPLETE! All DMP tables saved with Delta_Beta.\n")


############################################
# DMR (Differentially Methylated Region) Analysis
# Using DMRcate package, follows the same design as DMP analysis
############################################

library(DMRcate)
library(limma)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(ggplot2)

############################################
# STEP 1: Prepare data (use from DMP analysis)
############################################

cat("Setting up DMR analysis using existing data...\n")

# Verify objects exist
if (!exists("Mvals_sub")) stop("Run DMP analysis script first!")
if (!exists("design")) stop("Design matrix not found!")
if (!exists("pd")) stop("Sample data not found!")

cat("Data verified. Proceeding with DMR analysis...\n")
cat("M-values dimensions:", dim(Mvals_sub), "\n")
cat("Number of samples:", ncol(Mvals_sub), "\n")
cat("Sample types:", table(pd$Sample_Type), "\n")

############################################
# STEP 2: Annotate CpGs for DMRcate
############################################

cat("\nAnnotating CpGs for DMRcate...\n")

# Create annotation object for DMRcate, needs: chr, pos, gene info, etc.
myAnnotation <- cpg.annotate(
  datatype = "array",
  object = Mvals_sub,
  what = "M",
  arraytype = "EPICv2",
  analysis.type = "differential",
  design = design,
  coef = 2,  
  block = pd$Patient,
  correlation = corfit$consensus,
  fdr = 0.05  
)

cat("CpG annotation complete.\n")

############################################
# STEP 3: Find DMRs - IM vs Gastric
############################################

cat("\n=== Finding DMRs: IM vs Gastric ===\n")

# Find DMRs (lambda = smoothing bandwidth, C = minimum CpGs)
DMRs_IM_vs_Gas <- dmrcate(
  myAnnotation,
  lambda = 1000,      
  C = 2,             
  min.cpgs = 2,       
  pcutoff = 0.05      
)

# Extract results as dataframe
DMRs_IM_vs_Gas_results <- extractRanges(DMRs_IM_vs_Gas, genome = "hg38")

# Convert to data frame for easier manipulation
DMRs_IM_vs_Gas_df <- as.data.frame(DMRs_IM_vs_Gas_results)

cat("Found", nrow(DMRs_IM_vs_Gas_df), "DMRs for IM vs Gastric\n")

# View top DMRs
cat("\nTop 10 DMRs (IM vs Gastric):\n")
print(head(DMRs_IM_vs_Gas_df[order(DMRs_IM_vs_Gas_df$Stouffer, decreasing = TRUE), ], 10))

############################################
# STEP 4: Find DMRs - Duodenal vs Gastric
############################################

cat("\n=== Finding DMRs: Duodenal vs Gastric ===\n")

# Re-annotate for Duodenal comparison
myAnnotation_Duo <- cpg.annotate(
  datatype = "array",
  object = Mvals_sub,
  what = "M",
  arraytype = "EPICv2",
  analysis.type = "differential",
  design = design,
  coef = 3,  
  block = pd$Patient,
  correlation = corfit$consensus,
  fdr = 0.05
)

DMRs_Duo_vs_Gas <- dmrcate(
  myAnnotation_Duo,
  lambda = 1000,
  C = 2,
  min.cpgs = 2,
  pcutoff = 0.05
)

DMRs_Duo_vs_Gas_results <- extractRanges(DMRs_Duo_vs_Gas, genome = "hg38")
DMRs_Duo_vs_Gas_df <- as.data.frame(DMRs_Duo_vs_Gas_results)

cat("Found", nrow(DMRs_Duo_vs_Gas_df), "DMRs for Duodenal vs Gastric\n")

############################################
# STEP 5: Find DMRs - IM vs Duodenal
############################################

cat("\n=== Finding DMRs: IM vs Duodenal ===\n")

# Create contrast for IM vs Duodenal
contrast.matrix <- makeContrasts(
  IM_vs_Duo = Sample_TypeIM - Sample_TypeDuodenal,
  levels = design
)

# Fit the model with contrast
fit_IM_Duo <- lmFit(
  Mvals_sub,
  design,
  block = pd$Patient,
  correlation = corfit$consensus
)

fit_IM_Duo <- contrasts.fit(fit_IM_Duo, contrast.matrix)
fit_IM_Duo <- eBayes(fit_IM_Duo)

# Annotate for this contrast
myAnnotation_IM_Duo <- cpg.annotate(
  datatype = "array",
  object = Mvals_sub,
  what = "M",
  arraytype = "EPICv2",
  analysis.type = "differential",
  design = design,
  contrasts = TRUE,
  cont.matrix = contrast.matrix,
  coef = "IM_vs_Duo",
  block = pd$Patient,
  correlation = corfit$consensus,
  fdr = 0.05
)

DMRs_IM_vs_Duo <- dmrcate(
  myAnnotation_IM_Duo,
  lambda = 1000,
  C = 2,
  min.cpgs = 2,
  pcutoff = 0.05
)

DMRs_IM_vs_Duo_results <- extractRanges(DMRs_IM_vs_Duo, genome = "hg38")
DMRs_IM_vs_Duo_df <- as.data.frame(DMRs_IM_vs_Duo_results)

cat("Found", nrow(DMRs_IM_vs_Duo_df), "DMRs for IM vs Duodenal\n")

############################################
# STEP 6: Summary Statistics
############################################

cat("\n=== DMR SUMMARY ===\n")

create_dmr_summary <- function(dmr_df, comparison_name) {
  cat("\n", comparison_name, ":\n", sep = "")
  cat("  Total DMRs:", nrow(dmr_df), "\n")
  
  if (nrow(dmr_df) > 0) {
    cat("  Median width:", median(dmr_df$width), "bp\n")
    cat("  Mean CpGs per DMR:", round(mean(dmr_df$no.cpgs), 1), "\n")
    cat("  DMRs with 5+ CpGs:", sum(dmr_df$no.cpgs >= 5), "\n")
    cat("  DMRs with 10+ CpGs:", sum(dmr_df$no.cpgs >= 10), "\n")
    
    # Hypermethylated vs Hypomethylated
    if ("meandiff" %in% colnames(dmr_df)) {
      cat("  Hypermethylated:", sum(dmr_df$meandiff > 0), "\n")
      cat("  Hypomethylated:", sum(dmr_df$meandiff < 0), "\n")
    }
  }
}

create_dmr_summary(DMRs_IM_vs_Gas_df, "IM vs Gastric")
create_dmr_summary(DMRs_Duo_vs_Gas_df, "Duodenal vs Gastric")
create_dmr_summary(DMRs_IM_vs_Duo_df, "IM vs Duodenal")

############################################
# STEP 7: Add gene information and clean up
############################################

cat("\nFormatting results...\n")

# Function to clean and format DMR results
format_dmr_results <- function(dmr_df) {
  if (nrow(dmr_df) == 0) return(dmr_df)
  
  # Select and rename key columns
  result <- dmr_df[, c("seqnames", "start", "end", "width", "no.cpgs", 
                       "minfdr", "Stouffer", "maxdiff", "meandiff", 
                       "overlapping.genes")]
  
  colnames(result) <- c("chr", "start", "end", "width_bp", "num_cpgs",
                        "min_fdr", "stouffer_p", "max_diff", "mean_diff",
                        "genes")
  
  # Sort by Stouffer p-value
  result <- result[order(result$stouffer_p, decreasing = TRUE), ]
  
  return(result)
}

DMRs_IM_vs_Gas_clean <- format_dmr_results(DMRs_IM_vs_Gas_df)
DMRs_Duo_vs_Gas_clean <- format_dmr_results(DMRs_Duo_vs_Gas_df)
DMRs_IM_vs_Duo_clean <- format_dmr_results(DMRs_IM_vs_Duo_df)

############################################
# STEP 8: Save results
############################################

cat("\nSaving DMR results...\n")

write.csv(DMRs_IM_vs_Gas_clean, "DMRs_IM_vs_Gastric.csv", row.names = FALSE)
write.csv(DMRs_Duo_vs_Gas_clean, "DMRs_Duodenal_vs_Gastric.csv", row.names = FALSE)
write.csv(DMRs_IM_vs_Duo_clean, "DMRs_IM_vs_Duodenal.csv", row.names = FALSE)

cat("\n✅ DMR analysis complete!\n")
cat("Files saved:\n")
cat("  - DMRs_IM_vs_Gastric.csv\n")
cat("  - DMRs_Duodenal_vs_Gastric.csv\n")
cat("  - DMRs_IM_vs_Duodenal.csv\n")







# Load necessary plotting libraries
library(ggplot2)
library(dplyr)
library(gridExtra)
library(GenomicRanges)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

# ==============================================================================
# HELPER: Filter DMRs for Promoter Regions (TSS)
# ==============================================================================
# This function overlaps  DMR results with the array annotation to 
# keep only regions that contain probes falling in TSS1500 or TSS200.
# ==============================================================================

subset_dmrs_to_promoters <- function(dmr_df) {
  
  # 1. Check if annotation exists, load if not
  if (!exists("annEPICv2")) {
    message("Loading EPICv2 Annotation for filtering...")
    annEPICv2 <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
  }
  
  # 2. Identify Promoter Probes (TSS200 or TSS1500)
  promoter_probes <- annEPICv2[grep("TSS1500|TSS200", annEPICv2$UCSC_RefGene_Group), ]
  
  # 3. Create GRanges objects
  promoter_gr <- makeGRangesFromDataFrame(promoter_probes, 
                                          seqnames.field = "chr", 
                                          start.field = "pos", 
                                          end.field = "pos")
  
  dmr_gr <- makeGRangesFromDataFrame(dmr_df, keep.extra.columns = TRUE)
  
  # 4. Find Overlaps
  overlaps <- findOverlaps(dmr_gr, promoter_gr)
  
  # 5. Extract Gene Names & Filter
  hits_indices <- queryHits(overlaps)
  probe_indices <- subjectHits(overlaps)
  raw_genes <- promoter_probes$UCSC_RefGene_Name[probe_indices]
  
  gene_mapping <- data.frame(dmr_idx = hits_indices, gene = raw_genes, stringsAsFactors = FALSE) %>%
    group_by(dmr_idx) %>%
    summarise(
      Gene_List = paste(unique(unlist(strsplit(gene, ";"))), collapse = ", ")
    )
  
  promoter_dmrs <- dmr_df[unique(hits_indices), ]
  
  gene_vec <- rep(NA, nrow(dmr_df))
  gene_vec[gene_mapping$dmr_idx] <- gene_mapping$Gene_List
  
  promoter_dmrs$Gene <- gene_vec[unique(hits_indices)]
  promoter_dmrs$Gene[promoter_dmrs$Gene == "" | is.na(promoter_dmrs$Gene)] <- "Unknown"
  
  message(paste("Filtered:", nrow(dmr_df), "total DMRs ->", nrow(promoter_dmrs), "Promoter DMRs"))
  
  return(promoter_dmrs)
}

# ==============================================================================
# FUNCTION: Create Cascading DMR Bar Chart
# ==============================================================================

plot_cascading_dmrs <- function(dmr_df, title_text) {
  
  if(is.null(dmr_df) || nrow(dmr_df) == 0) {
    warning(paste("No DMRs found for", title_text))
    return(NULL)
  }
  
  # Prepare Data
  plot_data <- dmr_df %>%
    mutate(
      GeneDisplay = ifelse(nchar(Gene) > 25, paste0(substr(Gene, 1, 22), "..."), Gene),
      # Standard text label (no bolding tags)
      LocationLabel = paste0(GeneDisplay, " (", seqnames, ":", start, "-", end, ")"),
      Direction = ifelse(meandiff > 0, "Hypermethylated", "Hypomethylated")
    )
  
  # Filter Top 25/25
  top_hyper <- plot_data %>%
    filter(Direction == "Hypermethylated") %>%
    arrange(desc(abs(meandiff))) %>% 
    head(25)
  
  top_hypo <- plot_data %>%
    filter(Direction == "Hypomethylated") %>%
    arrange(desc(abs(meandiff))) %>% 
    head(25)
  
  final_data <- bind_rows(top_hyper, top_hypo)
  
  # Generate Plot
  p <- ggplot(final_data, aes(x = reorder(LocationLabel, meandiff), y = meandiff, fill = Direction)) +
    geom_bar(stat = "identity", width = 0.7) +
    coord_flip() + 
    scale_fill_manual(values = c("Hypermethylated" = "#D53E4F", "Hypomethylated" = "#3288BD")) +
    labs(
      title = title_text,
      subtitle = "Top 25 Promoter-Associated DMRs (Ranked by Delta Beta)",
      y = "Mean Methylation Difference (Delta Beta)",
      x = "Gene & Genomic Region",
      fill = "Status"
    ) +
    theme_bw() +
    theme(
      # Standard element_text avoids ggtext/xfun dependency
      axis.text.y = element_text(size = 8),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      legend.position = "bottom"
    )
  
  return(p)
}

# ==============================================================================
# EXECUTION
# ==============================================================================

# 1. IM vs Gastric
cat("Processing IM vs Gastric...\n")
promoter_DMRs_IM_Gas <- subset_dmrs_to_promoters(DMRs_IM_vs_Gas_df)
p1 <- plot_cascading_dmrs(promoter_DMRs_IM_Gas, "Promoter DMRs: IM vs Gastric")
ggsave("Plot_Promoter_DMR_IM_vs_Gastric.png", plot = p1, width = 10, height = 12, dpi = 300)

# 2. Duodenal vs Gastric
cat("Processing Duodenal vs Gastric...\n")
promoter_DMRs_Duo_Gas <- subset_dmrs_to_promoters(DMRs_Duo_vs_Gas_df)
p2 <- plot_cascading_dmrs(promoter_DMRs_Duo_Gas, "Promoter DMRs: Duodenal vs Gastric")
ggsave("Plot_Promoter_DMR_Duodenal_vs_Gastric.png", plot = p2, width = 10, height = 12, dpi = 300)

# 3. IM vs Duodenal
cat("Processing IM vs Duodenal...\n")
promoter_DMRs_IM_Duo <- subset_dmrs_to_promoters(DMRs_IM_vs_Duo_df)
p3 <- plot_cascading_dmrs(promoter_DMRs_IM_Duo, "Promoter DMRs: IM vs Duodenal")
ggsave("Plot_Promoter_DMR_IM_vs_Duodenal.png", plot = p3, width = 10, height = 12, dpi = 300)

cat("\n✅ Promoter-specific plots saved to working directory.\n")






############################################
# DMR-CENTRIC GENOMIC CONTEXT PLOTS
# Counts DMRs (not probes) per region
# Generates consolidated PDFs with Percent Labels
############################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(stringr)

# ==============================================================================
# 1. CLASSIFICATION FUNCTIONS
# ==============================================================================

classify_dmrs <- function(dmr_df, annotation_df, comparison_name) {
  
  if(nrow(dmr_df) == 0) return(NULL)
  
  ann_gr <- makeGRangesFromDataFrame(annotation_df, 
                                     seqnames.field = "chr", 
                                     start.field = "pos", 
                                     end.field = "pos",
                                     keep.extra.columns = TRUE)
  
  dmr_gr <- makeGRangesFromDataFrame(dmr_df, keep.extra.columns = TRUE)
  overlaps <- findOverlaps(dmr_gr, ann_gr)
  
  gene_col <- if("GencodeV41_Group" %in% colnames(annotation_df)) "GencodeV41_Group" else "UCSC_RefGene_Group"
  
  mapping <- data.frame(
    DMR_Index = queryHits(overlaps),
    Gene_Group = annotation_df[[gene_col]][subjectHits(overlaps)],
    Island_Rel = annotation_df$Relation_to_Island[subjectHits(overlaps)]
  )
  
  resolve_gene <- function(groups) {
    all_tags <- paste(groups, collapse = ";")
    splits <- unique(unlist(strsplit(all_tags, ";")))
    if ("TSS200" %in% splits) return("TSS200")
    if ("exon_1" %in% splits || "1stExon" %in% splits) return("1stExon")
    if ("TSS1500" %in% splits) return("TSS1500")
    if ("5UTR" %in% splits || "5'UTR" %in% splits) return("5'UTR")
    if ("3UTR" %in% splits || "3'UTR" %in% splits) return("3'UTR")
    if (any(grepl("exon_|Body|Exon", splits))) return("Body")
    return("IGR")
  }
  
  resolve_island <- function(groups) {
    all_tags <- paste(groups, collapse = ";")
    if (grepl("Island", all_tags)) return("Island")
    if (grepl("Shore", all_tags)) return("Shore")
    if (grepl("Shelf", all_tags)) return("Shelf")
    return("OpenSea")
  }
  
  dmr_classifications <- mapping %>%
    group_by(DMR_Index) %>%
    summarise(
      Gene_Context = resolve_gene(Gene_Group),
      Island_Context = resolve_island(Island_Rel)
    )
  
  dmr_df$DMR_Index <- 1:nrow(dmr_df)
  final_df <- merge(dmr_df, dmr_classifications, by = "DMR_Index", all.x = TRUE)
  final_df$Gene_Context[is.na(final_df$Gene_Context)] <- "IGR"
  final_df$Island_Context[is.na(final_df$Island_Context)] <- "OpenSea"
  final_df$Comparison <- comparison_name
  final_df$Status <- ifelse(final_df$meandiff > 0, "Hypermethylated", "Hypomethylated")
  
  return(final_df)
}

# ==============================================================================
# 2. PROCESS DATA
# ==============================================================================

cat("Processing and Classifying DMRs...\n")
res_1 <- classify_dmrs(DMRs_IM_vs_Gas_df, annEPICv2Sub, "IM vs Gastric")
res_2 <- classify_dmrs(DMRs_Duo_vs_Gas_df, annEPICv2Sub, "Duodenal vs Gastric")
res_3 <- classify_dmrs(DMRs_IM_vs_Duo_df, annEPICv2Sub, "IM vs Duodenal")

all_dmrs <- rbind(res_1, res_2, res_3)

all_dmrs$Gene_Context <- factor(all_dmrs$Gene_Context, 
                                levels = c("TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "3'UTR", "IGR"))

all_dmrs$Island_Context <- factor(all_dmrs$Island_Context, 
                                  levels = c("Island", "Shore", "Shelf", "OpenSea"))

# ==============================================================================
# 3. PLOTTING FUNCTION
# ==============================================================================

generate_consolidated_pdf <- function(data, context_col, filename, title_suffix) {
  
  cat(paste("Generating PDF:", filename, "...\n"))
  pdf(filename, width = 8, height = 7)
  
  comparisons <- unique(data$Comparison)
  
  for (comp in comparisons) {
    
    sub_data <- data %>% 
      filter(Comparison == comp) %>%
      group_by(Status, !!sym(context_col)) %>%
      summarise(Count = n(), .groups = 'drop') %>%
      mutate(Percentage = (Count / sum(Count)) * 100) %>%
      complete(Status, !!sym(context_col), fill = list(Count = 0, Percentage = 0))
    
    y_limit <- max(sub_data$Percentage, na.rm = TRUE) * 1.2
    if(y_limit < 15) y_limit <- 15
    
    p <- ggplot(sub_data, aes(x = !!sym(context_col), y = Percentage, fill = Status)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black") +
      # Add percentage labels above bars
      geom_text(aes(label = sprintf("%.1f%%", Percentage)), 
                position = position_dodge(width = 0.8), 
                vjust = -0.5, 
                size = 3.2, 
                fontface = "bold") +
      scale_fill_manual(values = c("Hypermethylated" = "#D53E4F", "Hypomethylated" = "#3288BD")) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, y_limit)) +
      scale_x_discrete(drop = FALSE) +
      theme_classic() +
      labs(title = comp,
           subtitle = paste("DMR Distribution:", title_suffix),
           x = "Genomic Context",
           y = "Percentage of Total DMRs (%)",
           fill = "") +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey30"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 13),
        legend.position = "top",
        panel.grid.major.y = element_line(color = "grey90", size = 0.2)
      )
    
    print(p)
  }
  
  dev.off()
}

# ==============================================================================
# 4. EXECUTE PLOTS
# ==============================================================================

# Create Consolidated PDF for Gene Context
generate_consolidated_pdf(all_dmrs, "Gene_Context", "DMR_Report_Gene_Context.pdf", "Gene Regions")

# Create Consolidated PDF for Island Context
generate_consolidated_pdf(all_dmrs, "Island_Context", "DMR_Report_Island_Context.pdf", "CpG Island Context")




