#This script is an attempt to replicate a previous PhD students in the work;
#I had only plots, with no scripts or raw data, to attempt to replicate this students work.


workdir= "setdirectory"
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
Sample_Data <- read.csv("setdirectory/IDATs/All Data/Sample_Data.csv")
myImport <- champ.import("setdirectory/IDATs/All Data", arraytype = "EPICv2")
myLoad <-  champ.filter(beta = myImport$beta,
                        pd=myImport$pd,
                        detP=myImport$detP,
                        beadcount=myImport$beadcount,
                        ProbeCutoff=0.1,
                        arraytype = "EPICv2")

###########################################
#ChAMP tells you the proportion of failed probes, its a good idea to make a few plots just showing that you've got good quality data
############################################
## Libraries
############################################
library(ggplot2)
library(dplyr)
library(stringr)
library(scales)
############################################
## Colours (AS PROVIDED)
############################################
annotation_colors <- list(
  Sample_Type = c(
    Gastric   = "#c85a5a",
    IM        = "#8e7cc3",
    Duodenal  = "#487cac"
  ),
  Patient = c(
    ID3 = "#e6ab02",
    ID4 = "#66a61e",
    ID6 = "#FFD92F",
    ID7 = "#1B9E77",
    ID8 = "#d16a25",
    ID9 = "#2c7c71"
  )
)
############################################
## Input data
############################################
failed_cpg <- data.frame(
  Sample = c("3Gas","3IM","3Duo","4Gas","4IM","4Duo","6Gas","6IM1","6IM2","6Duo",
             "7Gas","7IM1","7IM2","7Duo","8Gas","8IM","8Duo","9Gas","9IM1","9Duo"),
  FailedFraction = c(0.002121540,0.001982808,0.002280549,0.001454557,0.001452423,
                     0.001761903,0.001347840,0.001437482,0.001697873,0.001915576,
                     0.001711746,0.001596491,0.001536729,0.001403333,0.001327563,
                     0.001425743,0.002750105,0.001692537,0.001785381,0.002813069)
)
############################################
## Tidy metadata
############################################
failed_cpg <- failed_cpg %>%
  mutate(
    Patient = paste0("ID", str_extract(Sample, "^\\d+")),
    Sample_Type = case_when(
      str_detect(Sample, "Gas") ~ "Gastric",
      str_detect(Sample, "Duo") ~ "Duodenal",
      str_detect(Sample, "IM")  ~ "IM"
    ),
    Sample_Type = factor(Sample_Type, levels = c("Gastric","IM","Duodenal")),
    Patient = factor(Patient, levels = paste0("ID", c(3,4,6,7,8,9)))
  )
############################################
## Unified theme
############################################
theme_unified <- function() {
  theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 11),
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
}
############################################
# Individual samples (colored by phenotype, grouped by patient)
############################################
# Create x-axis positions with gaps between patients
failed_cpg_plot <- failed_cpg %>%
  arrange(Patient, Sample_Type) %>%
  mutate(
    patient_num = as.numeric(factor(Patient, levels = paste0("ID", c(3,4,6,7,8,9)))),
    x_pos = row_number() + (patient_num - 1) * 0.5
  )
# Calculate bracket positions for each patient
bracket_data <- failed_cpg_plot %>%
  group_by(Patient) %>%
  summarise(
    x_start = min(x_pos) - 0.4,
    x_end = max(x_pos) + 0.4,
    x_mid = mean(x_pos),
    .groups = "drop"
  ) %>%
  mutate(
    y_pos = -0.00015,
    color = annotation_colors$Patient[Patient]
  )

p_individual <- ggplot(failed_cpg_plot, aes(x = x_pos, y = FailedFraction, fill = Sample_Type)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  # Add brackets for each patient
  geom_segment(
    data = bracket_data,
    aes(x = x_start, xend = x_end, y = y_pos, yend = y_pos, color = Patient),
    inherit.aes = FALSE,
    linewidth = 1.2
  ) +
  geom_segment(
    data = bracket_data,
    aes(x = x_start, xend = x_start, y = y_pos, yend = y_pos + 0.00005, color = Patient),
    inherit.aes = FALSE,
    linewidth = 1.2
  ) +
  geom_segment(
    data = bracket_data,
    aes(x = x_end, xend = x_end, y = y_pos, yend = y_pos + 0.00005, color = Patient),
    inherit.aes = FALSE,
    linewidth = 1.2
  ) +
  geom_text(
    data = bracket_data,
    aes(x = x_mid, y = y_pos - 0.00008, label = Patient),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = annotation_colors$Sample_Type, name = "Sample Type") +
  scale_color_manual(values = annotation_colors$Patient, guide = "none") +
  scale_x_continuous(
    breaks = failed_cpg_plot$x_pos,
    labels = failed_cpg_plot$Sample,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.01),
    limits = c(-0.0003, 0.003),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Failed CpG Probe Fraction per Sample",
    x = "",
    y = "Failed CpG Fraction (%)"
  ) +
  theme_unified() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    plot.margin = margin(t = 10, r = 10, b = 40, l = 10)
  )
############################################
#️ Phenotype average (mean ± SD)
############################################
pheno_summary <- failed_cpg %>%
  group_by(Sample_Type) %>%
  summarise(
    Mean = mean(FailedFraction),
    SD   = sd(FailedFraction),
    .groups = "drop"
  )
p_phenotype <- ggplot(
  pheno_summary,
  aes(x = Sample_Type, y = Mean, fill = Sample_Type)
) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    width = 0.2,
    linewidth = 0.6
  ) +
  scale_fill_manual(values = annotation_colors$Sample_Type, name = "Sample Type") +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.01),
    limits = c(0, 0.003),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Failed CpG Probe Fraction by Phenotype",
    x = "Sample Type",
    y = "Failed CpG Fraction (%)"
  ) +
  theme_unified() +
  theme(legend.position = "none")

############################################
# Patient average (mean ± SD)
############################################
patient_summary <- failed_cpg %>%
  group_by(Patient) %>%
  summarise(
    Mean = mean(FailedFraction),
    SD   = sd(FailedFraction),
    .groups = "drop"
  )

p_patient <- ggplot(
  patient_summary,
  aes(x = Patient, y = Mean, fill = Patient)
) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    width = 0.2,
    linewidth = 0.6
  ) +
  scale_fill_manual(values = annotation_colors$Patient, name = "Patient") +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.01),
    limits = c(0, 0.003),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Failed CpG Probe Fraction by Patient",
    x = "Patient",
    y = "Failed CpG Fraction (%)"
  ) +
  theme_unified() +
  theme(legend.position = "none")

############################################
## Print plots
############################################
print(p_individual)
print(p_phenotype)
print(p_patient)

#Ok now back to normalizing and making QC plots
champ.QC(beta=myLoad$beta, pheno=myLoad$pd$Sample_Type)
myNorm <- champ.norm(myLoad$beta, method = "BMIQ", arraytype = "EPICv2")
write.table(myNorm, "myNorm_samples.txt", row.names = TRUE, sep = "\t")
champ.QC(beta=myNorm, pheno = myLoad$pd$Sample_Type)


###So we need to do SVD analysis to make sure Sample_Type, in this case phenotype,
# is correlating with PC1 or PC2. the ChAMP function makes a terrible heatmap, so find
# code beneath that does make two options.
svd_results <- champ.SVD(beta = myNorm,
          rgSet=NULL,
          pd=myLoad$pd,
          RGEffect=FALSE,
          PDFplot=FALSE,
          Rplot=TRUE,
          resultsDir="./CHAMP_SVDimages/")

############################################
## Libraries
############################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

############################################
## Extract SVD p-values from ChAMP output
############################################
# Assuming your ChAMP.SVD output is stored in a variable
# If you ran: svd_results <- champ.SVD(...)
# Then the p-value matrix is in svd_results
# For this example, I'll use the matrix you showed in the output
svd_pvalues <- matrix(
  c(0.0003299179, 0.98811200, 0.4568361, 0.67151812, 0.5095227, 0.4568361,
    0.0144350610, 0.27704119, 0.4568361, 0.15528393, 0.6837750, 0.4568361,
    0.1225293407, 0.15914617, 0.4568361, 0.01831564, 0.1770076, 0.4568361,
    0.3101084099, 0.38092120, 0.4568361, 0.02924219, 0.2475684, 0.4568361,
    0.0809688735, 0.05662005, 0.4568361, 0.21262733, 0.3842690, 0.4568361),
  nrow = 5, byrow = TRUE
)
colnames(svd_pvalues) <- c("Sample_Type", "Patient", "Sample_Well", "Slide", "Array", "ID")
rownames(svd_pvalues) <- paste0("PC", 1:5)

############################################
## Convert to long format for ggplot
############################################
svd_long <- svd_pvalues %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PC") %>%
  pivot_longer(cols = -PC, names_to = "Variable", values_to = "P_value") %>%
  mutate(
    PC = factor(PC, levels = paste0("PC", 1:5)),
    Variable = factor(Variable, levels = colnames(svd_pvalues)),
    # Create significance categories
    Significance = case_when(
      P_value < 0.001 ~ "p < 0.001",
      P_value < 0.01 ~ "p < 0.01",
      P_value < 0.05 ~ "p < 0.05",
      TRUE ~ "p ≥ 0.05"
    ),
    Significance = factor(Significance, 
                          levels = c("p < 0.001", "p < 0.01", "p < 0.05", "p ≥ 0.05")),
    # Negative log10 for better visualization
    NegLog10P = -log10(P_value)
  )

############################################
## Create custom SVD heatmap
############################################
p_svd <- ggplot(svd_long, aes(x = Variable, y = PC, fill = NegLog10P)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(P_value < 0.05, 
                               sprintf("%.3f", P_value),
                               "")),
            size = 3, fontface = "bold", color = "white") +
  scale_fill_viridis(
    option = "plasma",
    name = "-log10(p-value)",
    limits = c(0, max(svd_long$NegLog10P, na.rm = TRUE))
  ) +
  labs(
    title = "SVD Analysis: Association between PCs and Sample Characteristics",
    subtitle = "P-values shown for significant associations (p < 0.05)",
    x = "Sample Characteristics",
    y = "Principal Component"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )
############################################
## Alternative: Discrete significance levels
############################################
p_svd_discrete <- ggplot(svd_long, aes(x = Variable, y = PC, fill = Significance)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", P_value)),
            size = 3, fontface = "bold") +
  scale_fill_manual(
    values = c("p < 0.001" = "#440154",
               "p < 0.01" = "#31688e",
               "p < 0.05" = "#35b779",
               "p ≥ 0.05" = "#fde724"),
    name = "Significance"
  ) +
  labs(
    title = "SVD Analysis: Association between PCs and Sample Characteristics",
    subtitle = "Heatmap showing p-value significance levels",
    x = "Sample Characteristics",
    y = "Principal Component"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )
############################################
## Print plots
############################################
print(p_svd)
print(p_svd_discrete)

############################################
## Save plots
############################################
ggsave("SVD_heatmap_continuous.pdf", p_svd, width = 10, height = 6)
ggsave("SVD_heatmap_discrete.pdf", p_svd_discrete, width = 10, height = 6)



###########################################
# Now, we want to do some global views of the DNA methylation,
#dataset. This helps just see what the data looks like, where patterns are.
############################################
## Libraries
############################################
library(ggplot2)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(viridis)

############################################
## Colours (AS PROVIDED)
############################################
annotation_colors <- list(
  Sample_Type = c(
    Gastric   = "#c85a5a",
    IM        = "#8e7cc3",
    Duodenal  = "#487cac"
  ),
  Patient = c(
    ID3 = "#e6ab02",
    ID4 = "#66a61e",
    ID6 = "#FFD92F",
    ID7 = "#1B9E77",
    ID8 = "#d16a25",
    ID9 = "#2c7c71"
  )
)

############################################
## ASSUMES:
## myNorm: normalized beta matrix (probes x samples)
## myLoad$pd: sample metadata
############################################

############################################
## Precompute probe variance
############################################
probe_var <- apply(myNorm, 1, var, na.rm = TRUE)

############################################
## Function: PCA workflow
############################################
run_pca <- function(beta_mat, n_probes, file_suffix) {
  
  top_probes <- names(sort(probe_var, decreasing = TRUE)[1:n_probes])
  beta_var <- na.omit(beta_mat[top_probes, ])
  
  pca_res <- prcomp(t(beta_var), scale. = TRUE, center = TRUE)
  var_exp <- summary(pca_res)$importance[2, ] * 100
  
  pca_scores <- as.data.frame(pca_res$x)
  pca_scores$Sample_Type <- factor(myLoad$pd$Sample_Type,
                                   levels = c("Gastric", "IM", "Duodenal"))
  pca_scores$Patient <- factor(myLoad$pd$Patient,
                               levels = paste0("ID", c(3,4,6,7,8,9)))
  pca_scores$Sample <- rownames(pca_scores)
  
  p <- ggplot(pca_scores, aes(PC1, PC2,
                              color = Sample_Type,
                              shape = Patient)) +
    geom_point(size = 4, alpha = 0.8, stroke = 1.2) +
    scale_color_manual(values = annotation_colors$Sample_Type) +
    scale_shape_manual(values = c(15,16,17,18,7,8)) +
    labs(
      title = paste("PCA of Top", n_probes, "Variable Probes"),
      subtitle = "Unsupervised PCA",
      x = paste0("PC1 (", round(var_exp[1],1), "%)"),
      y = paste0("PC2 (", round(var_exp[2],1), "%)")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40")
    )
  
  ggsave(
    filename = paste0("PCA_top", file_suffix, "_probes.pdf"),
    plot = p,
    width = 10,
    height = 7
  )
  
  return(p)
}

############################################
## Function: Heatmap workflow
############################################
run_heatmap <- function(beta_mat, n_probes, file_suffix) {
  
  top_probes <- names(sort(probe_var, decreasing = TRUE)[1:n_probes])
  beta_hm <- na.omit(beta_mat[top_probes, ])
  
  annotation_col <- data.frame(
    Sample_Type = factor(myLoad$pd$Sample_Type,
                         levels = c("Gastric","IM","Duodenal")),
    Patient = factor(myLoad$pd$Patient,
                     levels = paste0("ID", c(3,4,6,7,8,9))),
    row.names = colnames(beta_hm)
  )
  
  pheatmap(
    beta_hm,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    scale = "row",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "complete",
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    show_rownames = FALSE,
    show_colnames = TRUE,
    fontsize_col = 9,
    border_color = NA,
    main = paste("Hierarchical Clustering of Top", n_probes, "Variable Probes"),
    filename = paste0("heatmap_top", file_suffix, "_probes.pdf"),
    width = 10,
    height = 12
  )
  
  # Also display interactively
  pheatmap(
    beta_hm,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    scale = "row",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "complete",
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    show_rownames = FALSE,
    show_colnames = TRUE,
    fontsize_col = 9,
    border_color = NA,
    main = paste("Hierarchical Clustering of Top", n_probes, "Variable Probes")
  )
}

############################################
## Run analyses
############################################
## PCA
pca_1k <- run_pca(myNorm, 1000, "1000")
pca_5k <- run_pca(myNorm, 5000, "5000")
## Heatmaps
run_heatmap(myNorm, 1000, "1000")
run_heatmap(myNorm, 5000, "5000")


############################################
## Additional: Scree plot to see variance
############################################
scree_data <- data.frame(
  PC = factor(paste0("PC", 1:10), levels = paste0("PC", 1:10)),
  Variance = var_explained[1:10]
)

p_scree <- ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_col(fill = "#487cac", color = "black", width = 0.7) +
  geom_text(aes(label = paste0(round(Variance, 1), "%")), 
            vjust = -0.5, size = 3.5, fontface = "bold") +
  labs(
    title = "PCA Scree Plot",
    subtitle = "Variance explained by each principal component",
    x = "Principal Component",
    y = "Variance Explained (%)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

print(p_scree)
ggsave("PCA_scree_plot.pdf", p_scree, width = 8, height = 6)

#############
#So still global methylation, but now i want to look at methylation variability in each phenotype to make violin plots.
#Also want to trace the top 100 variable CpGs across the plots
#############

################################################################################
## ANALYSIS: Methylation Variability (Stochasticity)
## Goal: Visualize "entropy" or disorder (SD) across phenotypes
################################################################################

library(matrixStats)
library(dplyr)
library(tidyr)
library(ggplot2)

# Ensure phenotype order is biologically logical: Gastric -> IM -> Duodenal
# (This order represents the metaplastic transition direction)
pheno_order <- c("Gastric", "IM", "Duodenal")

# Get sample indices for each group
idx_gas <- which(myLoad$pd$Sample_Type == "Gastric")
idx_im  <- which(myLoad$pd$Sample_Type == "IM")
idx_duo <- which(myLoad$pd$Sample_Type == "Duodenal")

# Calculate Standard Deviation (SD) for each probe within each group
# We use matrixStats::rowSds for speed on the large array
sd_gas <- rowSds(myNorm[, idx_gas])
sd_im  <- rowSds(myNorm[, idx_im])
sd_duo <- rowSds(myNorm[, idx_duo])

# Combine into a data frame
var_df <- data.frame(
  ProbeID = rownames(myNorm),
  Gastric = sd_gas,
  IM = sd_im,
  Duodenal = sd_duo
)

# Calculate Global Variability (sum of SDs or SD across all) to find top movers
var_df$Global_Var <- rowSds(myNorm) 

################################
#Identify Probes to Trace
################################

# A. Top 100 Most Variable Globally (The "Grey" lines)
top100_probes <- head(var_df[order(var_df$Global_Var, decreasing = TRUE), "ProbeID"], 100)

# B. IM-Specific Instability (The "Yellow" lines)
# Logic: High SD in IM, but relatively stable (Low SD) in Gastric and Duodenal
# Thresholds: Top 1% in IM, Bottom 50% in others (adjustable)
thresh_im_high <- quantile(var_df$IM, 0.99)
thresh_others_low <- quantile(c(var_df$Gastric, var_df$Duodenal), 0.50)

im_unstable_probes <- var_df %>%
  filter(IM > thresh_im_high & Gastric < thresh_others_low & Duodenal < thresh_others_low) %>%
  pull(ProbeID)

# Limit to top 50 if too many, to prevent overplotting
if(length(im_unstable_probes) > 50) im_unstable_probes <- head(im_unstable_probes, 50)

################################
#Reshape for plotting
################################

# Melt data to long format for ggplot
# We use dplyr::select to avoid conflicts with Bioconductor packages
var_long <- var_df %>%
  dplyr::select(ProbeID, Gastric, IM, Duodenal) %>%
  tidyr::pivot_longer(cols = c("Gastric", "IM", "Duodenal"), 
                      names_to = "Sample_Type", 
                      values_to = "SD") %>%
  dplyr::mutate(Sample_Type = factor(Sample_Type, levels = pheno_order))

# Create subsets for the traces
trace_grey <- var_long %>% dplyr::filter(ProbeID %in% top100_probes)
trace_yellow <- var_long %>% dplyr::filter(ProbeID %in% im_unstable_probes)

#####################################################
#Plot 1: Clean Methylation variability (violin plots)
######################################################

p_var_clean <- ggplot(var_long, aes(x = Sample_Type, y = SD, fill = Sample_Type)) +
  # Violin for distribution shape
  geom_violin(trim = FALSE, alpha = 0.8, color = NA) +
  # Boxplot for median/quartiles
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.8) +
  # Styling
  scale_fill_manual(values = annotation_colors$Sample_Type) +
  labs(title = "Global Methylation Variability Profile",
       subtitle = "Distribution of Probe Standard Deviations per Phenotype",
       y = "Methylation Variability (SD)",
       x = NULL) +
  theme_unified() +
  theme(legend.position = "none") 

######################################################
#Plot 2: Traced Variability (Violin + Lines)
######################################################

# We need a numeric X-axis to draw lines, but categorical labels
# Gastric=1, IM=2, Duodenal=3
p_var_trace <- ggplot(var_long, aes(x = Sample_Type, y = SD)) +
  # Background Violin (Greyed out to highlight traces)
  geom_violin(trim = FALSE, fill = "gray90", color = "gray80", alpha = 0.5) +
  
  # Trace: Top 100 Global (Grey)
  geom_line(data = trace_grey, 
            aes(group = ProbeID), 
            color = "gray40", alpha = 0.3, size = 0.5) +
  
  # Trace: IM Unstable (Yellow/Gold)
  geom_line(data = trace_yellow, 
            aes(group = ProbeID), 
            color = "#ffd700", alpha = 0.9, size = 0.8) +
  
  # Points for the yellow traces (optional, adds emphasis)
  geom_point(data = trace_yellow,
             fill = "#ffd700", color = "black", shape = 21, size = 1.5) +
  
  # Styling
  labs(title = "Trajectory of Methylation Instability",
       subtitle = "Grey: Top 100 variable globally | Gold: High variability specific to IM",
       y = "Methylation Variability (SD)",
       x = NULL) +
  theme_unified() +
  theme(panel.grid.major.x = element_blank())

######################################################
#Save Plots
######################################################

print(p_var_clean)
print(p_var_trace)

ggsave("Methylation_Variability_Violin_Clean.pdf", p_var_clean, width = 6, height = 6)
ggsave("Methylation_Variability_Trace.pdf", p_var_trace, width = 7, height = 6)


######################################################
#Alternative plotting options
######################################################

library(ggpubr) # Required for stat_compare_means

######################################################
#Statistical Comparisons
######################################################

# Define the comparisons you want to show on the plot
my_comparisons <- list( 
  c("Gastric", "IM"), 
  c("IM", "Duodenal"), 
  c("Gastric", "Duodenal") 
)

######################################################
#Plot 1: Clean Variability with Significance Stars
######################################################

p_var_stats <- ggplot(var_long, aes(x = Sample_Type, y = SD, fill = Sample_Type)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.8) +
  # Add T-test significance stars
  stat_compare_means(comparisons = my_comparisons, 
                     method = "t.test", 
                     label = "p.signif",
                     step.increase = 0.1) + 
  scale_fill_manual(values = annotation_colors$Sample_Type) +
  labs(title = "Global Methylation Variability Profile",
       subtitle = "Pairwise T-tests showing differences in probe-wise SD",
       y = "Methylation Variability (SD)",
       x = NULL) +
  theme_unified() +
  theme(legend.position = "none")

######################################################
#Analyze Top 100 Probes
######################################################

# Create a dataframe of just the top 100 probes
top_100_data <- var_df %>% dplyr::filter(ProbeID %in% top100_probes)

# Find which group has the maximum SD for each of these 100 probes
# 1 = Gastric, 2 = IM, 3 = Duodenal (based on columns 2,3,4)
top_100_data$Highest_Group <- colnames(top_100_data[2:4])[apply(top_100_data[2:4], 1, which.max)]

# Count occurrences
im_count <- sum(top_100_data$Highest_Group == "IM")

cat("------------------------------------------------------------\n")
cat("Analysis of Top 100 Most Variable CpGs:\n")
cat(paste0("Number of probes most variable in IM: ", im_count, "/100\n"))
cat("------------------------------------------------------------\n")

######################################################
#FInal Display
######################################################

print(p_var_stats)
# Re-run the trace plot from previous step to see updated visuals
print(p_var_trace)

ggsave("Methylation_Variability_with_Stats.pdf", p_var_stats, width = 7, height = 7)


# Try Log2 of the Variance
summary(log2(rowVars(myNorm)))

# Try Log2 of the SD
summary(log2(rowSds(myNorm)))


################################################################################
## ANALYSIS: M-Value Variability (The "Predecessor" Scale)
################################################################################

# 1. Convert Betas to M-values --------------------------------------------
# We add a small offset to avoid log(0) or log(inf)
m_vals <- log2(myNorm / (1 - myNorm))

# 2. Calculate Log2 Variance per Group ------------------------------------
# We use rowVars on M-values, then take log2
log2_var_gas <- log2(rowVars(m_vals[, idx_gas]))
log2_var_im  <- log2(rowVars(m_vals[, idx_im]))
log2_var_duo <- log2(rowVars(m_vals[, idx_duo]))

# Create Dataframe
var_df_m <- data.frame(
  ProbeID = rownames(myNorm),
  Gastric = log2_var_gas,
  IM = log2_var_im,
  Duodenal = log2_var_duo
)

# Replace -Inf (from zero variance) with the lowest non-infinite value
var_df_m[var_df_m == -Inf] <- min(var_df_m[is.finite(as.matrix(var_df_m))])

# 3. Reshape and Stats ----------------------------------------------------
var_long_m <- var_df_m %>%
  dplyr::select(ProbeID, Gastric, IM, Duodenal) %>%
  tidyr::pivot_longer(cols = -ProbeID, names_to = "Sample_Type", values_to = "Log2Var") %>%
  dplyr::mutate(Sample_Type = factor(Sample_Type, levels = pheno_order))

# Calculate how many of top 100 global variable (M-scale) are max in IM
var_df_m$Global_SD_M <- rowSds(m_vals)
top100_m <- head(var_df_m[order(var_df_m$Global_SD_M, decreasing = TRUE), "ProbeID"], 100)
top_100_data_m <- var_df_m %>% dplyr::filter(ProbeID %in% top100_m)
top_100_data_m$MaxGroup <- colnames(top_100_data_m[2:4])[apply(top_100_data_m[2:4], 1, which.max)]
im_peak_count <- sum(top_100_data_m$MaxGroup == "IM")

# 4. Plot 1: Clean M-Value Variability with Stats -------------------------
p1 <- ggplot(var_long_m, aes(x = Sample_Type, y = Log2Var, fill = Sample_Type)) +
  geom_violin(trim = TRUE, alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
  scale_fill_manual(values = annotation_colors$Sample_Type) +
  labs(title = "Global Epigenetic Dysregulation (M-value scale)",
       subtitle = paste0("Scale matches historical lab metrics (-8 to 4)"),
       y = "Log2(Variance of M-values)", x = NULL) +
  theme_unified() +
  coord_cartesian(ylim = c(-10, 5)) # Forcing the scale to match your predecessor

# 5. Plot 2: Tracing the Top 100 ------------------------------------------
trace_grey_m <- var_long_m %>% dplyr::filter(ProbeID %in% top100_m)

p2 <- ggplot(var_long_m, aes(x = Sample_Type, y = Log2Var)) +
  geom_violin(fill = "gray95", color = NA) +
  geom_line(data = trace_grey_m, aes(group = ProbeID), color = "gray40", alpha = 0.2) +
  # Highlighting probes that peak in IM
  geom_line(data = trace_grey_m %>% dplyr::filter(ProbeID %in% top_100_data_m$ProbeID[top_100_data_m$MaxGroup == "IM"]),
            aes(group = ProbeID), color = "#ffd700", alpha = 0.6, size = 0.7) +
  labs(title = "Trace of Top 100 Variable Probes",
       subtitle = paste(im_peak_count, "of the top 100 variable probes peak in IM"),
       y = "Log2(Variance of M-values)", x = NULL) +
  theme_unified() +
  coord_cartesian(ylim = c(-10, 5))

print(p1)
print(p2)

ggsave("GLobal with M values.pdf", p1, width = 7, height = 7)
ggsave("Global with M_yellow lines.pdf", p2, width = 7, height = 7)



####
#After discussions with what Will's Global variability plot could be,
#it seems the likely suspect is Genome-wide methylation variability was quantified 
#by calculating the median absolute deviation (MAD) of M-values for each CpG across 
#samples within each tissue phenotype. Variability values were log₂-transformed and visualised as violin plots.
############################################
## Libraries
############################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif) # For adding significance bars

############################################
## 1. Data Prep: Calculate Variability with Probe IDs
############################################
# Ensure beta is bounded
beta <- pmin(pmax(myNorm, 1e-6), 1 - 1e-6)
Mvals <- log2(beta / (1 - beta))

# Calculate MAD per group and KEEP PROBE IDs
mv_list <- list()
for (pheno in c("Gastric", "IM", "Duodenal")) {
  samples_pheno <- myLoad$pd$Sample_Type == pheno
  M_sub <- Mvals[, samples_pheno]
  
  # Calculate MAD and keep rownames (ProbeIDs)
  mv <- apply(M_sub, 1, mad, na.rm = TRUE)
  
  mv_list[[pheno]] <- data.frame(
    ProbeID = rownames(M_sub),
    MV = log2(mv + 1e-6),
    Sample_Type = pheno
  )
}

mv_df <- bind_rows(mv_list)

# Force Factor Order: Gastric -> IM -> Duodenal
mv_df$Sample_Type <- factor(mv_df$Sample_Type, levels = c("Gastric", "IM", "Duodenal"))

############################################
## 2. Statistical Testing (Paired T-tests)
############################################
# We need to pivot wide to perform paired tests easily
mv_wide <- mv_df %>%
  pivot_wider(names_from = Sample_Type, values_from = MV)

# Function to run paired t-test and format p-value
get_pval <- function(group1, group2, data) {
  res <- t.test(data[[group1]], data[[group2]], paired = TRUE)
  pval <- res$p.value * 3 # Bonferroni correction (x3 comparisons)
  
  # Format for display
  if(pval < 2.2e-16) return("p < 2.2e-16 ***")
  if(pval < 0.001) return("***")
  if(pval < 0.01) return("**")
  if(pval < 0.05) return("*")
  return("ns")
}

# Run the 3 specific comparisons
p_gas_im  <- get_pval("Gastric", "IM", mv_wide)
p_im_duo  <- get_pval("IM", "Duodenal", mv_wide)
p_gas_duo <- get_pval("Gastric", "Duodenal", mv_wide)

print(paste("Gastric vs IM:", p_gas_im))
print(paste("IM vs Duodenal:", p_im_duo))
print(paste("Gastric vs Duodenal:", p_gas_duo))

############################################
## 3. Prepare Plotting Data (Medians for Line)
############################################
summary_stats <- mv_df %>%
  group_by(Sample_Type) %>%
  summarise(Median_MV = median(MV, na.rm = TRUE))

############################################
## 4. Plot with Stats
############################################
p_final <- ggplot(mv_df, aes(x = Sample_Type, y = MV, fill = Sample_Type)) +
  # A. The Violins
  geom_violin(trim = TRUE, scale = "width", color = "black", linewidth = 0.4, alpha = 0.8) +
  
  # B. The Connecting Trend Line
  geom_line(data = summary_stats, aes(x = Sample_Type, y = Median_MV, group = 1), 
            inherit.aes = FALSE, color = "black", linewidth = 1, linetype = "dashed") +
  
  # C. The Median Points
  geom_point(data = summary_stats, aes(x = Sample_Type, y = Median_MV),
             inherit.aes = FALSE, color = "black", size = 3) +
  
  # D. Statistical Significance Bars
  # Note: y_position must be adjusted based on your actual data max values
  geom_signif(
    comparisons = list(c("Gastric", "IM")),
    annotations = p_gas_im,
    y_position = max(mv_df$MV) + 0.5, tip_length = 0.01, vjust = 0.4
  ) +
  geom_signif(
    comparisons = list(c("IM", "Duodenal")),
    annotations = p_im_duo,
    y_position = max(mv_df$MV) + 1.5, tip_length = 0.01, vjust = 0.4
  ) +
  geom_signif(
    comparisons = list(c("Gastric", "Duodenal")),
    annotations = p_gas_duo,
    y_position = max(mv_df$MV) + 2.5, tip_length = 0.01, vjust = 0.4
  ) +
  
  # Styling
  scale_fill_manual(values = annotation_colors) +
  labs(
    title = "Methylation Variability by Phenotype",
    subtitle = "Paired T-test (Bonferroni adj.) on CpG-wise variability (MAD)",
    x = "",
    y = "MV (log₂ MAD of M-values)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "none",
    axis.text = element_text(size = 11, color = "black"),
    panel.grid.major.x = element_blank()
  ) +
  # Expand Y axis to fit the significance bars
  coord_cartesian(ylim = c(min(mv_df$MV), max(mv_df$MV) + 3.5))

print(p_final)

ggsave("Methylation_Variability_Stats_Final.pdf", p_final, width = 7, height = 7)


############################################
## Libraries
############################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)

############################################
## Colours
############################################
annotation_colors <- c(
  Gastric   = "#c85a5a",
  IM        = "#8e7cc3",
  Duodenal  = "#487cac"
)

############################################
## 1. Calculate CpG-wise variability per phenotype
############################################
# Ensure beta is bounded
beta <- pmin(pmax(myNorm, 1e-6), 1 - 1e-6)
Mvals <- log2(beta / (1 - beta))

# Calculate MAD per CpG per phenotype
mv_list <- list()
for (pheno in c("Gastric", "IM", "Duodenal")) {
  samples_pheno <- myLoad$pd$Sample_Type == pheno
  M_sub <- Mvals[, samples_pheno]
  
  mv <- apply(M_sub, 1, mad, na.rm = TRUE)
  
  mv_list[[pheno]] <- data.frame(
    ProbeID = rownames(M_sub),
    MV = log2(mv + 1e-6),
    Sample_Type = pheno
  )
}

mv_df <- bind_rows(mv_list)
mv_df$Sample_Type <- factor(mv_df$Sample_Type, levels = c("Gastric", "IM", "Duodenal"))

############################################
## 2. Identify top 100 most variable CpGs (across all groups)
############################################
# Calculate overall variability per CpG (mean across phenotypes)
overall_var <- mv_df %>%
  group_by(ProbeID) %>%
  summarise(Mean_MV = mean(MV, na.rm = TRUE)) %>%
  arrange(desc(Mean_MV))

top_100_probes <- overall_var$ProbeID[1:100]

# Extract data for top 100 probes
top_100_data <- mv_df %>%
  filter(ProbeID %in% top_100_probes) %>%
  pivot_wider(names_from = Sample_Type, values_from = MV)

############################################
## 3. Identify IM-specific high variability CpGs (YELLOW)
############################################
# Criteria: Higher variability in IM compared to BOTH Gastric and Duodenal
im_specific <- top_100_data %>%
  filter(IM > Gastric & IM > Duodenal)

n_im_specific <- nrow(im_specific)
im_specific_probes <- im_specific$ProbeID

cat("\n=== IM-SPECIFIC HIGH VARIABILITY CpGs ===\n")
cat("Number of CpGs with highest variability in IM:", n_im_specific, "\n")
cat("(Criteria: IM variability > Gastric AND > Duodenal)\n\n")

############################################
## 4. Prepare data for plotting
############################################
# Add color category
top_100_long <- top_100_data %>%
  pivot_longer(cols = c(Gastric, IM, Duodenal), 
               names_to = "Sample_Type", 
               values_to = "MV") %>%
  mutate(
    Sample_Type = factor(Sample_Type, levels = c("Gastric", "IM", "Duodenal")),
    CpG_Category = ifelse(ProbeID %in% im_specific_probes, "IM-specific", "Other"),
    line_alpha = ifelse(ProbeID %in% im_specific_probes, 0.9, 0.3),
    line_size = ifelse(ProbeID %in% im_specific_probes, 0.8, 0.3)
  )

############################################
## 5. Calculate patient-level stats (for significance bars)
############################################
patient_mv_list <- list()

for (patient_id in unique(myLoad$pd$Patient)) {
  for (pheno in c("Gastric", "IM", "Duodenal")) {
    samples_idx <- myLoad$pd$Patient == patient_id & myLoad$pd$Sample_Type == pheno
    
    if (sum(samples_idx) > 0) {
      M_sub <- Mvals[, samples_idx, drop = FALSE]
      mv_per_cpg <- apply(M_sub, 1, mad, na.rm = TRUE)
      patient_mv <- mean(log2(mv_per_cpg + 1e-6), na.rm = TRUE)
      
      patient_mv_list[[length(patient_mv_list) + 1]] <- data.frame(
        Patient = patient_id,
        Sample_Type = pheno,
        Patient_MV = patient_mv
      )
    }
  }
}

patient_mv_df <- bind_rows(patient_mv_list)
patient_wide <- patient_mv_df %>%
  pivot_wider(names_from = Sample_Type, values_from = Patient_MV)

# Wilcoxon tests with error handling with bnferroni correction
get_wilcox_pval <- function(group1, group2, data) {
  tryCatch({
    res <- wilcox.test(data[[group1]], data[[group2]], paired = TRUE, exact = FALSE)
    pval <- res$p.value * 3  
    
    # Handle NA or invalid p-values
    if(is.na(pval) || is.nan(pval)) return("ns")
    
    if(pval < 0.001) return("***")
    if(pval < 0.01) return("**")
    if(pval < 0.05) return("*")
    return("ns")
  }, error = function(e) {
    warning(paste("Wilcoxon test failed for", group1, "vs", group2, ":", e$message))
    return("ns")
  })
}

p_gas_im  <- get_wilcox_pval("Gastric", "IM", patient_wide)
p_im_duo  <- get_wilcox_pval("IM", "Duodenal", patient_wide)
p_gas_duo <- get_wilcox_pval("Gastric", "Duodenal", patient_wide)

############################################
## 6. Summary stats for horizontal lines
############################################
summary_stats <- mv_df %>%
  group_by(Sample_Type) %>%
  summarise(Median_MV = median(MV, na.rm = TRUE))

############################################
## 7. Create plot
############################################
p_traced <- ggplot(mv_df, aes(x = Sample_Type, y = MV, fill = Sample_Type)) +
  # Background violins
  geom_violin(trim = TRUE, scale = "width", color = "black", 
              linewidth = 0.4, alpha = 0.8) +
  
  # Grey traces for top 100 (other CpGs)
  geom_line(data = filter(top_100_long, CpG_Category == "Other"),
            aes(x = as.numeric(Sample_Type), y = MV, group = ProbeID),
            color = "grey30", alpha = 0.3, linewidth = 0.3, inherit.aes = FALSE) +
  
  # Yellow traces for IM-specific CpGs
  geom_line(data = filter(top_100_long, CpG_Category == "IM-specific"),
            aes(x = as.numeric(Sample_Type), y = MV, group = ProbeID),
            color = "#FFD700", alpha = 0.9, linewidth = 0.8, inherit.aes = FALSE) +
  
  # Horizontal median lines
  geom_segment(data = summary_stats,
               aes(x = as.numeric(Sample_Type) - 0.2, 
                   xend = as.numeric(Sample_Type) + 0.2,
                   y = Median_MV, yend = Median_MV),
               color = "black", linewidth = 1, inherit.aes = FALSE) +
  
  # Statistical significance bars
  geom_signif(
    comparisons = list(c("Gastric", "IM")),
    annotations = p_gas_im,
    y_position = max(mv_df$MV) + 0.5, 
    tip_length = 0.02, 
    textsize = 4,
    vjust = 0.4
  ) +
  geom_signif(
    comparisons = list(c("IM", "Duodenal")),
    annotations = p_im_duo,
    y_position = max(mv_df$MV) + 1.5, 
    tip_length = 0.02,
    textsize = 4,
    vjust = 0.4
  ) +
  geom_signif(
    comparisons = list(c("Gastric", "Duodenal")),
    annotations = p_gas_duo,
    y_position = max(mv_df$MV) + 2.5, 
    tip_length = 0.02,
    textsize = 4,
    vjust = 0.4
  ) +
  
  # Add text annotation for yellow CpG count
  annotate("text", x = 2, y = min(mv_df$MV) - 0.3,
           label = paste0("n = ", n_im_specific, " IM-specific CpGs"),
           color = "#FFD700", fontface = "bold", size = 4) +
  
  # Styling
  scale_fill_manual(values = annotation_colors) +
  labs(
    title = "Methylation Variability with Top 100 Variable CpG Traces",
    subtitle = paste0("Paired Wilcoxon test (n = ", nrow(patient_wide), 
                      " patients). Yellow: IM-specific high variability CpGs"),
    x = "",
    y = "Methylation Variability (log₂ MAD of M-values)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    legend.position = "none",
    axis.text = element_text(size = 11, color = "black", face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(min(mv_df$MV) - 0.8, max(mv_df$MV) + 3.5))

print(p_traced)

############################################
## 8. Save outputs
############################################
ggsave("Methylation_Variability_with_Traces.pdf", p_traced, width = 9, height = 7)

# Save IM-specific CpGs for follow-up
write.csv(im_specific, "IM_Specific_Variable_CpGs.csv", row.names = FALSE)

cat("\n=== SAVED FILES ===\n")
cat("1. Plot: Methylation_Variability_with_Traces.pdf\n")
cat("2. IM-specific CpGs: IM_Specific_Variable_CpGs.csv\n")



############################################
#MAD Plots, plain MAD Plots for across sample, coloured by phenotype, and then just phenotype.
############################################

############################################
## Libraries
############################################
library(ggplot2)
library(dplyr)
library(tidyr)

############################################
## Colours
############################################
annotation_colors <- list(
  Sample_Type = c(
    Gastric   = "#c85a5a",
    IM        = "#8e7cc3",
    Duodenal  = "#487cac"
  ),
  Patient = c(
    ID3 = "#e6ab02",
    ID4 = "#66a61e",
    ID6 = "#FFD92F",
    ID7 = "#1B9E77",
    ID8 = "#d16a25",
    ID9 = "#2c7c71"
  )
)

############################################
## Calculate MAD and SE per sample - BETA VALUES
############################################
# Ensure beta is bounded
beta <- pmin(pmax(myNorm, 1e-6), 1 - 1e-6)

# Function to calculate MAD and its standard error
calc_mad_with_se <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  med <- median(x)
  abs_dev <- abs(x - med)
  mad_val <- median(abs_dev)
  
  # Standard error of MAD: SE = 1.4826 * sd(abs_dev) / sqrt(n)
  # The 1.4826 is a consistency constant for normal distribution
  se_mad <- 1.4826 * sd(abs_dev) / sqrt(n)
  
  return(c(MAD = mad_val, SE = se_mad))
}

# Calculate MAD and SE for each sample - Beta values
mad_beta_results <- apply(beta, 2, calc_mad_with_se)
mad_beta_results <- t(mad_beta_results)

# Create data frame
mad_beta_df <- data.frame(
  Sample = rownames(mad_beta_results),
  MAD = mad_beta_results[, "MAD"],
  SE = mad_beta_results[, "SE"],
  Sample_Type = myLoad$pd$Sample_Type,
  Patient = myLoad$pd$Patient,
  Value_Type = "Beta",
  stringsAsFactors = FALSE
)

############################################
## Calculate MAD and SE per sample - M VALUES
############################################
Mvals <- log2(beta / (1 - beta))

# Calculate MAD and SE for each sample - M values
mad_m_results <- apply(Mvals, 2, calc_mad_with_se)
mad_m_results <- t(mad_m_results)

# Create data frame
mad_m_df <- data.frame(
  Sample = rownames(mad_m_results),
  MAD = mad_m_results[, "MAD"],
  SE = mad_m_results[, "SE"],
  Sample_Type = myLoad$pd$Sample_Type,
  Patient = myLoad$pd$Patient,
  Value_Type = "M-value",
  stringsAsFactors = FALSE
)

############################################
## Combine both
############################################
mad_df <- bind_rows(mad_beta_df, mad_m_df)

# Factor ordering
mad_df$Sample_Type <- factor(mad_df$Sample_Type, 
                             levels = c("Gastric", "IM", "Duodenal"))
mad_df$Patient <- factor(mad_df$Patient, 
                         levels = paste0("ID", c(3, 4, 6, 7, 8, 9)))
mad_df$Value_Type <- factor(mad_df$Value_Type, levels = c("Beta", "M-value"))

# Order samples by Patient then Sample_Type for plotting
sample_order <- mad_beta_df %>%
  arrange(Patient, Sample_Type) %>%
  pull(Sample)

mad_df$Sample <- factor(mad_df$Sample, levels = sample_order)

############################################
## Plot 1: MAD per sample - BETA VALUES
############################################
p_mad_beta <- ggplot(mad_beta_df, aes(x = Sample, y = MAD, fill = Sample_Type)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = MAD - SE, ymax = MAD + SE),
                width = 0.3, linewidth = 0.5) +
  scale_fill_manual(values = annotation_colors$Sample_Type, name = "Sample Type") +
  labs(
    title = "Median Absolute Deviation per Sample (Beta Values)",
    subtitle = "Genome-wide methylation variability across all CpGs (error bars = SE)",
    x = "Sample",
    y = "MAD (Beta values)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

############################################
## Plot 2: MAD per sample - M VALUES
############################################
p_mad_m <- ggplot(mad_m_df, aes(x = Sample, y = MAD, fill = Sample_Type)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = MAD - SE, ymax = MAD + SE),
                width = 0.3, linewidth = 0.5) +
  scale_fill_manual(values = annotation_colors$Sample_Type, name = "Sample Type") +
  labs(
    title = "Median Absolute Deviation per Sample (M-values)",
    subtitle = "Genome-wide methylation variability across all CpGs (error bars = SE)",
    x = "Sample",
    y = "MAD (M-values)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

############################################
## Plot 3: Grouped by phenotype - BETA VALUES
############################################
p_beta_grouped <- ggplot(mad_beta_df, aes(x = Sample_Type, y = MAD, fill = Sample_Type)) +
  geom_point(size = 3, shape = 21, color = "black", stroke = 0.5, 
             position = position_jitter(width = 0.1, seed = 42)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = annotation_colors$Sample_Type, name = "Sample Type") +
  labs(
    title = "MAD by Tissue Phenotype (Beta Values)",
    subtitle = "Individual samples shown as points",
    x = "",
    y = "MAD (Beta values)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11, face = "bold"),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))

############################################
## Plot 4: Grouped by phenotype - M VALUES
############################################
p_m_grouped <- ggplot(mad_m_df, aes(x = Sample_Type, y = MAD, fill = Sample_Type)) +
  geom_point(size = 3, shape = 21, color = "black", stroke = 0.5, 
             position = position_jitter(width = 0.1, seed = 42)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = annotation_colors$Sample_Type, name = "Sample Type") +
  labs(
    title = "MAD by Tissue Phenotype (M-values)",
    subtitle = "Individual samples shown as points",
    x = "",
    y = "MAD (M-values)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11, face = "bold"),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))

############################################
## Print all plots
############################################
print(p_mad_beta)
print(p_mad_m)
print(p_beta_grouped)
print(p_m_grouped)

############################################
## Save plots
############################################
ggsave("MAD_per_sample_BETA.pdf", p_mad_beta, width = 10, height = 6)
ggsave("MAD_per_sample_MVAL.pdf", p_mad_m, width = 10, height = 6)
ggsave("MAD_by_phenotype_BETA.pdf", p_beta_grouped, width = 7, height = 6)
ggsave("MAD_by_phenotype_MVAL.pdf", p_m_grouped, width = 7, height = 6)

############################################
## Summary statistics
############################################
cat("\n=== MAD SUMMARY BY PHENOTYPE (Beta Values) ===\n")
print(mad_beta_df %>%
        group_by(Sample_Type) %>%
        summarise(
          Mean_MAD = mean(MAD),
          Median_MAD = median(MAD),
          SD_MAD = sd(MAD),
          Min_MAD = min(MAD),
          Max_MAD = max(MAD),
          N = n()
        ))

cat("\n=== MAD SUMMARY BY PHENOTYPE (M-values) ===\n")
print(mad_m_df %>%
        group_by(Sample_Type) %>%
        summarise(
          Mean_MAD = mean(MAD),
          Median_MAD = median(MAD),
          SD_MAD = sd(MAD),
          Min_MAD = min(MAD),
          Max_MAD = max(MAD),
          N = n()
        ))

# Save data
write.csv(mad_beta_df, "MAD_per_sample_BETA.csv", row.names = FALSE)
write.csv(mad_m_df, "MAD_per_sample_MVAL.csv", row.names = FALSE)



