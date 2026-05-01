# Organoid Signature Scorer
# Calculate stomach and colon signature scores from bulk RNA-seq TPM data

# Set working directory
setwd('setworkingdirectory')

# Load TPM data
tpm_data <- read.delim('salmon.merged.transcript_tpm.tsv', stringsAsFactors = FALSE)

# Define gene signatures
stomach_genes <- c(
  "PGA4", "GIF", "GKN1", "PGA3", "PGA5", "ATP4A", "GAST", "ATP4B", "LIPF", "CLDN18",
  "GAGE12G", "MUC5AC", "GAGE2E", "DAZ4", "GKN2", "TFF1", "TFF2", "DAZ2", "KCNE2", "NKX6-2",
  "MUC6", "C6orf58", "PGC", "TAAR1", "CCKBR", "DPCR1", "FUT9", "CHIA", "CYP2C18", "CCKAR",
  "VSIG1", "COL2A1", "CTSE", "NTS", "A4GNT", "ANXA10", "AQP5", "CAPN8", "FAM159B", "GHRL",
  "PSCA", "ALDH3A1", "BARX1", "EPS8L3", "ERN2", "ESRRB", "PAX6", "PSAPL1", "C2orf70", "HNF4A",
  "IGFALS", "IHH", "SDR16C5", "SOX21", "TPH1", "TRIM31", "ACER2", "AKR1C2", "AKR7A3", "FER1L6",
  "FOXA3", "ISL1", "LRRC31", "MYEOV", "NPAS1", "RFLNA", "SLC26A9", "SLC9A4", "SULT1B1",
  "ABC13-47488600E17.1", "ACSM6", "ADH1C", "AGR2", "AKR1B10", "AMTN", "ANKRD22", "ANO7",
  "ARL14", "B3GNT6", "B4GALNT3", "BCAS1", "BPIFB1", "CA2", "CA9", "CAPN13", "CAPN9", "CFC1",
  "CLDN23", "CLIC6", "CRYBA2", "CXCL17", "CYP2S1", "DRD5", "EPS8L1", "FA2H", "FAM177B",
  "FAM83E", "FEV", "FOXA2", "FOXQ1", "FRMD1", "GALNT5", "GALNT6", "GIPR", "GPR25", "GUCA1C",
  "HAP1", "HTR1B", "INSM1", "KCNK16", "KPNA7", "LA16c-312E8.5", "LGALS9B", "LGALS9C", "LIME1",
  "LINC00675", "LIPH", "MBOAT4", "MFSD4A", "MIA", "MUC1", "MYRF", "NEUROD1", "NKX2-2", "NKX6-3",
  "NMUR2", "NPW", "NQO1", "OASL", "ONECUT3", "OVOL2", "PIK3C2G", "PLA2G10", "PRSS22", "PTGDR2",
  "RAB27B", "RASSF6", "REP15", "RFX6", "RNF223", "RP11-599B13.6", "S100P", "SCGN", "SHH",
  "SLC5A5", "SLC9A2", "SMIM6", "SOSTDC1", "SPTSSB", "SST", "SSTR1", "SULT1C2", "SYTL2",
  "SYTL5", "TESC", "TMED6", "TMEM211", "TMEM238", "TRIM50", "TRIM74", "TRNP1", "UNC5CL",
  "UPK1B", "VILL", "VSIG2", "ZSCAN4"
)

colon_genes <- c(
  "INSL5", "AQP8", "MEP1A", "PRAC1", "ISX", "CA1", "TBX10", "BTNL3", "MUC12", "TMIGD1",
  "KRTAP13-2", "MS4A12", "NOX1", "CD177", "GUCA2A", "PYY", "CHST5", "GUCY2C", "MYO1A", "NAT2",
  "NXPE1", "SDHD", "LYPD8", "MOGAT2", "SLC39A5", "AC011513.3", "C10orf99", "PIGY", "CHP2",
  "GLRA2", "HHLA2", "NR1I2", "ATOH1", "CDHR5", "CDX1", "LRRC19", "NXPE4", "REG4", "TMEM236",
  "CDH17", "CLCA1", "CLRN3", "DHRS11", "GAL3ST2", "GPA33", "HOXD13", "CEACAM7", "EPS8L3",
  "ERN2", "FABP1", "FSIP1", "GALNT8", "LRRC26", "PHGR1", "PPP1R14D", "SH2D7", "VIL1", "AMN",
  "BTNL8", "CASP5", "CDX2", "CEACAM6", "CLCA4", "EFNA2", "HAVCR1", "HNF4A", "IHH", "LINC01207",
  "MISP", "MOGAT3", "MUC13", "NEU4", "SLC17A4", "SPINK4", "TPH1", "TRIM31", "URAD", "ZG16",
  "B3GALT5", "FOXA3", "IL22RA1", "LRRC31", "NOS2", "SATB2", "SULT1B1", "TRPM5", "VIP",
  "AC009133.22", "ADAMDEC1", "AIFM3", "ATP10B", "B3GALT1", "B3GNT6", "B4GALNT2", "BCL2L15",
  "BEST2", "BEST4", "C15orf48", "C2orf72", "CA2", "CA4", "CA7", "CDC42EP5", "CEACAM1", "CEACAM5",
  "CES3", "CH17-360D5.1", "CKMT1B", "CLDN23", "CLDN3", "CTC-273B12.7", "DHRS9", "ENTPD8",
  "FAM3D", "FCGBP", "FFAR4", "FOXD2", "FRMD1", "FUT3", "FXYD3", "GCNT3", "GPR15", "GUCA2B",
  "HEPACAM2", "HOXD12", "HSD11B2", "ITLN1", "KLK15", "KRT20", "LEFTY1", "LGALS4", "LINC00675",
  "MAB21L2", "MUC4", "NOXO1", "NPY4R", "NRARP", "NXPE2", "OTOP2", "PADI2", "PIGR", "PIGZ",
  "PKIB", "PLA2G10", "PLA2G2A", "REP15", "RETNLB", "RP11-599B13.6", "RXFP4", "SLC22A18AS",
  "SLC26A2", "SLC26A3", "SLC9A2", "SLC9A3", "ST6GALNAC1", "TFF3", "TMEM171", "TPSG1", "TRABD2A",
  "TRIM15", "TRPM6", "TSPAN1", "TSPAN8", "UGT1A10", "UGT1A8", "UGT2B17"
)

# Function to calculate signature scores
calculate_signature_scores <- function(tpm_data, gene_col = 1) {
  # Convert gene column to uppercase for matching
  gene_names <- toupper(tpm_data[[gene_col]])
  
  # Get sample columns (all columns except gene column)
  sample_cols <- setdiff(names(tpm_data), names(tpm_data)[gene_col])
  
  # Initialize results data frame
  results <- data.frame(
    sample = sample_cols,
    stomach_score = numeric(length(sample_cols)),
    colon_score = numeric(length(sample_cols)),
    aggregate_score = numeric(length(sample_cols)),
    stomach_genes_found = integer(length(sample_cols)),
    colon_genes_found = integer(length(sample_cols)),
    stringsAsFactors = FALSE
  )
  
  # Calculate scores for each sample
  for (i in seq_along(sample_cols)) {
    sample <- sample_cols[i]
    
    # Extract TPM values for this sample
    tpm_values <- tpm_data[[sample]]
    
    # Get stomach gene values
    stomach_idx <- which(gene_names %in% toupper(stomach_genes))
    stomach_tpm <- tpm_values[stomach_idx]
    stomach_tpm <- stomach_tpm[stomach_tpm > 0 & !is.na(stomach_tpm)]
    stomach_log2 <- log2(stomach_tpm)
    
    # Get colon gene values
    colon_idx <- which(gene_names %in% toupper(colon_genes))
    colon_tpm <- tpm_values[colon_idx]
    colon_tpm <- colon_tpm[colon_tpm > 0 & !is.na(colon_tpm)]
    colon_log2 <- log2(colon_tpm)
    
    # Calculate scores
    results$stomach_score[i] <- ifelse(length(stomach_log2) > 0, 
                                       mean(stomach_log2), 
                                       NA)
    
    results$colon_score[i] <- ifelse(length(colon_log2) > 0, 
                                     mean(colon_log2), 
                                     NA)
    
    # Calculate aggregate score: average of inverted stomach and colon values
    inverted_stomach <- -1 * stomach_log2
    all_values <- c(inverted_stomach, colon_log2)
    results$aggregate_score[i] <- ifelse(length(all_values) > 0,
                                         mean(all_values),
                                         NA)
    
    results$stomach_genes_found[i] <- length(stomach_log2)
    results$colon_genes_found[i] <- length(colon_log2)
  }
  
  return(results)
}

# Calculate scores
scores <- calculate_signature_scores(tpm_data_clean, gene_col = 1)

# View results
print(scores)

# Save results
write.csv(scores, "organoid_signature_scores.csv", row.names = FALSE)

# Print summary
cat("\n===========================================\n")
cat("Signature Scoring Complete!\n")
cat("===========================================\n")
cat("Signature Gene Counts:\n")
cat("  Stomach genes:", length(stomach_genes), "\n")
cat("  Colon genes:", length(colon_genes), "\n")
cat("\nScoring Methodology:\n")
cat("  - Stomach Score: Average of log2(TPM) for stomach signature genes\n")
cat("  - Colon Score: Average of log2(TPM) for colon signature genes\n")
cat("  - Aggregate Score: Average of inverted stomach log2(TPM) and colon log2(TPM)\n")
cat("    * Higher aggregate = more colon-like\n")
cat("    * Lower aggregate = more stomach-like\n")
cat("\nResults saved to: organoid_signature_scores.csv\n")
cat("===========================================\n")


#PLOTTING
library(ggplot2)
library(reshape2)

# Reshape for ggplot
scores_long <- melt(scores[, c("sample", "stomach_score", "colon_score")], id.vars = "sample")

ggplot(scores_long, aes(x = sample, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal(base_size = 14) +
  labs(x = "Sample", y = "Mean log2(TPM)", fill = "Signature",
       title = "Stomach vs Colon Signature Scores per Sample") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(scores, aes(x = stomach_score, y = colon_score, label = sample)) +
  geom_point(size = 3, color = "steelblue") +
  geom_text(vjust = -0.7, size = 3.5) +
  theme_minimal(base_size = 14) +
  labs(x = "Stomach Signature (mean log2 TPM)",
       y = "Colon Signature (mean log2 TPM)",
       title = "Organoid Identity Space: Stomach vs Colon") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50")

ggplot(scores, aes(x = reorder(sample, aggregate_score), y = aggregate_score, fill = aggregate_score)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "darkred", high = "darkgreen") +
  theme_minimal(base_size = 14) +
  labs(x = "Sample", y = "Aggregate Score",
       title = "Aggregate Signature Score (Colon vs Stomach)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cor.test(scores$stomach_score, scores$colon_score)

ggplot(scores, aes(x = stomach_score, y = colon_score)) +
  geom_point(size = 3, color = "purple") +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  theme_minimal(base_size = 14) +
  labs(title = "Correlation Between Stomach and Colon Signature Scores",
       x = "Stomach Score", y = "Colon Score")

scores$group <- sub("[0-9]+", "", scores$sample)  # extract prefix before digits

ggplot(scores, aes(x = group, y = aggregate_score, fill = group)) +
  geom_boxplot() +
  theme_minimal(base_size = 14) +
  labs(title = "Aggregate Signature by Organoid Type",
       x = "Organoid Type", y = "Aggregate Score") +
  scale_fill_brewer(palette = "Set2")







