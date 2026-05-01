###Illumina-EPICv2 Full Analysis

workdir= "setworkingdirectory"
setwd(workdir)

install.packages("ggsignif")
install.packages("gplots")


BiocManager::install("DNAmArray")
BiocManager::install("sva")
BiocManager::install("DNAcopy")
BiocManager::install("impute")
BiocManager::install("wateRmelon")
BiocManager::install("ggfortify")
BiocManager::install("irlba")
BiocManager::install("devtools")
BiocManager::install("fastcluster")
BiocManager::install("ggpubr")
BiocManager::install("qqman")
BiocManager::install("methylGSA")
BiocManager::install("GCally")




library(devtools)
httr::set_config(httr::config( ssl_verifypeer = 0L))
install_github("molepi/DNAmArray") ##for master
install_github("molepi/DNAmArray", ref="R-3.3.0") ##for other branches

library(ChAMP)
library(knitr)
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


install.packages("setdirectory/ChAMP", repos = NULL, type = "source")


### ChAMP
#Check the CSV, the sentrix ID often gets converted into scientifc notation, and it needs to be in text, putting a single coma fixes.
#this command isn't working at current: myLoad <- champ.load("setdirectory/Methylation/207695700005", arraytype = "EPICv2")
myImport <- champ.import("setdirectory/Methylation/207695700005", arraytype = "EPICv2")
myLoad <-  champ.filter(beta = myImport$beta,
                          pd=myImport$pd,
                          detP=myImport$detP,
                          beadcount=myImport$beadcount,
                          ProbeCutoff=0.1,
                          arraytype = "EPICv2")

CpG.GUI(arraytype="EPICv2")

champ.QC(beta=myLoad$beta, pheno=myLoad$pd$Sample_Type)
QC.GUI(beta=myLoad$beta, pheno=myLoad$pd$Sample_Type, arraytype = "EPICv2")

##Normalization
myNorm <- champ.norm(myLoad$beta, method = "BMIQ", arraytype = "EPICv2")
myNorm2 <- champ.norm(myLoad$beta, method = "PBC", arraytype = "EPICv2")
# myNorm[is.na(myNorm)] <- median(myNorm, na.rm=TRUE)
# myNorm2[is.na(myNorm2)] <- median(myNorm2, na.rm=TRUE)

write.table(myNorm, "myNorm_samples.txt", row.names = TRUE, sep = "\t")

champ.QC(beta=myNorm, pheno = myLoad$pd$Sample_Type)
champ.QC(beta=myNorm2, pheno = myLoad$pd$Sample_Type)

QC.GUI(beta=myNorm, pheno = myLoad$pd$Sample_Type, arraytype = "EPICv2")
QC.GUI(beta=myNorm2, pheno = myLoad$pd$Sample_Type, arraytype = "EPICv2")

# Find the 1000 most variable CpGs based on standard deviation
most_variable_cpgs <- names(sort(apply(myNorm, 1, sd), decreasing = TRUE))[1:1000]

# Extract the table for the 1000 most variable CpGs
table_1000_cpgs <- myNorm[most_variable_cpgs, ]

# Now you can write this table to a file (e.g., a tab-separated values file)
write.table(table_1000_cpgs, "table_1000_variable_CpGs.txt", row.names = TRUE, sep = "\t")

# If you want to view library(ggplot2)
#the table in R, you can print the first few rows using the head() function
head(table_1000_cpgs)

##Making graph of average B value per condition
# Assuming 'myNorm' is your normalized beta matrix
beta_df <- as.data.frame(myNorm)
# Convert rownames to a column to keep track of CpG sites
beta_df$CpG_ID <- rownames(beta_df)
# Extract sample types from column names
sample_labels <- colnames(beta_df)[-ncol(beta_df)]  # Exclude the CpG_ID column
# Extracting broad sample categories like "Gastric", "IM", "Duodenal"
sample_types <- gsub("\\d+", "", sample_labels)  # Remove patient numbers
# Create a dataframe mapping samples to their sample types
sample_types_df <- data.frame(Sample = sample_labels, Sample_Type = sample_types)
# Reshape beta values to long format and join with sample types
long_beta_df <- beta_df %>%
  pivot_longer(cols = -CpG_ID, names_to = "Sample", values_to = "Beta_Value") %>%
  left_join(sample_types_df, by = "Sample")
# Calculate average beta values for each CpG site and sample type
avg_beta_df <- long_beta_df %>%
  group_by(Sample_Type) %>%
  summarise(Average_Beta = mean(Beta_Value, na.rm = TRUE), .groups = 'drop')
# Generate the plot using ggplot2
ggplot(avg_beta_df, aes(x = Sample_Type, y = Average_Beta, fill = Sample_Type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Sample Type", y = "Average Beta Value", 
       title = "Average Beta Values for Gastric, IM, and Duodenal Samples") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#Boxplot rather than barchart
sample_types_df <- data.frame(Sample = sample_labels, Sample_Type = sample_types)
# Reshape beta values to long format and join with sample types
long_beta_df <- beta_df %>%
  pivot_longer(cols = -CpG_ID, names_to = "Sample", values_to = "Beta_Value") %>%
  left_join(sample_types_df, by = "Sample")
# Generate the boxplot using ggplot2
p <- ggplot(long_beta_df, aes(x = Sample_Type, y = Beta_Value, fill = Sample_Type)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Sample Type", y = "Beta Value", 
       title = "Distribution of Beta Values for Gastric, IM, and Duodenal Samples") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# Add p-values to the plot using pairwise comparisons (Wilcoxon test by default)
p + stat_compare_means(method = "wilcox.test", 
                       comparisons = list(c("Gas", "IM"), c("Gas", "Duo"), c("IM", "Duo")),
                       label = "p.signif",
                       tip.length = 0.02)
                     


###########################
#DMP Generation
###########################

tmp_pheno <- myLoad$pd$Sample_Type
tmp_pheno <- gsub(" ","_",tmp_pheno)
myDMP1 <- champ.DMP(beta = myNorm,pheno=tmp_pheno, arraytype = "EPICv2", compare.group = c("Gas", "IM"))
DMP.GUI(DMP=myDMP1[[1]], beta=myNorm, pheno=myLoad$pd$Sample_Type)

myDMP2 <- champ.DMP(beta = myNorm,pheno=tmp_pheno, arraytype = "EPICv2", compare.group = c("IM", "Duo"))
DMP.GUI(DMP=myDMP2[[1]], beta=myNorm, pheno=myLoad$pd$Sample_Type)

myDMP3 <- champ.DMP(beta = myNorm,pheno=tmp_pheno, arraytype = "EPICv2", compare.group = c("Gas", "Duo"))
DMP.GUI(DMP=myDMP3[[1]], beta=myNorm, pheno=myLoad$pd$Sample_Type)

write.table(myDMP1, "myDMPs1_Gas_vs_IM.txt", row.names = TRUE, sep = "\t")
write.table(myDMP2, "myDMPs2_Duo_vs_IM.txt", row.names = TRUE, sep = "\t")
write.table(myDMP3, "myDMPs3_Gas_vs_Duo.txt", row.names = TRUE, sep = "\t")

##Trying to plot something with my DMP.
fdr_threshold <- 0.05
significant_dmp <- myDMP1$Gas_to_IM[myDMP1$Gas_to_IM$adj.P.Val < fdr_threshold, ]
# Define the number of top significant results to display
top_n <- 50
# Select top N by magnitude of methylation difference
top_dmp_by_diff <- significant_dmp[order(abs(significant_dmp$deltaBeta), decreasing = TRUE), ][1:top_n, ]
# Create a plot using ggplot2
ggplot(top_dmp_by_diff, aes(x = reorder(gene, deltaBeta), y = deltaBeta)) +
  geom_bar(stat = "identity", aes(fill = ifelse(deltaBeta > 0, "Increased", "Decreased")), 
           position = "dodge") +
  theme_minimal() +
  labs(x = "Gene", y = "Methylation Difference (Delta Beta)", 
       title = "Top Genes with Highest Methylation Differences") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip()

##Plotting volcano plots based on DMPs
dmp_data1 <- myDMP1$Gas_to_IM
# Calculate -log10(p-value) for the y-axis
dmp_data1 <- dmp_data1 %>%
  mutate(neg_log10_pvalue = -log10(P.Value),
         Significance = ifelse(adj.P.Val < 0.05 & abs(logFC) > 0.5, "Significant", "Not Significant"))
# Create the volcano plot with more restrictive thresholds
volcano_plot <- ggplot(dmp_data1, aes(x = logFC, y = neg_log10_pvalue)) +
  geom_point(aes(color = Significance), alpha = 0.5, size = 1) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  labs(x = "Log2 Fold Change", y = "-Log10(p-value)", title = "Volcano Plot of DMPs (Gas vs IM)") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue")
# Display the plot
print(volcano_plot)

##############################
#DMR Generation
##############################

myDMR1 <- champ.DMR(beta=myNorm, pheno=tmp_pheno,method="Bumphunter", arraytype = "EPICv2", compare.group = c("Gas", "IM"))
myDMR2 <- champ.DMR(beta=myNorm, pheno=tmp_pheno,method="Bumphunter", arraytype = "EPICv2", compare.group = c("Duo", "IM"))
myDMR3 <- champ.DMR(beta=myNorm, pheno=tmp_pheno,method="Bumphunter", arraytype = "EPICv2", compare.group = c("Gas", "Duo"))



write.table(myDMR1, "myDMRs1_Gas_vs_IM.txt", row.names = TRUE, sep = "\t")
write.table(myDMR2, "myDMRs2_Duo_vs_IM.txt", row.names = TRUE, sep = "\t")
write.table(myDMR3, "myDMRs3_Gas_vs_Duo.txt", row.names = TRUE, sep = "\t")

DMR.GUI(DMR=myDMR1, beta=myNorm, runDMP= FALSE, pheno=myLoad$pd$Sample_Type, arraytype = "EPICv2")
DMR.GUI(DMR=myDMR2[[1]], beta=myNorm, pheno=myLoad$pd$Sample_Type, )
DMR.GUI(DMR=myDMR2[[1]], beta=myNorm, pheno=myLoad$pd$Sample_Type)

myDMRall <- champ.DMR(beta=myNorm, pheno=tmp_pheno,method="Bumphunter", arraytype = "EPICv2")
DMR.GUI(DMR=myDMRall, beta=myNorm, runDMP= TRUE, pheno=myLoad$pd$Sample_Type, arraytype = "EPICv2", compare = c("Gas","IM"))

#######################
#Minfi + DNAmarray can be used rather than EPICv2
#######################

targets <- read.metharray.sheet("setdirectory/Methylation/207695700005")
RGSet <- read.metharray.exp(targets = targets,force=TRUE,extended=TRUE)
# setting the right array type and annotation will do the work:
RGSet@annotation <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
RGSet
#Filtering probes
RGSet_filtered = DNAmArray::probeFiltering(RGSet)
betas <- getBeta(RGSet, type="Illumina")
#betadensityPlots
ggbg <- function() {
  points(0, 0, pch=16, cex=1e6, col="grey90")
  grid(col="white", lty=1)
}

par(mar=c(4,4,3,2), mgp=c(2.5,1,0), 
    cex.main=1.5, font.main="1", 
    fg="#6b6b6b", col.main="#4b4b4b")

densityPlot(RGSet, 
            main="Beta density plot", 
            xlab="Beta values", 
            panel.first=ggbg()) 
##PC Plot
pc <- prcomp_irlba(t(betas), n=6)
summary(pc)
# Create the principal components plot
p <- autoplot(pc, data = targets, colour = "Sample_Type", main = "Principal Components Plot")
# Add sample names as labels using ggrepel
p <- p + geom_text_repel(aes(label = Sample_Name), size = 3, show.legend = FALSE, max.overlaps = 20)
# Display the plot
print(p)
# Save the plot as a PDF file
ggsave("PC_plot.pdf", plot = p, width = 8, height = 6)

##QC
MSet <- preprocessRaw(RGSet)
qc <- getQC(MSet)
addQC(MSet,qc)
MSet=fixMethOutliers(MSet, K = -3, verbose = FALSE)

pdf("QC_Report_MSet.pdf")
par(mar=c(3,10,3,10),cex.axis=0.7)
plotQC(qc)
densityPlot(MSet, sampGroups = MSet$Sample_Type)
densityBeanPlot(MSet, sampNames = MSet$Sample_Name)
dev.off()

pdf("QC_Report_RGSet.pdf")
par(mar=c(3,10,3,10),cex.axis=0.7)
plotQC(qc)
densityPlot(MSet, sampGroups = RGSet$Sample_Group)
densityBeanPlot(MSet, sampNames = RGSet$Sample_Name)
dev.off()

#QC report
qcReport(RGSet, pdf= "QC_Report_v2_RGSet.pdf")
par(mar=c(3,10,3,10),cex.axis=0.7)
phenoData <- pData(RGSet)
manifest <- getManifest(RGSet)

#Sex Prediction
sex <- getSex(mapToGenome(RGSet))
pdf("Gender_sample_names.pdf")
ggplot(as.data.frame(sex), aes(x = xMed, y = yMed, fill = predictedSex)) +
  geom_point(shape = 21, size = 3) +
  geom_text(aes(label = targets$Sample_Name), vjust = 2) +
  theme_light() +
  scale_fill_manual(values = c("pink", "blue"), labels = unique(sex$predictedSex)) +
  xlab("X chr, median total intensity (log2)") +
  ylab("Y chr, median total intensity (log2)") +
  guides(fill = guide_legend(title = "Predicted Sex")) +
  theme(legend.position = "top")
dev.off()

#Map to genome
RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)

#Normalization- Functional Normalization
pc <- screeplot(RGSet)
RGset.funnorm <- preprocessFunnorm(RGSet,sex=sex$predictedSex,nPCs=6)
mybeta <- reduce(RGset.funnorm,RG,what="beta")
write.table(mybeta, "mybeta.txt", row.names = TRUE, sep = "\t")  
most_variable_cpgs_mybeta <- names(sort(apply(mybeta, 1, sd), decreasing = TRUE))[1:1000]
table_1000_cpgs_mybeta <- mybeta[most_variable_cpgs_mybeta, ]
# Now you can write this table to a file (e.g., a tab-separated values file)
write.table(table_1000_cpgs_mybeta, "table_1000_variable_CpGs_mybeta.txt", row.names = TRUE, sep = "\t")

#MDS plot for RAW
pdf("MDS_RAW.pdf", width = 10, height = 8)
mdsPlot(getM(MSet), 
        numPositions = 1000, 
        sampGroups = targets$Sample_Type,
        sampNames = targets$Sample_Name,
        legendPos = "bottomleft",
        main = "RAW")
dev.off()
#MDS plot for funorm
pdf("MDS_Funnorm.pdf", width = 10, height = 8)
mdsPlot(getM(RGset.funnorm ), 
        numPositions = 2000, 
        sampGroups = targets$Sample_Type,
        sampNames = targets$Sample_Name,
        legendPos = "bottomleft",
        main = "preprocessFunnorm")
dev.off()

# PCA plot and heatmap
beta <- getBeta(RGset.funnorm)
beta[is.na(beta)] <- median(beta, na.rm=TRUE)
pc <- prcomp_irlba(t(beta), n=6)
summary(pc)
plot(pc$x)
# Making the plot for the correlations of PCs
library(pheatmap)
# Plot the heat map
df <- apply(targets, 2, function(x) as.numeric(factor(x)))
# if there is na in the beta values; fill them with median value
df[is.na(df)] <- median(df, na.rm=TRUE)
keep <- apply(df, 2, sd) > 0    
cxy <- cor(pc$x, scale(df[, keep]))   
pheatmap(cxy[,c(2:3,4,5)], cluster_rows=FALSE)

#Plot dendogram
# Change column names to sequential numbers
colnames(beta) <- 1:ncol(beta)
library(fastcluster)
d <- dist(t(beta),method="euclidean")
fit <- hclust(d, method="average")
Interesting <- sample(colnames(beta),1)

colorLeafs <- function(x) {
  if (is.leaf(x) && attr(x, "label") %in% Interesting) {
    attr(x, "nodePar") <- list(lab.col="red", pch=NA)
  }
  return(x)
} 
dd <- dendrapply(as.dendrogram(fit), colorLeafs)
op <- par(mar=c(10,4,4,2))
plot(dd)



d <- dist(t(beta),method="euclidean")
fit <- hclust(d, method="average")
Interesting <- sample(colnames(beta),1)

colorLeafs <- function(x) {
  if (is.leaf(x)) {
    label <- attr(x, "label")
    if (label %in% colnames(beta)[1:23]) {
      attr(x, "nodePar") <- list(lab.col = "blue", pch = NA)
    } else if (label %in% colnames(beta)[24:41]) {
      attr(x, "nodePar") <- list(lab.col = "red", pch = NA)
    }
  }
  return(x)
}

dd <- dendrapply(as.dendrogram(fit), colorLeafs)
op <- par(mar = c(10, 4, 4, 2))
plot(dd)

