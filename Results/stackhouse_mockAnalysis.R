## Yeh Lab Mock Analysis
## Author: Christian Stackhouse
## Version: 0.0.1
## Date: 10.1.2020
## Contact: ctstackh@uab.edu

#### Instructions ####
# - TCGA_PAAD_plus.RData contains 181 PDAC patient samples of RNA-seq data.
# - ex: unlogged TPM values; sampInfo: sample information; featInfo: gene information.
# - Please only consider 150 samples with non-NA sample information in the Grade column for all the analysis.
# - Please submit the code for analysis and result figures.
# - Please push your results as a separate folder with your name to this same repo.

# Set WD and Import Data
setwd("C:/Users/Christian/Desktop/YehLab/MockAnalysis")
load("TCGA_PAAD_plus.RData")
aurka <- read.delim(file = "AURKA_sig_genes.txt", sep = "\t", header = F)
basal <- read.delim(file = "Moffitt-basal25.txt", sep = "\t", header = F)
classical <- read.delim(file = "Moffitt-classical25.txt", sep = "\t", header = F)

exp <- TCGA_PAAD_plus$ex
rownames(exp) <- TCGA_PAAD_plus$featInfo$SYMBOL
na_list <- is.na(TCGA_PAAD_plus$sampInfo$Grade)
subtype <- data.frame(TCGA_PAAD_plus$sampInfo$mRNA.Moffitt.clusters.All.150.Samples.1basal.2classical)
subtype <- na.omit(subtype)
colnames(subtype) <- "Subtype"
exp <- exp[, !na_list]
exp <- log2(exp+1)

#### Q1 ####
# Perform a consensus clustering for 150 samples using genes in Moffitt-basal25.txt (orange) and Moffitt-classical25.txt (blue). 
# Use K=2 for both column and row. Plot a heatmap to show the result. The resultant two clusters of samples are basal-like (orange)
# and classical (blue) subtypes respectively. 
library(BiocManager)
BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)
library(ALL)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(grid)
library(tidyr)

exp_b <- exp[rownames(exp)%in%basal$V1,]
exp_c <- exp[rownames(exp)%in%classical$V1,]
exp_new <- rbind(exp_b,exp_c)

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)

subtype[] <- lapply(subtype, as.character)
subtype <- subtype %>% mutate(Subtype = case_when(
    subtype == 1 ~ "Basal",
    subtype == 2 ~ "Classical"
))

row.names(subtype) <- colnames(exp_new)

pdf("Q1 HC of Samples.pdf")

pheatmap(exp_new, cluster_cols = T, cluster_rows = T, scale = "none", main = "HC of Samples", 
         color = myColor, show_rownames = T, annotation_col = subtype, annotation_colors = list(Subtype=c(Basal="orange",Classical="blue")),
         cutree_rows = 2, cutree_cols = 2, show_colnames = F)
dev.off()


#### Q3 ####
# Identify all the differentially expressed (DE) genes between basal-like and classical samples using DESeq2 
# Plot MA plot and volcano plot to show the results. Cut-offs to use: p=0.05, log2(fold-change)=1. 
# Read counts are in file TCGA_PAAD_plus.cnt.txt.gz. Patterns of sample ID need to be used to match to samples in .RData.
library(DESeq2)

exp <- TCGA_PAAD_plus$ex
rownames(exp) <- TCGA_PAAD_plus$featInfo$SYMBOL
na_list <- is.na(TCGA_PAAD_plus$sampInfo$Grade)
colnames(subtype) <- "Subtype"
exp <- exp[, !na_list]
exp <- round(exp,0)

dds <- DESeqDataSetFromMatrix(countData = exp,                                  
                             colData = coldata,
                             design= ~ Subtype)
dds <- DESeq(dds)

res <- results(dds, alpha = 0.05)

pdf("MA-plot.pdf")

plotMA(res,ylim=c(-1,1))

dev.off()

res1 <- res
pVal <- -log10(res1[,6])
res1 <- cbind(res1,pVal)
res1 <- res1[,c(-1,-3,-4,-5,-6)]
logFC <- res1[,1]

thresh <- -log10(.05)
lfc <- 1

pdf("Volcano_Plot.pdf")

volcanoData <- cbind(logFC,pVal)
colnames(volcanoData) <- c("LogFC", "-LogPValue")
plot(volcanoData, pch=1, cex = 0.5) 
abline(h=thresh, col="red", lty=2)
abline(v=c(lfc,-lfc), lty=2)

dev.off()

#### Q2 ####
# Test if each of the genes in AURKA_sig_genes.txt is differentially expressed between basal-like and classical subtypes. 
# Use boxplot to illustrate the results with p-values labeled on the figures.

write.table(res,"DESeq_results.txt", sep = "\t", quote = F)
res <- read.delim("DESeq_results.txt", sep = "\t", header = T)

for (i in 1:length(aurka$V1))
{
  
goi <- aurka[i,1]
pvalue <- res[rownames(res)%in%goi,6]

data <- plotCounts(dds, gene=goi, intgroup=c("Subtype"), returnData=TRUE)

pdf(paste0(goi," ","Boxplots.pdf"))
print(ggplot(data, aes(x=Subtype, y=count)) +
  scale_y_log10() + 
  geom_point(position=position_jitter(width=.1,height=0), size=3) +
  geom_boxplot() +
  ggtitle(paste0(goi," ","P_Val=",pvalue)))
dev.off()

i == i+1

}


#### Q4 ####
# Is subtype assocaited with the clinical variables (sex, age, gender and race)? Use a table with p-values to show the results.

info <- data.frame(ID = TCGA_PAAD_plus$sampInfo$Tumor.Sample.ID,
                   gender = TCGA_PAAD_plus$sampInfo$Gender,
                   age = TCGA_PAAD_plus$sampInfo$Age.at.initial.pathologic.diagnosis,
                   race = TCGA_PAAD_plus$sampInfo$Race,
                   subtype = TCGA_PAAD_plus$sampInfo$SSC.subtype)
info <- info[!na_list,]

gender <- kruskal.test(info$subtype,info$gender)
prob_gender <- gender$p.value  ## 0.642

age <- kruskal.test(info$subtype,info$age)
prob_age <- age$p.value  ## 0.621

race <- kruskal.test(info$subtype,info$race)
prob_race <- race$p.value  ## 0.334

grid.table()

# Save and finish

save.image("stackhouse_mockAnalysis.RData")
