##Script designed to perform DE analysis

#--------------------------------------Packages--------------------------------
##Load the required packages
library(ggplot2)
library(tidyverse)
library(edgeR)
library(DESeq2)
library(PCAtools)
library(marray)
library(pheatmap)
source("functions.R")

#--------------------------------------Data importation--------------------------------
##In this case, count matrix is provided as a .txt file
counts <- read.table("../results/counts.txt", he = T)

#--------------------------------------Data exploration--------------------------------
##First, genes with low abundance must be excluded. Select genes with at least 3 cpm in at leats 2 samples
keep <- rowSums(cpm(counts, log = T) > 3) > 2
##Visualize how many genes passed the filtering criteria
table(keep)
##Cut the genes which passed filtering criteria in the original matrix
counts <- counts[keep, ]

##Create a factor object to indicate the name of the experimental conditions
groups <- factor(sub("..$", "", names(counts))) ##"..$" indicates substitute the two last characters with "" (nothing)
table(groups)

##Store the counts and the groups in an edgeR object
edgeRlist <- DGEList(counts = counts, 
                     group = groups, 
                     genes = rownames(counts))
##Calculate the normalization factors using the TMM method
edgeRlist <- calcNormFactors(edgeRlist, method = "TMM")
##Visualize the normalization factors
edgeRlist$samples
##Plot the results using absolute vs relative expression in each sample
pdf("../results/MD_plots.pdf", height = 7, width = 10)
par(mfrow = c(2, 6)) ##Generate a frame to store 6 plots in 2 rows and 3 columns
for (i in c(1:12)) {
  print(plotMD(cpm(edgeRlist, log = T), column = i))
  grid(col = "blue")
  abline(h = 0, col = "red", lty = 2, lwd = 2)
}
dev.off()
##Inspect replicates by performing a PCA analysis
pca <- pca(cpm(edgeRlist$counts, log = T))
##Plot the results
png("../results/PCA.png", height = 700, width = 800)
biplot(pca, lab = colnames(edgeRlist$counts), pointSize = 15, title = "PCA", labSize = 10)
dev.off()

#--------------------------------------Differential expression analysis--------------------------------
##Get the design matrix
design <- model.matrix(~0+edgeRlist$samples$group)
##Write the colnames of the design matrix by using the levels of the experimental groups
colnames(design) <- levels(edgeRlist$samples$group)
design
##Estimate data dispersion
edgeRlist <- estimateDisp(edgeRlist, design = design, robust = T)
##Visualize the dispersion levels
png("../results/data_dispersion.png", height = 700, width = 800)
plotBCV(edgeRlist)
dev.off()
##Construct the contrast matrix. In this case we are going to compare poison treated cells vs CT
contrast <- makeContrasts(
  "Poison" = "A_Verafinib - A_Control",
  levels = edgeRlist$design
)
contrast

##Adjust data to a negative bi-nomial generalized linear model and using the trended dispersion
fit <- glmQLFit(edgeRlist, design = design, robust = T,  dispersion = edgeRlist$trended.dispersion)
##Test the null hypothesis in which lfc of genes in the poison group are equal to zero respect to CT
qlf <- glmQLFTest(fit, contrast = contrast[, "Poison"])
##Visualize how many genes rejected the null hypothesis with a FDR threshold of 0.05
degPoison_vs_CT <- decideTestsDGE(qlf, p.value = 0.05, adjust.method = "BH", lfc = 0)
table(degPoison_vs_CT)

##Store the results in a df
DEGPoison_vs_CT <- DEGResults(qlf)
##Create a volcano plot by adding a new column according to |lfc| > 0 and FDR < 0.05 values
DEGPoison_vs_CT <- edgeResults(DEGPoison_vs_CT, logfc = 0, padj = 0.05)
##Get the volcano plot
png("../results/Volcano_plot.png", height = 600, width = 550)
volcano_edgeR(DEGPoison_vs_CT, lfc = 0, padj = 0.05)+
  xlim(c(-5, 5))
dev.off()

##Get all significant DE genes
significant_genes <- DEGPoison_vs_CT %>% filter(logFC > 0 & FDR < 0.05 | logFC < 0 & FDR < 0.05)
significant_ids <- significant_genes$genes
##Get the cpm values of significant DE genes
significant_cpm <- cpm(edgeRlist$counts, log = T) 
##Cut the significant DE genes 
significant_cpm <- significant_cpm[significant_ids, ]
##Create a nice color palette (You are able to change the colors anytime!)
PinkBlue <- maPalette(low = rgb(1, 0, 1), high = rgb(0, 0, 1), mid = rgb(0.6, 0.5, 1))
##Plot a heatmap
png("../results/Heatmap.png", height = 600, width = 800)
pheatmap(significant_cpm, 
         color = PinkBlue, 
         border_color = NA, 
         show_rownames = F, 
         scale = "row", 
         angle_col = 0)
dev.off()

##Save the results
write.table(DEGPoison_vs_CT, "../results/All_DEG_Verafinib.txt", sep = "\t", quote = F)
write.table(significant_genes, "../results/Significant_DEG_Verafinib.txt", sep = "\t", quote = F)
