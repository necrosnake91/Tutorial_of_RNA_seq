##Script designed to perform DE analysis

#--------------------------------------Packages--------------------------------
##Load the required packages
packages <- c("tidyverse", "tximeta", "here", "DESeq2", "apeglm", "PCAtools", "Glimma", "RColorBrewer", "pheatmap", "clusterProfiler", "biomaRt", "enrichplot", "fgsea")
lapply(packages, library, character.only = T)
source("functions.R")
#--------------------------------------Data importation--------------------------------
##Read the metadata file
coldata <- read.table("../results/salmon_quants/metadata.txt", he = T, sep = "\t", stringsAsFactors = T)
##Add the path to quant files
coldata <- mutate(coldata, files = file.path(here("results/salmon_quants"), paste0(coldata$key, "_quant"), "quant.sf")) 
##Check if path to quant files is ok
coldata <- mutate(coldata, exist = file.exists(coldata$files))
##Create a new column consisting of key values  plus _quant
coldata <- mutate(coldata, names = paste0(coldata$key, "_quant"))
##Set rownames of coldata by using the unique id (treatment + replicate)
rownames(coldata) <- coldata$unique_id
##Read the files using tximeta
se <-  tximeta(coldata)
##Get the counts by sumarizing at gene level
gse <- summarizeToGene(se)
##For the shorthest path unlock the next line and run it
#load("../data/SummarizedExperimentObj.rdata")
#-------------------------------------Data exploration---------------------------------
##Create the DESeq object using the  information stored in the tximeta summarizedExperiment
dds <- DESeqDataSet(se = gse, 
                    design = ~ Treatment)
##Relevel the control group
levels(dds$Treatment)
dds$Treatment <- relevel(dds$Treatment, "siRNA_control")
##Remove genes with low abundance (low counts) using the following threshold
keep <- rowSums(counts(dds) >= 1) >= 3 ##Keep genes showing at least 1 raw count in at least 3 samples
##Visualize how many genes passed  the filtering criteria
table(keep)
##Cut the genes from the dds object
dds <-  dds[keep, ]
##Perform PCA analysis using the PCA tools package
rld <- rlog(dds, blind = F) ##Normalize the raw count values by using the regularized-logarithm transformation
pca <- pca(mat = assay(rld), metadata = colData(rld), scale = T) ##Create a PCA object
biplot(pca, lab = rownames(colData(rld)), colby = "Treatment") ##Get the biplot of the two main components
##Visualize the data using  an interactive MDS plot
glimmaMDS(dds)
##If necessary, batch effects must be reduced by aggregating a new column in coldata, called as batch. You can use the run ID or date of run
##for each sample. Then, when creating the dds object use: DESeqDataSet(se = gse, design = ~ batch + Treatment)
#-------------------------------------Differential expression analysis---------------------------------
##Run the DESeq function to normalize counts and adjust the data to the negative bi-nomial model
dds <- DESeq(dds)
##Get the results from the differential expression analysis
res <- results(dds, lfcThreshold = 1, alpha = 0.05) ##As default, we are testing for genes showing |lfc| > 0 and padj < 0.1
summary(res)
res <- as.data.frame(res)
volcanoplotR(res, logfc = 1, p.adj = 0.05)
##Compare the results from siRNA_NRF2 vs siRNA_control treatment
res_sirna <- results(dds, lfcThreshold = 0, alpha = 0.01)
summary(res_sirna)
##Visualize the data using a MA-plot
tiff("../results/Maplot.tiff", res = 300, height = 8, width = 10, units = "in")
plotMA(res_sirna)
dev.off()
##Shrunk the logFC for genes with low counts
resultsNames(dds)
res_shrink <- lfcShrink(dds, res = res_sirna, coef = "Treatment_siRNA_NRF2_vs_siRNA_control", type = "apeglm")
plotMA(res_shrink)
res_shrink <- as.data.frame(res_shrink)
res_sirna <- as.data.frame(res_sirna)
volcanoplotR(res_shrink, logfc = 0, p.adj = 0.01)
##Get the list of differentialy expressed genes
deg <- dplyr::filter(res_shrink, log2FoldChange > 0 & padj < 0.01 |
                log2FoldChange < 0 & padj < 0.01)
##Get the normalized counts matrix
norm_counts <- counts(dds, normalized = T)
##Plot the heatmap
##Obtain the annotation for the columns of the heatmap
annotation_col <- data.frame(coldata[1:6, c(2, 4)])
##Select a nice palette
RdBlu <- rev(brewer.pal(n= 10, name = "RdBu"))
pheatmap(norm_counts[rownames(deg), 1:6], scale = "row", 
         border_color = NA, show_rownames = F, clustering_distance_rows = "euclidean",  
         clustering_distance_cols = "euclidean", clustering_method = "single", show_colnames = F, 
         annotation_col = annotation_col, 
         color = RdBlu)
#-------------------------------------Annotation analysis---------------------------------
#-------------------------------------ORA---------------------------------
##Convert the ensembl gene ids into entrezgene and hgnc symbol
##Retrieve information from ensembl database
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
##Get the list of up- and down regulated genes differentially expressed
up_DEG  <- dplyr::filter(deg, log2FoldChange > 0) %>% rownames_to_column(var = "ensembl_gene_id")
down_DEG <- dplyr::filter(deg, log2FoldChange < 0) %>% rownames_to_column(var = "ensembl_gene_id")
##Covert the ids
up_DEG <- id_converter(mart = ensembl, ##Object with ensembl db information
            input = up_DEG, 
            attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "gene_biotype"), ##Name of the attributes to get in the output
            filter = "ensembl_gene_id") ##Name of the attribute (column) of the input to match with db
down_DEG <- id_converter(mart = ensembl,
                         input = down_DEG, 
                         attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "gene_biotype"), 
                         filter = "ensembl_gene_id")
##Read the gmt files for GOBP database
go <- read.gmt("../data/c5.go.bp.v7.5.1.entrez.gmt")
##Perform ORA
ego_up <- enricher(gene = up_DEG$entrezgene_id, ##Select he column of entrez gene ids from the converted list
                   TERM2GENE = go, ##This is the background information. In this case use the GO db
                   pAdjustMethod = "BH", ##For multiple testing, adjust p values by using the  Benjamini-Hoechberg method
                   pvalueCutoff = 0.05) ##Cutoff to keep pathways or gene sets showing  padj < 0.05
ego_down <- enricher(gene = down_DEG$entrezgene_id, 
                     TERM2GENE = go, 
                     pAdjustMethod = "BH", 
                     pvalueCutoff = 0.05)
##Visualize the results
dotplot(ego_up, showCategory = 15)
barplot(ego_up, showCategory = 15)
ego_up <- pairwise_termsim(ego_up)
treeplot(ego_up)
#-------------------------------------GSEA---------------------------------
##For GSEA analysis, use the results obtained after the differential expression analysis
res_shrink <- rownames_to_column(res_shrink, var = "ensembl_gene_id")
##Convert ensembl gene ids from  results data frame
res_shrink <- id_converter(mart = ensembl,  
                           input = res_shrink, 
                           attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "gene_biotype"), 
                           filter = "ensembl_gene_id")
##Calculate the metric to rank the genes. Metric = -log10(padj)*log2FoldChange
res_shrink <- mutate(res_shrink, stat = -log10(res_shrink$padj)*res_shrink$log2FoldChange) %>%
  arrange(desc(stat)) %>% ##Sort the data frame respect to the stat metric
  drop_na() ##Eliminate NA values
##For GSEA input create a named vector using the stat value
gsea_list <- res_shrink$stat
names(gsea_list) <- res_shrink$entrezgene_id ##Name the elements using the entrezgene id
##Load the database again
go <- gmtPathways("../data/c5.go.bp.v7.5.1.entrez.gmt")
##Perform GSEA
GSEA_res <- fgseaMultilevel(go, 
                            gsea_list, 
                            minSize = 15, ##Exclude gene sets with less than 15 genes
                            maxSize = 500) ##Exclude gene sets with more than 500 genes
##Tidy the results and add a column depending on the value of padj
GSEA_res <- as.tibble(GSEA_res) %>%
  arrange(desc(NES)) %>%
  mutate(Significance = ifelse(padj < 0.05, "Significant", "NS"))
##Plot the results
barplot_GSEA(GSEA_res, Head = 20, Tail = 20)
