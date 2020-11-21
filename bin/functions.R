DEGResults <- function(qlf) {
  ##This function returns a dataframe with all DEG
  ##qlf: Object obatined from the generalized linear model
  qlf <- topTags(qlf, n = Inf)
  qlf <- as.data.frame(qlf) 
  return(qlf)
}

edgeResults <- function(dge.obj, logfc, padj) {
  ##This function adds a new column (T or F) according to the FDR and LFC of each gene in edgeR list of DEG
  ##dge.obj: List with DEG
  volc <- dge.obj %>%
    mutate(condition = ifelse((dge.obj$logFC > logfc) & (dge.obj$FDR < padj), "Over-expressed",
                              ifelse((dge.obj$logFC < -logfc) & (dge.obj$FDR < padj), "Sub-expressed",
                                     ifelse((dge.obj$logFC > logfc) & (dge.obj$FDR > padj), "NS",
                                            ifelse((dge.obj$logFC < -logfc) & (dge.obj$FDR > padj), "NS",
                                                   ifelse((dge.obj$logFC < logfc) & (dge.obj$FDR > padj), "NS", "NS"))))))
  return(volc)
}

volcano_edgeR <- function(volc, lfc, padj) {
  ##The function returns the volcano plot of the dataset
  ##volcano_data: Data frame containing the data to create the volcano plot
  volcano_plot <- ggplot(volc)+
    geom_point(aes(x = logFC, y = -log10(FDR), color = condition))+
    scale_color_manual(name = "Condition",
                       labels = paste(c("NS", "Over-expressed", "Sub-expressed"), c(sum(volc$condition == "NS"), sum(volc$condition == "Over-expressed"), sum(volc$condition == "Sub-expressed"))),
                       values = c("#6e6d6e","#d84b47","#66c343"))+
    geom_vline(aes(xintercept = lfc), linetype = "dashed")+
    geom_vline(aes(xintercept = -lfc), linetype = "dashed")+
    geom_hline(aes(yintercept = -log10(padj)), linetype = "dashed")+
    theme_set(theme_bw())+
    theme(plot.title = element_text(face = "bold", size = 18), 
          axis.title = element_text(size = 18),
          legend.title = element_text(face = "bold", size = 15),
          legend.text = element_text(size = 15), 
          legend.position = "bottom")
  return(volcano_plot)
}

DESEqResults <- function(dge.obj, logfc, FDR) {
  ##This function adds a new column (T or F) according to the FDR and LFC of each gene in edgeR list of DEG
  ##dge.obj: List with DEG
  volc <- dge.obj %>%
    mutate(condition = ifelse((dge.obj$log2FoldChange > logfc) & (dge.obj$padj < FDR), "Over-expressed",
                              ifelse((dge.obj$log2FoldChange < -logfc) & (dge.obj$padj < FDR), "Sub-expressed",
                                     ifelse((dge.obj$log2FoldChange > logfc) & (dge.obj$padj > FDR), "NS",
                                            ifelse((dge.obj$log2FoldChange < -logfc) & (dge.obj$padj > FDR), "NS",
                                                   ifelse((dge.obj$log2FoldChange < logfc) & (dge.obj$padj > FDR), "NS", "NS")))))) %>%
    drop_na()
  return(volc)
}

volcano_DESeq2 <- function(volc, lfc, FDR) {
  ##The function returns the volcano plot of the dataset
  ##volcano_data: Data frame containing the data to create the volcano plot
  volcano_plot <- ggplot(volc)+
    geom_point(aes(x = log2FoldChange, y = -log10(padj), color = condition))+
    scale_color_manual(name = "Condition",
                       labels = paste(c("NS", "Over-expressed", "Sub-expressed"), c(sum(volc$condition == "NS"), sum(volc$condition == "Over-expressed"), sum(volc$condition == "Sub-expressed"))),
                       values = c("#6e6d6e","#d84b47","#66c343"))+
    geom_vline(aes(xintercept = lfc), linetype = "dashed")+
    geom_vline(aes(xintercept = -lfc), linetype = "dashed")+
    geom_hline(aes(yintercept = -log10(FDR)), linetype = "dashed")+
    theme_set(theme_bw())+
    theme(plot.title = element_text(face = "bold", size = 18), 
          axis.title = element_text(size = 18),
          legend.title = element_text(face = "bold", size = 15),
          legend.text = element_text(size = 15), 
          legend.position = "bottom")
  return(volcano_plot)
}
