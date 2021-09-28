DEGResults <- function(qlf) {
  ##This function returns a dataframe with all DEG
  ##qlf: Object obatined from the generalized linear model
  qlf <- topTags(qlf, n = Inf)
  qlf <- as.data.frame(qlf) 
  return(qlf)
}

volcanoplotR <- function(dge.obj, logfc, p.adj, colup = "#d84b47", coldown = "#82bf56", colns = "#6e6d6e") {
  ##This function adds a new column (T or F) according to the FDR and LFC of each gene in edgeR list of DEG
  ##dge.obj: List with DEG
  ##logfC: logFC threshold used for the differential expression test
  ##p.adj: p.adj or FDR threshold to obtain significant genes
  ##Updated 15-abr-2021 by Rodolfo Chavez
  #Rename the lfc and padj columns from DESeq2
  #Color update 2-ago-2021 by Rodolfo Chavez
  #Code update 21-sep-2021 by Rodolfo Chavez
  if("log2FoldChange" %in% colnames(dge.obj) & "padj" %in% colnames(dge.obj)) {
    dge.obj <- dge.obj %>% dplyr::rename(logFC = log2FoldChange, FDR = padj)
  }
  #Add an extra column according to the logfc and p.adj thresholds
  dge.obj <- dge.obj %>%
    mutate(condition = ifelse((dge.obj$logFC > logfc) & (dge.obj$FDR < p.adj), "Over-expressed",
                              ifelse((dge.obj$logFC < -logfc) & (dge.obj$FDR < p.adj), "Sub-expressed",
                                     ifelse((dge.obj$logFC > logfc) & (dge.obj$FDR > p.adj), "NS",
                                            ifelse((dge.obj$logFC < -logfc) & (dge.obj$FDR > p.adj), "NS",
                                                   ifelse((dge.obj$logFC < logfc) & (dge.obj$FDR > p.adj), "NS", "NS"))))))
  ##Eliminate NA values
  if(TRUE %in% unique(is.na(dge.obj$FDR))) {
    dge.obj <- dge.obj %>% drop_na()
  }
  #Plot the results
  volcano_plot <- ggplot(dge.obj)+
    geom_point(aes(x = logFC, y = -log10(FDR), color = condition))+
    scale_color_manual(name = "Condition",
                       labels = paste(c("NS", "Over-expressed", "Sub-expressed"), c(sum(dge.obj$condition == "NS"), sum(dge.obj$condition == "Over-expressed"), sum(dge.obj$condition == "Sub-expressed"))),
                       values = c(colns, colup, coldown))+
    geom_vline(aes(xintercept = logfc), linetype = "dashed")+
    geom_vline(aes(xintercept = -logfc), linetype = "dashed")+
    geom_hline(aes(yintercept = -log10(p.adj)), linetype = "dashed")+
    theme(plot.title = element_text(face = "bold", size = 18), 
          axis.title = element_text(size = 18),
          legend.title = element_text(face = "bold", size = 15),
          legend.text = element_text(size = 15), 
          legend.position = "bottom")
  return(volcano_plot)
}

id_converter <- function(mart, input, attributes, filter) {
  ##The function return a data frame containing the columns of entrez id's for
  ##ensembl gene id or transcript gene id. Takes list of unique, shared or all
  ##significant DEGS
  #Mart: Data set obtained from BioMart
  #input: List of unique, shared or all DEG
  #attributes: columns to be included in the output
  #filter: Feature for filtering and match data
  #Get the list with the required attributes
  ##Updated 15-abr-2021 by Rodolfo Chavez
  ##Updated 24-ago-2021 by Rodolfo Chavez
  list_output <- getBM(attributes = attributes, #Attributes to obtain
                       #Name of the dataset
                       mart = mart, 
                       #Filter the list using ensembl_gene_id as input
                       filters = filter, 
                       #Input list
                       values = input)
  #Code the entrez id as character
  if("entrezgene_id" %in% colnames(list_output)) {
    list_output$entrezgene_id <- as.character(list_output$entrezgene_id)
  }
  #Fusion the input list with the list of new ids
  list_output <- full_join(input, list_output, by = filter)
  #Eliminate duplicated ids
  if(TRUE %in% duplicated(list_output[, 1])){
    list_output <- list_output[!duplicated(list_output[, 1]), ]
  }
  return(list_output)
}

barplot_GSEA <- function(fgsea_object, Head, Tail, col1 = "red", col2 = "gray") {
  ##This function creates a dotplot showing the most significant pathways splitted by condition
  #fgsea_object: Tibble containing the results of the significant pathways
  #Head: Number of pathways to be displayed in the first strings
  #Tail: Number of pathways to be displayed in the last strings 
  ##Updated by Rodolfo Chavez 5-ago-2021
  #Get the dotplot
  barplot <- ggplot(rbind(head(fgsea_object, Head), tail(fgsea_object, Tail)), aes(x = NES, y= reorder(pathway, NES), fill = Significance))+
    geom_col()+
    scale_fill_manual(values = c("Significant"=col1,
                                 "NS"=col2), 
                      #labels = c("Activated", "Supressed"), 
                      #name = "Condition")
    ) +
    labs(x = "NES", y = "Pathway")+
    theme(plot.title = element_text(face = "bold", size = 18, hjust = 1), 
          axis.title = element_text(face = "bold", size = 18),
          legend.title = element_text(face = "bold", size = 15), 
          legend.text = element_text(size = 15))+
    scale_y_discrete(labels = function(x) str_wrap(x, width = 50))
  return(barplot)
}
