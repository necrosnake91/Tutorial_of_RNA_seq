## Tutorial for Differential Expression analysis

Hello everyone! 

Welcome to this GitHub repository. This repository was designed to provide you a wide view for the Differential Expression (DE) analysis using bioconductor packages **[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)** and **[DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html)**.


## Required packages

Before running the analysis, please install the following R packages:

* __tidyverse__
* __pheatmap__
* __edgeR__ ```BiocManager::install("edgeR")```
* __DESeq2__ ```BiocManager::install("DESeq2")```
* __PCAtools__ ```BiocManager::install("PCAtools")```
* __marray__ ```BiocManager::install("marray")```

Make sure that you have installed the last version of BioConductor. If not, run the next lines on R:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.12")
```

## Folder content

This is the distribution of folders and files in this repository:

```
+-- Tutorial for Differential Expression analysis
|		+--bin/
|		+--results/
```
### ```bin/```
This folder stores the required scripts to run the DE analysis in R:

* ```Differential_expression_analysis.R``` Script for performing DE analysis using edgeR package.
*  ```functions.R``` This script contains some functions to recodificate the results (edgeR or DESeq2) of DE genes into a dataframe and to create volcano plots.

### ```results/```
This folder will store all the plots, in .png format, created by running 
[```Differential_expression_analysis.R```](https://github.com/necrosnake91/Tutorial_of_RNA_seq/blob/main/bin/Differential_expression_analysis.R) script. Also, the count matrix ```counts.txt``` is stored in this folder.

Additionally, you will find ```Differential_expression_analysis_tutorial.html``` file which contains the presentation.

I suggest you take a read to the user's manual for each package in order to find useful information to perform DE analysis on your data. 