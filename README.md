## Tutorial for Differential Expression analysis

Hello everyone! 

Welcome to this GitHub repository. This repository was designed to provide you a wide view for the Differential Expression (DE) analysis using bioconductor packages **[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)** and **[DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html)**.

I suggest you take a read to the user's manual for each package in order to find useful information regarding to the critical parameters to perform DE analysis on your data. 

## Required packages

Before running the analysis, please install the next packages in R:

* __tidyverse__
* __pheatmap__
* __edgeR__ ```BiocManager::install("edgeR")```
* __DESeq2__ ```BiocManager::install("DESeq2")```
* __PCAtools__ ```BiocManager::install("PCAtools")```
* __marray__ ```BiocManager::install("marray")```

## Folder content

In this repository you will find the next folders and files:

```
+-- Tutorial for Differential Expression analysis
|		+--bin/
|		+--results/
```
###```bin/```
This folder contains the required scripts to run the DE analysis in R:

* ```Differential_expression_analysis.R``` Script for performing DE analysis using edgeR package.
*  ```funciones.R``` This script contains some functions to recodificate the results (edgeR or DESeq2) obtained from DE analysis into a dataframe and to create volcano plots.

###```results/```
This folder will store all the plots created by running 
```Differential_expression_analysis.R``` script in .png format. Please, after cloning this repository, save the provided count matrix ```counts.txt``` in this folder.