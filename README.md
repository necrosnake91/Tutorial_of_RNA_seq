## Tutorial for Differential Expression analysis

Hello everyone! 

Welcome to this GitHub repository. This repository was designed to provide you a wide view for the Differential Expression (DE) analysis using bioconductor packages **[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)** and **[DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html)**.


## Required packages

Before running the analysis, please install the following R packages from CRAN:

* tidyverse
* RColorBrewer
* pheatmap
* here

Also, make sure that you have installed the last version of BioConductor. If not, run the next lines on R:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.13")
```
and install:

* tximeta
* DESeq2
* PCAtools
* apeglm
* Glimma
* clusterProfiler
* enrichplot
* biomaRt
* fgsea

## Folder content

This is the distribution of folders and files in this repository:

```
+-- Tutorial for Differential Expression analysis
|		+--bin/
|		+--data/
|		+--Output/
|		+--slides/
|		+--reference/
```
### `bin/`
This folder stores the required scripts to run the DE analysis in R:

`1. QC_analysis`. Script to perform quality control of raw reads using **FastQC** and **MultiQC**.

`2. Read_cleaning`. Script for clean the reads and remove the Iluminna universal adapter sequence from raw reads using **Cutadapt**.

`3. Pseudoalignment`. This script contains the code for performing pseudoalignment and abundance estimation using **Salmon**

`4. Differential_expression_analysis.R`. This script performs differential expression testing using the algorithm of **DESeq2**.

`functions.R`. Script containing useful functions.

### `data`
Folder for storing the raw counts in `fastq.gz` format. After getting the raw counts, please copy them in this folder

### `Output`
This folder stores the results derived from:

* FastQC
* MultiQC
* Salmon_quants

### `slides`
This folder contains the Rmd and html files for the slides and the cheatsheet. Follow this [link](https://github.com/necrosnake91/Tutorial_of_RNA_seq/blob/main/slides/Slides.html#1) to visualize the Slides. The cheatsheet is ubicated [here](https://github.com/necrosnake91/Tutorial_of_RNA_seq/blob/main/slides/Differential_expression_analysis_tutorial.html).

### `reference`
The reference transcriptome downloaded from the [GENECODE](https://www.gencodegenes.org/human/) site is stored here in `fasta` format.

I suggest you take a read to the user's manual for each package in order to find useful information to perform DE analysis on your data. 