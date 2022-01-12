#!bin/bash
##########################################################################
#This script is designed for analyze the quality of raw reads
##########################################################################

##Create directories to store output files
mkdir -p ../results/fastqc_output
mkdir -p ../results/multiquc_output

##Run fastqc
fastqc ../data/fastq_files/*.fastq.gz -o ../results/fastqc_output

##Run multiqc
multiqc ../results/fastqc_output/*.zip -o ../results/multiquc_output