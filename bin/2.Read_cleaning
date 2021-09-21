#!bin/bash
##########################################################################
#This script performs trimming of illumina adapters
##########################################################################

##Create a folder to store trimmed data
mkdir -p ../data/trimmed

##Run cutadapt to trim illumina universal adapters
for i in {1..12} ; do
R1="samp${i}_R1.fastq.gz"
R2="samp${i}_R2.fastq.gz"
R1_o="samp${i}_clean_R1.fastq.gz"
R2_o="samp${i}_clean_R2.fastq.gz"
printf "\n"
echo "Processing sample samp${i}"
cutadapt -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -o ../data/trimmed/$R1_o -p ../data/trimmed/$R2_o ../data/fastq_files/$R1 ../data/fastq_files/$R2 
done
