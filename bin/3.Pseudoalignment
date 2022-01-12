#!bin/bash
##########################################################################
#This script is designed for perform pseudoalignment and counting of reads
##########################################################################

##Create the output directory for storing salmon quants
mkdir -p ../results/salmon_quants

##Create the index
salmon index -t ../reference/Homo_sapiens.GRCh38.cdna.all.fa -i ../reference/hsa_v38_gencode

##Perform pseudoalignment
for i in {1..12}; do
R1="samp${i}_clean_R1.fastq.gz"
R2="samp${i}_clean_R2.fastq.gz"
printf "\n"
echo "Processing sample samp${i}"
salmon quant -i ../reference/hsa_v38_gencode -l A -1 ../data/trimmed/$R1 -2 ../data/trimmed/$R2 -p 8 --validateMappings -o ../results/salmon_quants/samp${i}_quant
done
