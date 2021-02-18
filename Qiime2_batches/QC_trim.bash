#This script needs be run under data directory, output is in QC_data at upper folder
#The file names must be in form of name_R1.fastq.gz name_R2.fastq.gz
#How to use: bash step1_QC_pandaseq.bash

mkdir ../QC_data

ls *gz|cut -f 1 -d "_" |uniq|while read line; do mv "$line"*_R1_001.fastq.gz "$line"_R1.fastq.gz;done

ls *gz|cut -f 1 -d "_" |uniq|while read line; do mv "$line"*_R2_001.fastq.gz "$line"_R2.fastq.gz;done

ls *gz|cut -f 1 -d "_" |uniq|while read line; do java -jar /opt/software/Trimmomatic/0.36-Java-1.8.0_92/trimmomatic-0.36.jar PE -phred33 "$line"_R1.fastq.gz "$line"_R2.fastq.gz ../QC_data/"$line"_R1.fastq.gz ../QC_data/"$line".qcup_R1.fastq.gz ../QC_data/"$line"_R2.fastq.gz ../QC_data/"$line".qcup_R2.fastq.gz ILLUMINACLIP:/opt/software/Trimmomatic/0.36-Java-1.8.0_92/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36; done

rm ../QC_data/*qcup*

for i in ../QC_data/*gz; do zcat $i|wc -l >> ../QC_data/seq_num.txt; done
