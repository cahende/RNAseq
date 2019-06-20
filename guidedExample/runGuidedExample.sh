#!/bin/bash

# RNAsequencing pipeline example 
# tools required: NCBI utilities, SRA toolkit, trimmomatic, STAR splice-aware aligner, samtools, bamtools, picard, rSubread, and edgeR

# download test data from SRA - Aedes aegypt transcrptome infected with Zika virus - control v infected at 2 days post infection (Hughes et al. 2017) 
esearch -db sra -query PRJNA399504 | efetch -format runinfo > PRJNA399504.csv 
mkdir data
cd data
fastq-dump --split-files -O zikv1 -X 10000 SRR5955021 
fastq-dump --split-files -O zikv2 -X 10000 SRR5955020
fastq-dump --split-files -O zikv3 -X 10000 SRR5955019
fastq-dump --split-files -O ctrl1 -X 10000 SRR5955018
fastq-dump --split-files -O ctrl2 -X 10000 SRR5955017
fastq-dump --split-files -O ctrl3 -X 10000 SRR5955016
cd ..

# trim adapter sequences and low quality bases from reads using trimmomatic in paired-end mode
trimmomatic PE -phred33 data/zikv1/*1.fastq data/zikv1/*2.fastq data/zikv1/zikv1_R1_paired.fastq data/zikv1/zikv1_R1_unpaired.fastq data/zikv1/zikv1_R2_paired.fastq data/zikv1/zikv1_R2_unpaired.fastq ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

trimmomatic PE -phred33 data/zikv2/*1.fastq data/zikv2/*2.fastq data/zikv2/zikv2_R1_paired.fastq data/zikv2/zikv2_R1_unpaired.fastq data/zikv2/zikv2_R2_paired.fastq data/zikv2/zikv2_R2_unpaired.fastq ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

trimmomatic PE -phred33 data/zikv3/*1.fastq data/zikv3/*2.fastq data/zikv3/zikv3_R1_paired.fastq data/zikv3/zikv3_R1_unpaired.fastq data/zikv3/zikv3_R2_paired.fastq data/zikv3/zikv3_R2_unpaired.fastq ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

trimmomatic PE -phred33 data/ctrl1/*1.fastq data/ctrl1/*2.fastq data/ctrl1/ctrl1_R1_paired.fastq data/ctrl1/ctrl1_R1_unpaired.fastq data/ctrl1/ctrl1_R2_paired.fastq data/ctrl1/ctrl1_R2_unpaired.fastq ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

trimmomatic PE -phred33 data/ctrl2/*1.fastq data/ctrl2/*2.fastq data/ctrl2/ctrl2_R1_paired.fastq data/ctrl2/ctrl2_R1_unpaired.fastq data/ctrl2/ctrl2_R2_paired.fastq data/ctrl2/ctrl2_R2_unpaired.fastq ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

trimmomatic PE -phred33 data/ctrl3/*1.fastq data/ctrl3/*2.fastq data/ctrl3/ctrl3_R1_paired.fastq data/ctrl3/ctrl3_R1_unpaired.fastq data/ctrl3/ctrl3_R2_paired.fastq data/ctrl3/ctrl3_R2_unpaired.fastq ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# download reference genome (Aedes aegypti liverpool strain chromosomes) and index for STAR aligner
# also download GTF
mkdir genomes
wget -O genomes/AagL5.fa.gz https://www.vectorbase.org/download/aedes-aegypti-lvpagwgchromosomesaaegl5fagz
wget -O genomes/AagL5.gtf.gz https://www.vectorbase.org/download/aedes-aegypti-lvpagwgbasefeaturesaaegl52gtfgz
gunzip -d genomes/*.gz 

STAR --runMode genomeGenerate --genomeDir genomes --genomeFastaFiles genomes/AagL5.fa --sjdbGTFfile genomes/AagL5.gtf --sjdbOverhang 74

# align trimmed and paired reads to reference genome using STAR
STAR --genomeDir genomes --readFilesIn data/zikv1/zikv1_R1_paired.fastq data/zikv1/zikv1_R2_paired.fastq --sjdbGTFfile genomes/AagL5.gtf --sjdbOverhang 74 --outFileNamePrefix data/zikv1/

STAR --genomeDir genomes --readFilesIn data/zikv2/zikv2_R1_paired.fastq data/zikv2/zikv2_R2_paired.fastq --sjdbGTFfile genomes/AagL5.gtf --sjdbOverhang 74 --outFileNamePrefix data/zikv2/

STAR --genomeDir genomes --readFilesIn data/zikv3/zikv3_R1_paired.fastq data/zikv3/zikv3_R2_paired.fastq --sjdbGTFfile genomes/AagL5.gtf --sjdbOverhang 74 --outFileNamePrefix data/zikv3/

STAR --genomeDir genomes --readFilesIn data/ctrl1/ctrl1_R1_paired.fastq data/ctrl1/ctrl1_R2_paired.fastq --sjdbGTFfile genomes/AagL5.gtf --sjdbOverhang 74 --outFileNamePrefix data/ctrl1/

STAR --genomeDir genomes --readFilesIn data/ctrl2/ctrl2_R1_paired.fastq data/ctrl2/ctrl2_R2_paired.fastq --sjdbGTFfile genomes/AagL5.gtf --sjdbOverhang 74 --outFileNamePrefix data/ctrl2/

STAR --genomeDir genomes --readFilesIn data/ctrl3/ctrl3_R1_paired.fastq data/ctrl3/ctrl3_R2_paired.fastq --sjdbGTFfile genomes/AagL5.gtf --sjdbOverhang 74 --outFileNamePrefix data/ctrl3/

# sort and index STAR aligned output
samtools view -Sb data/zikv1/Aligned.out.sam | samtools sort -o data/zikv1/zikv1.sort.bam 
samtools index data/zikv1/zikv1.sort.bam data/zikv1/zikv1.sort.bam.bai

samtools view -Sb data/zikv2/Aligned.out.sam | samtools sort -o data/zikv2/zikv2.sort.bam
samtools index data/zikv2/zikv2.sort.bam data/zikv2/zikv2.sort.bam.bai

samtools view -Sb data/zikv3/Aligned.out.sam | samtools sort -o data/zikv3/zikv3.sort.bam
samtools index data/zikv3/zikv3.sort.bam data/zikv3/zikv3.sort.bam.bai

samtools view -Sb data/ctrl1/Aligned.out.sam | samtools sort -o data/ctrl1/ctrl1.sort.bam
samtools index data/ctrl1/ctrl1.sort.bam data/ctrl1/ctrl1.sort.bam.bai

samtools view -Sb data/ctrl2/Aligned.out.sam | samtools sort -o data/ctrl2/ctrl2.sort.bam
samtools index data/ctrl2/ctrl2.sort.bam data/ctrl2/ctrl2.sort.bam.bai

samtools view -Sb data/ctrl3/Aligned.out.sam | samtools sort -o data/ctrl3/ctrl3.sort.bam
samtools index data/ctrl3/ctrl3.sort.bam data/ctrl3/ctrl3.sort.bam.bai

# fix mate pairs
picard FixMateInformation INPUT=data/zikv1/zikv1.sort.bam OUTPUT=data/zikv1/zikv1.sort.fix.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

picard FixMateInformation INPUT=data/zikv2/zikv2.sort.bam OUTPUT=data/zikv2/zikv2.sort.fix.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

picard FixMateInformation INPUT=data/zikv3/zikv3.sort.bam OUTPUT=data/zikv3/zikv3.sort.fix.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

picard FixMateInformation INPUT=data/ctrl1/ctrl1.sort.bam OUTPUT=data/ctrl1/ctrl1.sort.fix.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

picard FixMateInformation INPUT=data/ctrl2/ctrl2.sort.bam OUTPUT=data/ctrl2/ctrl2.sort.fix.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

picard FixMateInformation INPUT=data/ctrl3/ctrl3.sort.bam OUTPUT=data/ctrl3/ctrl3.sort.fix.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

# filter mapped and paired reads
bamtools filter -isMapped true -in data/zikv1/zikv1.sort.fix.bam -out data/zikv1/zikv1.sort.fix.map.bam

bamtools filter -isMapped true -in data/zikv2/zikv2.sort.fix.bam -out data/zikv2/zikv2.sort.fix.map.bam

bamtools filter -isMapped true -in data/zikv3/zikv3.sort.fix.bam -out data/zikv3/zikv3.sort.fix.map.bam

bamtools filter -isMapped true -in data/ctrl1/ctrl1.sort.fix.bam -out data/ctrl1/ctrl1.sort.fix.map.bam

bamtools filter -isMapped true -in data/ctrl2/ctrl2.sort.fix.bam -out data/ctrl2/ctrl2.sort.fix.map.bam

bamtools filter -isMapped true -in data/ctrl3/ctrl3.sort.fix.bam -out data/ctrl3/ctrl3.sort.fix.map.bam

# remove duplicate reads
picard MarkDuplicates INPUT=data/zikv1/zikv1.sort.fix.map.bam OUTPUT=data/zikv1/zikv1.sort.fix.map.dup.bam VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000 METRICS_FILE=data/zikv1/zikv1.duplog.txt

picard MarkDuplicates INPUT=data/zikv2/zikv2.sort.fix.map.bam OUTPUT=data/zikv2/zikv2.sort.fix.map.dup.bam VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000 METRICS_FILE=data/zikv2/zikv2.duplog.txt

picard MarkDuplicates INPUT=data/zikv3/zikv3.sort.fix.map.bam OUTPUT=data/zikv3/zikv3.sort.fix.map.dup.bam VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000 METRICS_FILE=data/zikv3/zikv3.duplog.txt

picard MarkDuplicates INPUT=data/ctrl1/ctrl1.sort.fix.map.bam OUTPUT=data/ctrl1/ctrl1.sort.fix.map.dup.bam VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000 METRICS_FILE=data/ctrl1/ctrl1.duplog.txt

picard MarkDuplicates INPUT=data/ctrl2/ctrl2.sort.fix.map.bam OUTPUT=data/ctrl2/ctrl2.sort.fix.map.dup.bam VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000 METRICS_FILE=data/ctrl2/ctrl2.duplog.txt

picard MarkDuplicates INPUT=data/ctrl3/ctrl3.sort.fix.map.bam OUTPUT=data/ctrl3/ctrl3.sort.fix.map.dup.bam VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000 METRICS_FILE=data/ctrl3/ctrl3.duplog.txt

# filter for mapping quality
bamtools filter -mapQuality '>=20' -length '50' -in data/zikv1/zikv1.sort.fix.map.dup.bam -out data/zikv1/zikv1.final.bam 

bamtools filter -mapQuality '>=20' -length '50' -in data/zikv2/zikv2.sort.fix.map.dup.bam -out data/zikv2/zikv2.final.bam

bamtools filter -mapQuality '>=20' -length '50' -in data/zikv3/zikv3.sort.fix.map.dup.bam -out data/zikv3/zikv3.final.bam

bamtools filter -mapQuality '>=20' -length '50' -in data/ctrl1/ctrl1.sort.fix.map.dup.bam -out data/ctrl1/ctrl1.final.bam

bamtools filter -mapQuality '>=20' -length '50' -in data/ctrl2/ctrl2.sort.fix.map.dup.bam -out data/ctrl2/ctrl2.final.bam

bamtools filter -mapQuality '>=20' -length '50' -in data/ctrl3/ctrl3.sort.fix.map.dup.bam -out data/ctrl3/ctrl3.final.bam

# perform differential expression analysis using edgeR and rSubread
Rscript localExample.R
