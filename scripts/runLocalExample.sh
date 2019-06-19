#!/usr/bin/bash

# Download test data from SRA - Aedes aegypt transcrptome infected with Zika virus - control v infected at 2 days post infection (Hughes et al. 2017) 
esearch -db sra -query PRJNA399504 | efetch -format runinfo > PRJNA399504.csv 
mkdir data
mv data
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

trimmomatic PE -phred33 data/ctrl3/*1.fastq zikv1/*2.fastq data/zikv1/zikv1_R1_paired.fastq data/zikv1/zikv1_R1_unpaired.fastq data/zikv1/zikv1_R2_paired.fastq data/zikv1/zikv1_R2_unpaired.fastq ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# download reference genome (Aedes aegypti liverpool strain chromosomes) and index for STAR aligner
# also download GTF
mkdir genomes
wget -O genomes/AagL5.fa.gz https://www.vectorbase.org/download/aedes-aegypti-lvpagwgchromosomesaaegl5fagz
wget -O genomes/AagL5.gtf.gz https://www.vectorbase.org/download/aedes-aegypti-lvpagwgbasefeaturesaaegl52gtfgz
gunzip -d geomes/*.gz 

STAR --runMode genomeGenerate --genomeDir genomes --genomeFastaFiles genomes/AagL5.fa --sjdbGTFfile genomes/AagL5.gtf --sjdbOverhang 74

# align trimmed and paired reads to reference genome using STAR
STAR --genomeDir genomes --readFilesIn data/zikv1/zikv1_R1_paired.fastq data/zikv1/zikv1_R2_paired.fastq --sjdbGTFfile genomes/AagL5.gtf --sjdbOverhang 74 --outFileNamePrefix data/zikv1/starAlign/

STAR --genomeDir genomes --readFilesIn data/zikv2/zikv2_R1_paired.fastq data/zikv2/zikv2_R2_paired.fastq --sjdbGTFfile genomes/AagL5.gtf --sjdbOverhang 74 --outFileNamePrefix data/zikv2/starAlign/

STAR --genomeDir genomes --readFilesIn data/zikv3/zikv3_R1_paired.fastq data/zikv3/zikv3_R2_paired.fastq --sjdbGTFfile genomes/AagL5.gtf --sjdbOverhang 74 --outFileNamePrefix data/zikv3/starAlign/

STAR --genomeDir genomes --readFilesIn data/ctrl1/ctrl1_R1_paired.fastq data/ctrl1/ctrl1_R2_paired.fastq --sjdbGTFfile genomes/AagL5.gtf --sjdbOverhang 74 --outFileNamePrefix data/ctrl1/starAlign/

STAR --genomeDir genomes --readFilesIn data/ctrl2/ctrl2_R1_paired.fastq data/ctrl2/ctrl2_R2_paired.fastq --sjdbGTFfile genomes/AagL5.gtf --sjdbOverhang 74 --outFileNamePrefix data/ctrl2/starAlign/

STAR --genomeDir genomes --readFilesIn data/ctrl3/ctrl3_R1_paired.fastq data/ctrl3/ctrl3_R2_paired.fastq --sjdbGTFfile genomes/AagL5.gtf --sjdbOverhang 74 --outFileNamePrefix data/ctrl3/starAlign/

# sort and index STAR aligned output
samtools view -Sb zikv1/starAlign/Aligned.out.sam | samtools sort -o zikv1/zikv1.sort.bam 
samtools index zikv1/zikv1.sort.bam zikv1/zikv1.sort.bam.bai

samtools view -Sb zikv2/starAlign/Aligned.out.sam | samtools sort -o zikv2/zikv2.sort.bam
samtools index zikv2/zikv2.sort.bam zikv2/zikv2.sort.bam.bai

samtools view -Sb zikv3/starAlign/Aligned.out.sam | samtools sort -o zikv3/zikv3.sort.bam
samtools index zikv3/zikv3.sort.bam zikv3/zikv3.sort.bam.bai

samtools view -Sb ctrl1/starAlign/Aligned.out.sam | samtools sort -o ctrl1/ctrl1.sort.bam
samtools index ctrl1/ctrl1.sort.bam ctrl1/ctrl1.sort.bam.bai

samtools view -Sb ctrl2/starAlign/Aligned.out.sam | samtools sort -o ctrl2/ctrl2.sort.bam
samtools index ctrl2/ctrl2.sort.bam ctrl2/ctrl2.sort.bam.bai

samtools view -Sb ctrl3/starAlign/Aligned.out.sam | samtools sort -o ctrl3/ctrl3.sort.bam
samtools index ctrl3/ctrl3.sort.bam ctrl3/ctrl3.sort.bam.bai

# fix mate pairs
picard FixMateInformation INPUT=zikv1/zikv1.sort.bam OUTPUT=zikv1/zikv1.sort.fix.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

picard FixMateInformation INPUT=zikv2/zikv2.sort.bam OUTPUT=zikv2/zikv2.sort.fix.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

picard FixMateInformation INPUT=zikv3/zikv3.sort.bam OUTPUT=zikv3/zikv3.sort.fix.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

picard FixMateInformation INPUT=ctrl1/ctrl1.sort.bam OUTPUT=ctrl1/ctrl1.sort.fix.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

picard FixMateInformation INPUT=ctrl2/ctrl2.sort.bam OUTPUT=ctrl2/ctrl2.sort.fix.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

picard FixMateInformation INPUT=ctrl3/ctrl3.sort.bam OUTPUT=ctrl3/ctrl3.sort.fix.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

# filter mapped and paired reads
bamtools filter -isMapped true -in zikv1/zikv1.sort.fix.bam -out zikv1/zikv1.sort.fix.map.bam

bamtools filter -isMapped true -in zikv2/zikv2.sort.fix.bam -out zikv2/zikv2.sort.fix.map.bam

bamtools filter -isMapped true -in zikv3/zikv3.sort.fix.bam -out zikv3/zikv3.sort.fix.map.bam

bamtools filter -isMapped true -in ctrl1/ctrl1.sort.fix.bam -out ctrl1/ctrl1.sort.fix.map.bam

bamtools filter -isMapped true -in ctrl2/ctrl2.sort.fix.bam -out ctrl2/ctrl2.sort.fix.map.bam

bamtools filter -isMapped true -in ctrl3/ctrl3.sort.fix.bam -out ctrl3/ctrl3.sort.fix.map.bam

# remove duplicate reads
picard MarkDuplicates INPUT=zikv1/zikv1.sort.fix.map.bam OUTPUT=zikv1/zikv1.sort.fix.map.dup.bam VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000

picard MarkDuplicates INPUT=zikv2/zikv2.sort.fix.map.bam OUTPUT=zikv2/zikv2.sort.fix.map.dup.bam VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000

picard MarkDuplicates INPUT=zikv3/zikv3.sort.fix.map.bam OUTPUT=zikv3/zikv3.sort.fix.map.dup.bam VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000

picard MarkDuplicates INPUT=ctrl1/ctrl1.sort.fix.map.bam OUTPUT=ctrl1/ctrl1.sort.fix.map.dup.bam VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000

picard MarkDuplicates INPUT=ctrl2/ctrl2.sort.fix.map.bam OUTPUT=ctrl2/ctrl2.sort.fix.map.dup.bam VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000

picard MarkDuplicates INPUT=ctrl3/ctrl3.sort.fix.map.bam OUTPUT=ctrl3/ctrl3.sort.fix.map.dup.bam VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000

# filter for mapping quality
bamtools filter -mapQuality '>=20' -length '100' -in zikv1/zikv1.sort.fix.map.dup.bam -out zikv1/zikv1.final.bam 

bamtools filter -mapQuality '>=20' -length '100' -in zikv2/zikv2.sort.fix.map.dup.bam -out zikv2/zikv2.final.bam

bamtools filter -mapQuality '>=20' -length '100' -in zikv3/zikv3.sort.fix.map.dup.bam -out zikv3/zikv3.final.bam

bamtools filter -mapQuality '>=20' -length '100' -in ctrl1/ctrl1.sort.fix.map.dup.bam -out ctrl1/ctrl1.final.bam

bamtools filter -mapQuality '>=20' -length '100' -in ctrl2/ctrl2.sort.fix.map.dup.bam -out ctrl2/ctrl2.final.bam

bamtools filter -mapQuality '>=20' -length '100' -in ctrl3/ctrl3.sort.fix.map.dup.bam -out ctrl3/ctrl3.final.bam

# perform differential expression analysis using edgeR and rSubread
Rscript scripts/rScripts/localExample.R
