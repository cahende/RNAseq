#!/usr/bin/bash

# Download test data from SRA - Aedes aegypt transcrptome infected with Zika virus - control v infected at 2 days post infection (Hughes et al. 2017) 
esearch -db sra -query PRJNA399504 | efetch -format runinfo > PRJNA399504.csv 
fastq-dump --split-files -O zikvInf1 -X 10000 SRR5955021 
fastq-dump --split-files -O zikvInf2 -X 10000 SRR5955020
fastq-dump --split-files -O zikvInf3 -X 10000 SRR5955019
fastq-dump --split-files -O ctrl1 -X 10000 SRR5955018
fastq-dump --split-files -O ctrl2 -X 10000 SRR5955017
fastq-dump --split-files -O ctrl3 -X 10000 SRR5955016

# trim adapter sequences and low quality bases from reads using trimmomatic in paired-end mode
trimmomatic PE -phred33 zikv1/*1.fastq zikv1/*2.fastq zikv1/zikv1_R1_paired.fastq zikv1/zikv1_R1_unpaired.fastq zikv1/zikv1_R2_paired.fastq zikv1/zikv1_R2_unpaired.fastq ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

trimmomatic PE -phred33 zikv2/*1.fastq zikv2/*2.fastq zikv2/zikv2_R1_paired.fastq zikv2/zikv2_R1_unpaired.fastq zikv2/zikv2_R2_paired.fastq zikv2/zikv2_R2_unpaired.fastq ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

trimmomatic PE -phred33 zikv3/*1.fastq zikv3/*2.fastq zikv3/zikv3_R1_paired.fastq zikv3/zikv3_R1_unpaired.fastq zikv3/zikv3_R2_paired.fastq zikv3/zikv3_R2_unpaired.fastq ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

trimmomatic PE -phred33 ctrl1/*1.fastq ctrl1/*2.fastq ctrl1/ctrl1_R1_paired.fastq ctrl1/ctrl1_R1_unpaired.fastq ctrl1/ctrl1_R2_paired.fastq ctrl1/ctrl1_R2_unpaired.fastq ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

trimmomatic PE -phred33 ctrl2/*1.fastq ctrl2/*2.fastq ctrl2/ctrl2_R1_paired.fastq ctrl2/ctrl2_R1_unpaired.fastq ctrl2/ctrl2_R2_paired.fastq ctrl2/ctrl2_R2_unpaired.fastq ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

trimmomatic PE -phred33 ctrl3/*1.fastq zikv1/*2.fastq zikv1/zikv1_R1_paired.fastq zikv1/zikv1_R1_unpaired.fastq zikv1/zikv1_R2_paired.fastq zikv1/zikv1_R2_unpaired.fastq ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

rule trim_reads:
    input:
        expand("data/processedData/trimmed_reads/{{sample}}_{read}_paired.fastq.gz", read=["R1", "R2"])
    output:
        "data/processedData/trimmed_reads/{sample}_R1_paired.fastq.gz",
        "data/processedData/trimmed_reads/{sample}_R1_unpaired.fastq.gz",
        "data/processedData/trimmed_reads/{sample}_R2_paired.fastq.gz",
        "data/processedData/trimmed_reads/{sample}_R2_unpaired.fastq.gz"
    log: "logs/trim/{sample}.trim.log"
    shell:
        "module load trimmomatic fastqc;"
        "fastqc {input};"
        "java -jar {config[TRIMMOMATIC]} PE -phred33 -trimlog {log} {input} {output} ILLUMINACLIP:{config[ADAPTERS]}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

rule decompress_trimmed_reads:
    input:
        R1="data/processedData/trimmed_reads/{sample}_R1_paired.fastq.gz",
        R2="data/processedData/trimmed_reads/{sample}_R2_paired.fastq.gz"
    output:
        R1="data/processedData/trimmed_reads/{sample}_R1_paired.fastq",
        R2="data/processedData/trimmed_reads/{sample}_R2_paired.fastq"
    shell:
        "gunzip -dc {input.R1} > {output.R1};"
        "gunzip -dc {input.R2} > {output.R2}"

rule star_map:
    input:
        expand("data/processedData/trimmed_reads/{{sample}}_{read}_paired.fastq", read=["R1", "R2"])
    output:
        directory("data/processedData/aligned_reads/star_output/{sample}/")
    log: "logs/{sample}.map_and_bam.log"
    shell:
        "module load star fastqc;"
        "fastqc {input};"
        "STAR --runThreadN 8 --genomeDir genomes --readFilesIn {input} --sjdbGTFfile {config[GENOME_ANNOTATION]} --sjdbOverhang 74 --outFileNamePrefix {output}"

rule sort_bam:
    input:
        "data/processedData/aligned_reads/star_output/{sample}/"
    output:
        "data/processedData/aligned_reads/sorted/{sample}.PE.star.sorted.bam"
    log: "logs/{sample}.sort_bam.log"
    shell:
        "module load samtools fastqc;"
        "fastqc {input};"
        "samtools view -Sb {input}/Aligned.out.sam | samtools sort -o {output} --threads 8"

rule index_bam:
    input:
        "data/processedData/aligned_reads/sorted/{sample}.PE.star.sorted.bam"
    output:
        "data/processedData/aligned_reads/sorted/{sample}.PE.star.sorted.bam.bai"
    log: "logs/{sample}.index_bam.log"
    shell:
        "module load samtools fastqc;"
        "fastqc {input};"
        "samtools index {input} {output}"

rule fix_mate_pairs:
    input:
        "data/processedData/aligned_reads/sorted/{sample}.PE.star.sorted.bam"
    output:
        "data/processedData/aligned_reads/fix_mate_pairs/{sample}.PE.star.sorted.fixed.bam"
    log: "logs/{sample}.fix_mate_pairs.log"
    shell:
        "module load picard fastqc;"
        "fastqc {input};"
        "java -jar {config[PICARD]} FixMateInformation INPUT={input} OUTPUT={output} SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true"

rule filter_mapped_and_paired_reads:
    input:
        "data/processedData/aligned_reads/fix_mate_pairs/{sample}.PE.star.sorted.fixed.bam"
    output:
        "data/processedData/aligned_reads/mapped_and_paired_filter/{sample}.PE.star.sorted.fixed.filtered.bam"
    log: "logs/{sample}.filter_mapped_and_paired_reads.log"
    shell:
        "module load bamtools fastqc;"
        "fastqc {input};"
        "bamtools filter -isMapped true -in {input} -out {output}"

rule remove_duplicate_reads:
    input:
        "data/processedData/aligned_reads/mapped_and_paired_filter/{sample}.PE.star.sorted.fixed.filtered.bam"
    output:
        "data/processedData/aligned_reads/duplicate_removal/{sample}.PE.star.sorted.fixed.filtered.postdup.bam"
    log: "logs/{sample}.remove_duplicate_reads.log"
    shell:
        "module load picard fastqc;"
        "fastqc {input};"
        "java -jar {config[PICARD]} MarkDuplicates INPUT={input} OUTPUT={output} VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000 METRICS_FILE={log}"

rule add_read_groups:
    input:
        "data/processedData/aligned_reads/duplicate_removal/{sample}.PE.star.sorted.fixed.filtered.postdup.bam"
    output:
        "data/processedData/aligned_reads/read_group/{sample}.PE.star.sorted.fixed.filtered.postdup.RG.bam"
    log: "logs/{sample}.add_read_groups.log"
    shell:
        "module load picard fastqc;"
        "fastqc {input};"
        "java -jar {config[PICARD]} AddOrReplaceReadGroups INPUT={input} OUTPUT={output} RGLB={wildcards.sample}.PE RGPL=Illumina RGPU=Group1 RGSM={wildcards.sample}.PE"

rule quality_filter_reads:
    input:
        "data/processedData/aligned_reads/read_group/{sample}.PE.star.sorted.fixed.filtered.postdup.RG.bam"
    output:
        "data/processedData/aligned_reads/quality_filter/{sample}.PE.star.sorted.fixed.filtered.postdup.RG.passed.bam"
    log: "logs/{sample}.quality_filter_reads.log"
    shell:
        "module load bamtools fastqc;"
        "fastqc {input};"
        "bamtools filter -mapQuality '>=20' -length '75' -in {input} -out {output};"
        "fastqc {output}"

rule rSubread:
    input:
        "data/processedData/aligned_reads/quality_filter/{sample}.PE.star.sorted.fixed.filtered.postdup.RG.passed.bam"
    output:
        "data/processedData/read_counts/{sample}.readCounts.txt"
    log: "logs/{sample}.rSubread.log"
    script:
        "scripts/rScripts/rSubreadFeatureCounts.R"

rule gene_level_differential_expression:
    input:
        expand("data/processedData/read_counts/{sample}.readCounts.txt", sample=config["SAMPLES"])
    output:
        "data/processedData/top_hits/glmControlVsInfectedPValue.txt",
        "data/processedData/top_hits/glmControlVsInfectedLogFC.txt"
    log: "logs/edgeR.log"
    script:
        "scripts/rScripts/edgeR.R"

