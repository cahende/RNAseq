#!/usr/bin/env/python
shell.prefix("set -o pipefail; ")
shell.executable("/bin/bash")
shell.prefix("source ~/.bashrc; ")
configfile: "config.yaml"

rule all:
    input:
        expand("data/processedData/read_counts/{sample}.readCounts.txt", sample=config["SAMPLES"]),
        "data/processedData/top_hits/glmControlVsInfectedPValue.txt",
        "data/processedData/top_hits/glmControlVsInfectedLogFC.txt"

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
        "fastqc {input};"
        "trimmomatic PE -phred33 -trimlog {log} {input} {output} ILLUMINACLIP:{config[ADAPTERS]}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

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
        "fastqc {input};"
        "STAR --runThreadN 8 --genomeDir genomes --readFilesIn {input} --sjdbGTFfile {config[GENOME_ANNOTATION]} --sjdbOverhang 74 --outFileNamePrefix {output}"

rule sort_bam:
    input:
        "data/processedData/aligned_reads/star_output/{sample}/"
    output:
        "data/processedData/aligned_reads/sorted/{sample}.PE.star.sorted.bam"
    log: "logs/{sample}.sort_bam.log"
    shell:
        "fastqc {input};"
        "samtools view -Sb {input}/Aligned.out.sam | samtools sort -o {output} --threads 8"

rule index_bam:
    input:
        "data/processedData/aligned_reads/sorted/{sample}.PE.star.sorted.bam"
    output:
        "data/processedData/aligned_reads/sorted/{sample}.PE.star.sorted.bam.bai"
    log: "logs/{sample}.index_bam.log"
    shell:
        "fastqc {input};"
        "samtools index {input} {output}"

rule fix_mate_pairs:
    input:
        "data/processedData/aligned_reads/sorted/{sample}.PE.star.sorted.bam"
    output:
        "data/processedData/aligned_reads/fix_mate_pairs/{sample}.PE.star.sorted.fixed.bam"
    log: "logs/{sample}.fix_mate_pairs.log"
    shell:
        "fastqc {input};"
        "picard FixMateInformation INPUT={input} OUTPUT={output} SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true"

rule filter_mapped_and_paired_reads:
    input:
        "data/processedData/aligned_reads/fix_mate_pairs/{sample}.PE.star.sorted.fixed.bam"
    output:
        "data/processedData/aligned_reads/mapped_and_paired_filter/{sample}.PE.star.sorted.fixed.filtered.bam"
    log: "logs/{sample}.filter_mapped_and_paired_reads.log"
    shell:
        "fastqc {input};"
        "bamtools filter -isMapped true -in {input} -out {output}"

rule remove_duplicate_reads:
    input:
        "data/processedData/aligned_reads/mapped_and_paired_filter/{sample}.PE.star.sorted.fixed.filtered.bam"
    output:
        "data/processedData/aligned_reads/duplicate_removal/{sample}.PE.star.sorted.fixed.filtered.postdup.bam"
    log: "logs/{sample}.remove_duplicate_reads.log"
    shell:
        "fastqc {input};"
        "picard MarkDuplicates INPUT={input} OUTPUT={output} VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000 METRICS_FILE={log}"

rule add_read_groups:
    input:
        "data/processedData/aligned_reads/duplicate_removal/{sample}.PE.star.sorted.fixed.filtered.postdup.bam"
    output:
        "data/processedData/aligned_reads/read_group/{sample}.PE.star.sorted.fixed.filtered.postdup.RG.bam"
    log: "logs/{sample}.add_read_groups.log"
    shell:
        "fastqc {input};"
        "picard AddOrReplaceReadGroups INPUT={input} OUTPUT={output} RGLB={wildcards.sample}.PE RGPL=Illumina RGPU=Group1 RGSM={wildcards.sample}.PE"

rule quality_filter_reads:
    input:
        "data/processedData/aligned_reads/read_group/{sample}.PE.star.sorted.fixed.filtered.postdup.RG.bam"
    output:
        "data/processedData/aligned_reads/quality_filter/{sample}.PE.star.sorted.fixed.filtered.postdup.RG.passed.bam"
    log: "logs/{sample}.quality_filter_reads.log"
    shell:
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


