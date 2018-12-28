#!/bin/sh
# properties = {"type": "single", "rule": "tophat_map", "local": false, "input": ["data/processedData/trimmed_reads/C_D14_Rep_1_R1_paired.fastq.gz", "data/processedData/trimmed_reads/C_D14_Rep_1_R2_paired.fastq.gz"], "output": ["data/processedData/aligned_reads/tophat_output/C_D14_Rep_1/"], "wildcards": {"sample": "C_D14_Rep_1"}, "params": {}, "log": ["logs/C_D14_Rep_1.map_and_bam.log"], "threads": 1, "resources": {}, "jobid": 160, "cluster": {}}
cd /gpfs/scratch/cah422/stephensiTranscriptome-november2018-mayvInf && \
/storage/home/cah422/.local/.bin/miniconda3/bin/python \
-m snakemake data/processedData/aligned_reads/tophat_output/C_D14_Rep_1/ logs/C_D14_Rep_1.map_and_bam.log --snakefile /gpfs/scratch/cah422/stephensiTranscriptome-november2018-mayvInf/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /gpfs/scratch/cah422/stephensiTranscriptome-november2018-mayvInf/.snakemake/tmp.h4yzfsdn data/processedData/trimmed_reads/C_D14_Rep_1_R1_paired.fastq.gz data/processedData/trimmed_reads/C_D14_Rep_1_R2_paired.fastq.gz --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --no-hooks --nolock --timestamp --mode 2  --allowed-rules tophat_map  && touch "/gpfs/scratch/cah422/stephensiTranscriptome-november2018-mayvInf/.snakemake/tmp.h4yzfsdn/160.jobfinished" || (touch "/gpfs/scratch/cah422/stephensiTranscriptome-november2018-mayvInf/.snakemake/tmp.h4yzfsdn/160.jobfailed"; exit 1)

