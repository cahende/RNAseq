#!/bin/sh
# properties = {"type": "single", "rule": "tophat_map", "local": false, "input": ["data/processedData/trimmed_reads/C_D2_Rep_2_R1_paired.fastq.gz", "data/processedData/trimmed_reads/C_D2_Rep_2_R2_paired.fastq.gz"], "output": ["data/processedData/aligned_reads/tophat_output/C_D2_Rep_2/"], "wildcards": {"sample": "C_D2_Rep_2"}, "params": {}, "log": ["logs/C_D2_Rep_2.map_and_bam.log"], "threads": 1, "resources": {}, "jobid": 155, "cluster": {}}
cd /gpfs/scratch/cah422/stephensiTranscriptome-november2018-mayvInf && \
/storage/home/cah422/.local/.bin/miniconda3/bin/python \
-m snakemake data/processedData/aligned_reads/tophat_output/C_D2_Rep_2/ logs/C_D2_Rep_2.map_and_bam.log --snakefile /gpfs/scratch/cah422/stephensiTranscriptome-november2018-mayvInf/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /gpfs/scratch/cah422/stephensiTranscriptome-november2018-mayvInf/.snakemake/tmp.h4yzfsdn data/processedData/trimmed_reads/C_D2_Rep_2_R1_paired.fastq.gz data/processedData/trimmed_reads/C_D2_Rep_2_R2_paired.fastq.gz --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --no-hooks --nolock --timestamp --mode 2  --allowed-rules tophat_map  && touch "/gpfs/scratch/cah422/stephensiTranscriptome-november2018-mayvInf/.snakemake/tmp.h4yzfsdn/155.jobfinished" || (touch "/gpfs/scratch/cah422/stephensiTranscriptome-november2018-mayvInf/.snakemake/tmp.h4yzfsdn/155.jobfailed"; exit 1)

