#!/bin/sh
# properties = {"type": "single", "rule": "tophat_map", "local": false, "input": ["data/processedData/trimmed_reads/V_D2_Rep_3_R1_paired.fastq.gz", "data/processedData/trimmed_reads/V_D2_Rep_3_R2_paired.fastq.gz"], "output": ["data/processedData/aligned_reads/tophat_output/V_D2_Rep_3/"], "wildcards": {"sample": "V_D2_Rep_3"}, "params": {}, "log": ["logs/V_D2_Rep_3.map_and_bam.log"], "threads": 1, "resources": {}, "jobid": 154, "cluster": {}}
cd /gpfs/scratch/cah422/stephensiTranscriptome-november2018-mayvInf && \
/storage/home/cah422/.local/.bin/miniconda3/bin/python \
-m snakemake data/processedData/aligned_reads/tophat_output/V_D2_Rep_3/ logs/V_D2_Rep_3.map_and_bam.log --snakefile /gpfs/scratch/cah422/stephensiTranscriptome-november2018-mayvInf/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /gpfs/scratch/cah422/stephensiTranscriptome-november2018-mayvInf/.snakemake/tmp.h4yzfsdn data/processedData/trimmed_reads/V_D2_Rep_3_R1_paired.fastq.gz data/processedData/trimmed_reads/V_D2_Rep_3_R2_paired.fastq.gz --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --no-hooks --nolock --timestamp --mode 2  --allowed-rules tophat_map  && touch "/gpfs/scratch/cah422/stephensiTranscriptome-november2018-mayvInf/.snakemake/tmp.h4yzfsdn/154.jobfinished" || (touch "/gpfs/scratch/cah422/stephensiTranscriptome-november2018-mayvInf/.snakemake/tmp.h4yzfsdn/154.jobfailed"; exit 1)

