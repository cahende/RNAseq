#!/bin/sh
# properties = {"type": "single", "rule": "rSubread", "local": false, "input": ["data/processedData/aligned_reads/quality_filter/C_D14_Rep_1.PE.star.sorted.fixed.filtered.postdup.RG.passed.bam"], "output": ["data/processedData/read_counts/C_D14_Rep_1.readCounds.txt"], "wildcards": {"sample": "C_D14_Rep_1"}, "params": {}, "log": ["logs/C_D14_Rep_1.quality_filter_reads.log"], "threads": 1, "resources": {}, "jobid": 17, "cluster": {}}
cd /gpfs/scratch/cah422/stephensiTranscriptome-november2018-mayvInf && \
/storage/home/cah422/.local/.bin/miniconda3/bin/python \
-m snakemake data/processedData/read_counts/C_D14_Rep_1.readCounds.txt logs/C_D14_Rep_1.quality_filter_reads.log --snakefile /gpfs/scratch/cah422/stephensiTranscriptome-november2018-mayvInf/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /gpfs/scratch/cah422/stephensiTranscriptome-november2018-mayvInf/.snakemake/tmp.e40ni_48 data/processedData/aligned_reads/quality_filter/C_D14_Rep_1.PE.star.sorted.fixed.filtered.postdup.RG.passed.bam --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --no-hooks --nolock --timestamp --mode 2  --allowed-rules rSubread  && touch "/gpfs/scratch/cah422/stephensiTranscriptome-november2018-mayvInf/.snakemake/tmp.e40ni_48/17.jobfinished" || (touch "/gpfs/scratch/cah422/stephensiTranscriptome-november2018-mayvInf/.snakemake/tmp.e40ni_48/17.jobfailed"; exit 1)

