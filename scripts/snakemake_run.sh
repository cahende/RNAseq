#!/bin/bash

# -------------------------------------------------- #
#Run snakemake! Make sure you are in the same directory as Snakefile to execute
# -------------------------------------------------- #

#run snakemake
conda activate bioinfo
snakemake --unlock
snakemake --ri 

exit;
