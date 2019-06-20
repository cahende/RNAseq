#!/bin/bash

# -------------------------------------------------- #
# format reference genome
# make sure you are in the same directory as the reference genome and gtf file and that only one of each are present
# -------------------------------------------------- #

#load tools
conda activate bioinfo

#generate star index
STAR --runMode genomeGenerate --genomeDir . --genomeFastaFiles *.fa --sjdbGTFfile *.gtf --sjdbOverhang 74

exit;
