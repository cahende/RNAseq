#!/storage/home/cah422/work/bin/R-3.5.1/bin

# -------------------------------------------------------------------------------------- #
#RNAseq - Get read counts from filtered and mapped reads - Rsubread (Bioconductor Package)
# -------------------------------------------------------------------------------------- #

#Load packages
library(Rsubread)

#Begin writing output to file
sink(snakemake@log[[1]])

#Run feature counts on quality filtered RNAseq reads mapped to reference genome
counts <- featureCounts(snakemake@input[[1]], 
              annot.ext=snakemake@config[["GENOME_ANNOTATION"]],
              genome=snakemake@config[["GENOME"]],
              isGTFAnnotationFile=TRUE,
              isPairedEnd=TRUE)

write.table(x=data.frame(counts$annotation[ , c("GeneID")], counts$counts, stringsAsFactors=FALSE), file=snakemake@output[[1]], quote=FALSE, sep="\t",r ow.names=FALSE, col.names=FALSE)
