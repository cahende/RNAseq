#!/storage/home/cah422/work/bin/R-3.5.1/bin

# -------------------------------------------------------------------------------------- #
#RNAseq - Combine read counts into single file and create a DGEList - edgeR (Bioconductor Package)
# -------------------------------------------------------------------------------------- #

#Load packages
library(edgeR)

#Begin writing output to file
sink(snakemake@log[[1]])

#Use readDGE in edgeR to create a DGEList containing all read counts for all samples in one file
group <- factor(c(snakemake@config[["GROUP_ID"]]))
counts <- readDGE(snakemake@input, columns=c(1,2), group=group)

#Filter out lowly expressed genes
keep <- rowSums(cpm(counts)>1) >= 2
counts <- counts[keep, , keep.lib.sizes=FALSE]

#Correct for RNA composition
counts <- calcNormFactors(counts)

#Estimate common dispersion from mean
design <- model.matrix(~group)
counts <- estimateDisp(counts, design)

#GLM to look for differentially expressed genes
fit <- glmQLFit(counts, design)
qlf <- glmQLFTest(fit, coef=2)
tagsPValue <- topTags(qlf, n=10000, sort.by="PValue")
tagsFoldChange <- topTags(qlf, n=10000, sort.by"logFC")

#Export top hits
write.table(x=tagsPValue, file=snakemake@output[[1]], quote=FALSE,sep="\t", row.names=TRUE, col.names=TRUE)
write.table(x=tagsFoldChange, file=snakemake@output[[2]], quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

#End writing output to file
sink()
