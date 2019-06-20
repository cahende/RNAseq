# -------------------------------------------------------------------------------------- #
#RNAseq - Get read counts from filtered and mapped reads - Rsubread (Bioconductor Package)
# -------------------------------------------------------------------------------------- #

#Load packages
library(Rsubread)
library(edgeR)

#Run feature counts on quality filtered RNAseq reads mapped to reference genome
counts <- featureCounts("data/zikv1/zikv1.final.bam", 
                        annot.ext="genomes/AagL5.gtf",
                        genome="genomes/AagL5.fa",
                        isGTFAnnotationFile=TRUE,
                        isPairedEnd=TRUE)
write.table(x=data.frame(counts$annotation[ , c("GeneID")], counts$counts, stringsAsFactors=FALSE), 
            file="data/zikv1/zikv1.readCounts.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

counts <- featureCounts("data/zikv2/zikv2.final.bam", 
                        annot.ext="genomes/AagL5.gtf",
                        genome="genomes/AagL5.fa",
                        isGTFAnnotationFile=TRUE,
                        isPairedEnd=TRUE)
write.table(x=data.frame(counts$annotation[ , c("GeneID")], counts$counts, stringsAsFactors=FALSE), 
            file="data/zikv2/zikv2.readCounts.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

counts <- featureCounts("data/zikv3/zikv3.final.bam", 
                        annot.ext="genomes/AagL5.gtf",
                        genome="genomes/AagL5.fa",
                        isGTFAnnotationFile=TRUE,
                        isPairedEnd=TRUE)
write.table(x=data.frame(counts$annotation[ , c("GeneID")], counts$counts, stringsAsFactors=FALSE), 
            file="data/zikv3/zikv3.readCounts.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

counts <- featureCounts("data/ctrl1/ctrl1.final.bam", 
                        annot.ext="genomes/AagL5.gtf",
                        genome="genomes/AagL5.fa",
                        isGTFAnnotationFile=TRUE,
                        isPairedEnd=TRUE)
write.table(x=data.frame(counts$annotation[ , c("GeneID")], counts$counts, stringsAsFactors=FALSE), 
            file="data/ctrl1/ctrl1.readCounts.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

counts <- featureCounts("data/ctrl2/ctrl2.final.bam", 
                        annot.ext="genomes/AagL5.gtf",
                        genome="genomes/AagL5.fa",
                        isGTFAnnotationFile=TRUE,
                        isPairedEnd=TRUE)
write.table(x=data.frame(counts$annotation[ , c("GeneID")], counts$counts, stringsAsFactors=FALSE), 
            file="data/ctrl2/ctrl2.readCounts.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

counts <- featureCounts("data/ctrl3/ctrl3.final.bam", 
                        annot.ext="genomes/AagL5.gtf",
                        genome="genomes/AagL5.fa",
                        isGTFAnnotationFile=TRUE,
                        isPairedEnd=TRUE)
write.table(x=data.frame(counts$annotation[ , c("GeneID")], counts$counts, stringsAsFactors=FALSE), 
            file="data/ctrl3/ctrl3.readCounts.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

#Use readDGE in edgeR to create a DGEList containing all read counts for all samples in one file
group <- factor(c(1,1,1,2,2,2))
counts <- readDGE(c("data/zikv1/zikv1.readCounts.txt", 
                    "data/zikv2/zikv2.readCounts.txt", 
                    "data/zikv3/zikv3.readCounts.txt",
                    "data/ctrl1/ctrl1.readCounts.txt", 
                    "data/ctrl2/ctrl2.readCounts.txt", 
                    "data/ctrl3/ctrl3.readCounts.txt"), 
                  columns=c(1,2), group=group)

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

#Isolate all genes and P-value in a named (gene ID) vector with numeric values (P-value)
qlfTable <- qlf[["table"]]
qlfNamedVector <- qlfTable$PValue
names(qlfNamedVector) <- row.names(qlfTable)

#Identify interesting genes
for (i in 1:nrow(qlfTable)) {
  if (qlfTable$PValue[i] <= 0.05 & (qlfTable$logFC[i] > 2 | qlfTable$logFC[i] < -2)) {
    qlfTable$interesting[i] <- TRUE
  }
  else {
    qlfTable$interesting[i] <- FALSE
  }
}

#Write table to file
write.table(qlfTable, sep = ",", col.names = NA, row.names = TRUE, 
            file = "diffExp.csv")
