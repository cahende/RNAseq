#!/$HOME/work/bin/R-3.5.1/bin/

# -------------------------------------------------------------------------------------- #
#RNAseq - Get read counts from filtered and mapped reads - Rsubread (Bioconductor Package)
# -------------------------------------------------------------------------------------- #

#install bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")
BiocManager::install("Rsubread", version = "3.8", lib = "/storage/home/cah422/work/bin/R-3.5.1/library/")

#Test with load packages
library(Rsubread)
