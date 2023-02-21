# This is based on RNA seq analysis in R
#https://bioinformatics-core-shared-training.github.io/RNAseq-R/rna-seq-preprocessing.nb.html

# clear workspace
rm(list = ls())

#libraries 
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)

path <- "data/combined/allSamples/"
# load data counts
seqdata <- read.csv(paste(path,"summarized.csv", sep=""))

# remove first column
countdata <- seqdata[,-1]

# set rownames
rownames(countdata) <- seqdata[,1]

# discard genes with lower reads than 10 in almost all samples
keep <- rowSums(countdata>10)>=2
counts.keep <- countdata[keep,]

# discard samples with total count lower than 1 million
keepSamples <- colSums(countdata)>1000000
counts.keep <- counts.keep[, keepSamples]


# generate DGEList, which is an edgeR-specific object
dgeObj <- DGEList(counts.keep)

# calculate TMM
y <- calcNormFactors(dgeObj)
# Get log2 counts per million
logcounts <- cpm(y,log=TRUE)
write.csv(logcounts, paste(path,"logcounts.csv", sep=""))


# perform normalized MDS plot
mds <- plotMDS(y)
write.csv(mds$cmdscale.out, paste(path,"mds.csv", sep=""))
distance_mat <- mds$distance.matrix
write.csv(distance_mat, paste(path,"distance.csv", sep=""))

