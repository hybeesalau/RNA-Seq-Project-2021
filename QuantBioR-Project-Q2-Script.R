#Ibrahim Odunayo Salaudeen QBio Assignmnet

#####  #####   #####
  #    #    #  #   #
  #    #####   ####
  #    #    #  #   #
#####  #####   #    #

#########################
# Install and load DESeq2
#########################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("apeglm")

library(DESeq2)





#########################
# Read in data (the count matrix)
#########################


cts <- as.matrix(read.table("Ibrahim_Odunayo_Salaudeen.CNT.tsv", header=T, row.names=1, sep="\t", check.names = F))

cts <- cts[rownames(cts) != "7SK", ];  # get rid of some genes which are clear outlines

cts.int <- round(cts)





#########################
#DESeq Sample Condition Design with Knockouts(KO) and Wild-Types(WT)  OR generate a data frame of sample conditions
#########################
coldata <- read.table("SampleCondition.txt", header=T, row.names=1, sep="\t")




#########################
#Prefiltering - generate DESeq object
#########################

dds <- DESeqDataSetFromMatrix(countData = cts.int, colData = coldata, design= ~ condition)
dds
#class: DESeqDataSet 
#dim: 11642 8 
#metadata(1): version
#assays(1): counts
#rownames(11642): 0610009B22Rik 0610009O20Rik ... mt-Rnr1 mt-Rnr2
#rowData names(0):
#  colnames(8): Case1 Case2 ... Ctrl3 Ctrl4
#colData names(1): condition


#########################
# Alternatively, pre-filtering can also be done based on counts
#########################

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]



#########################
# Perform DESeq2
#########################

dds <- DESeq(dds)


############
# get sizeFactor based on normalization
############
dds$sizeFactor



############
# get normalized read count
############
dds.norm <- counts(dds, normalize=T)

resultsNames(dds) # lists the coefficients
res <- results(dds, name="condition_WT_vs_KO", alpha=0.05)
res <- res[is.na(res$padj) == 0,]
resOrdered <- res[order(res$pvalue),]



deg.FC2.padj05 <- resOrdered[abs(resOrdered$log2FoldChange)>1 & resOrdered$padj < 0.05,]
write.table(deg.FC2.padj05, file="DESeq2.DEG.FC2.padj05.tsv", sep="\t", quote=F)

####################
# Plot M-A plot
####################

#TO AVOID::::::: Error in plot.new() : figure margins too large,  I added this line of code ---> par(mar=c(1, 1, 1, 1))
par(mar=c(1, 1, 1, 1))
DESeq2::plotMA(res,ylim=c(-2,2), alpha=0.05, main = "M-A plot")

####################
# Plot volcano plot
####################

plot(res$log2FoldChange, -log10(res$padj), main="Volcano plot", 
     xlab="Effect size: log2(fold-change)", 
     ylab="-log10(adjusted p-value)", pch=20, cex=0.6) 

alpha <- 0.05
abline(h=-log10(alpha), col="brown") 
abline(v=0) 
abline(v=c(-1,1), col="brown") 

################
# draw heatmap by pheatmap package
################
install.packages("pheatmap")
library(pheatmap)

#sampleType.df <- as.data.frame(coldata[,"condition"])
#rownames(sampleType.df) <- colnames(dds)
#colnames(sampleType.df) <- "condition"

degName <- rownames(deg.FC2.padj05)

####
# Raw read counts
####
# We can get raw read counts from input count matrix
pheatmap(cts.int[rownames(cts.int) %in% degName,], annotation_col=coldata)
pheatmap(cts.int[rownames(cts.int) %in% degName,], scale="row", annotation_col=coldata)

# We can also get raw read counts from the dds object using counts() function
pheatmap(counts(dds)[rownames(dds) %in% degName,], annotation_col=coldata)
pheatmap(counts(dds)[rownames(dds) %in% degName,], scale="row", annotation_col=coldata)

####
# Normalized read counts
####
# We have to use counts() function to get the DESeq2 normalized read counts from dds object

dds.norm <- counts(dds, normalize=T)

pheatmap(dds.norm[rownames(dds.norm) %in% degName,], annotation_col=coldata)

pheatmap(dds.norm[rownames(dds.norm) %in% degName,], scale="row", annotation_col=coldata)

# We then can transform the normalized counts to log2-tranSformed value for heatmap

dds.norm.log <- log(dds.norm+1)/log(2)
pheatmap(dds.norm.log[rownames(dds.norm.log) %in% degName,],  annotation_col=coldata)
pheatmap(dds.norm.log[rownames(dds.norm.log) %in% degName,], scale="row", annotation_col=coldata)

