rm(list=ls())

#####################################################
# ATAC-seq 
#
# DESeq2 anlaysis to identify differentially 
# accessible peaks between MDR/MG High vs MDR/MG Low
# (in vivo or in vitro)
#####################################################

library("DESeq2")

######################

reg_data="ATAC_inVivo_"
#reg_data="ATAC_inVitro_"

# Set working directory including count matrix generated on consensus peak list (in vivo or in vitro)
setwd("data_analysis/coverage_merged_macs2_peaks/DESeq2_High_Low_samples")

list.files()


# count matrix
infile=as.data.frame(read.table("read_counts_matrix.reg_IDs.txt",sep="\t", header=TRUE, row.names=1, dec=".", fill = TRUE))

# coldata file where rows of coldata correspond to columns of count matrix
coldata=as.data.frame(read.table("coldata.txt",sep="\t", header=TRUE, row.names=1, dec="."))


COUNTS <- infile
head(COUNTS)

samples <- colnames(COUNTS)
genes <- rownames(COUNTS)
COUNTS <- apply(as.matrix(COUNTS),2,as.numeric)
rownames(COUNTS) <- genes
head(COUNTS)

Data <- newSeqExpressionSet(counts=COUNTS)

## Filtering absent genes
detectionLimit<-1
filter <- rowSums(cpm(exprs(Data))> detectionLimit) > 0 # cpm = function from edgeR and means 'counts-per-million'
nGenes <- length(filter)
nAbove <- sum(filter)
PpData <- newSeqExpressionSet(counts=exprs(Data)[filter,],featureData=fData(Data)[filter,],phenoData=pData(Data))
head(exprs(PpData))
dim(PpData)


data<-exprs(PpData)

infile<-data
dim(infile)



dds <- DESeqDataSetFromMatrix(countData = infile, colData = coldata, design = ~ tissue)

# Run the DEseq
dds <- DESeq(dds)



# Specify the conditions to be compared:
tissue1="Tumor_high"
tissue2="Tumor_low"
print(paste(tissue1,tissue2,sep="_vs_"))

res <- DESeq2::results(dds, independentFiltering=T, cooksCutoff=T,contrast=c("tissue",tissue1,tissue2))

write.table(res, file=paste(reg_data,"DESeq2.",paste(tissue1,tissue2,sep="_vs_"),".txt",sep=""),quote=FALSE,sep="\t",row.names=T,col.names=T)

A<-as.data.frame(res)
  
  

######################################################################
# Normalize matrix 
# - using  Variance stabilizing Transformation 
######################################################################

vsd <- varianceStabilizingTransformation(dds, blind=TRUE,fitType="parametric")

vstMat<-assay(vsd)
dim(vstMat)

write.table(vstMat, file=paste(reg_data,"_DESeq2normalized-vst.txt",sep="_"),quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)


# Median centering genes normalization
med.vst<-apply(vstMat,1,median)

vst_medianCenteredGenes<-sweep(vstMat,1,med.vst)
head(vst_medianCenteredGenes)

write.table(file=paste(reg_data,"_DESeq2normalized-vst.MedianCenteredGenes.txt",sep=""),vst_medianCenteredGenes,quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE) 




