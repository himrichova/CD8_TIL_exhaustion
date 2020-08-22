rm(list=ls())


###########################################
# RNA - MDR/MG High vs Low, mm10
###########################################

require(DESeq2)


########################
# set working directory

setwd("data_analysis/star_gene_counts/DESeq2_analysis_7samples/")

# Coldata file where rows of coldata correspond to columns of countData (count_matrix)
coldata=as.data.frame(read.table("coldata.txt",sep="\t", header=TRUE, row.names=1, dec="."))


# Count matrix
infile=as.data.frame(read.table("../count_matrix.txt",sep="\t", header=TRUE, row.names=1, dec=".", fill = TRUE))


COUNTS <- infile
head(COUNTS)
dim(COUNTS)

samples <- colnames(COUNTS)
genes <- rownames(COUNTS)
COUNTS <- apply(as.matrix(COUNTS),2,as.numeric)
rownames(COUNTS) <- genes

Data <- newSeqExpressionSet(counts=COUNTS)



## Filtering absent genes
Filter=0
detectionLimit <- Filter
filter <- rowSums(cpm(exprs(Data))> detectionLimit) >= 0 # cpm = function from edgeR and means 'counts-per-million'
nGenes <- length(filter)
nAbove <- sum(filter)
PpData <- newSeqExpressionSet(counts=exprs(Data)[filter,],featureData=fData(Data)[filter,],phenoData=pData(Data))
head(exprs(PpData))
dim(PpData)


data<-exprs(PpData)

infile<-data
dim(infile)


dds <- DESeqDataSetFromMatrix(countData = infile, colData = coldata, design = ~ condition)


# Run the DEseq
dds <- DESeq(dds)


# specify the conditions to compare:
tissue1="Tumor_high"
tissue2="Tumor_low"

print(paste(tissue1,tissue2,sep="_vs_"))

res <- DESeq2::results(dds, independentFiltering=T, cooksCutoff=T,contrast=c("condition",tissue1,tissue2))

A<-as.data.frame(res)
  
write.table(A, file=paste("DESeq2_RNAseq.",paste(tissue1,tissue2,sep="_vs_"),".Filt_",Filter,".txt",sep=""),quote=FALSE,sep="\t",row.names=T,col.names=T)



mart_annot=read.table("/Volumes/data/groups/lab_winter/himrichova/resources/genomes/mm10/mart_export_mouse_human.txt",sep="\t", header=TRUE, dec=".", fill = TRUE)

A_mart <- merge(A, mart_annot, by.x="row.names", by.y="Mouse_gene_stable_ID")


write.table(unique(A_mart), file=paste("DESeq2_RNAseq.",paste(tissue1,tissue2,sep="_vs_"),".Filt_",Filter,".mart_annot.txt",sep=""),quote=FALSE,sep="\t",row.names=F,col.names=T)
  
  
  

 