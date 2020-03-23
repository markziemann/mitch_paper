library(plyr)
library("DESeq2")

# Curate sample sheet
s<-c(0,0,0,0,1,1,1,1)
rn<-c("ENCFF082ICE","ENCFF467SQA","ENCFF620TAH","ENCFF778BJF","ENCFF744NRK","ENCFF054BTH","ENCFF130DRZ","ENCFF165EET")
samplesheet<-as.data.frame(s)
rownames(samplesheet)<-rn
colnames(samplesheet)="dex"

#download counts
download.file("https://www.encodeproject.org/files/ENCFF082ICE/@@download/ENCFF082ICE.tsv","ENCFF082ICE.tsv")
download.file("https://www.encodeproject.org/files/ENCFF467SQA/@@download/ENCFF467SQA.tsv","ENCFF467SQA.tsv")
download.file("https://www.encodeproject.org/files/ENCFF620TAH/@@download/ENCFF620TAH.tsv","ENCFF620TAH.tsv")
download.file("https://www.encodeproject.org/files/ENCFF778BJF/@@download/ENCFF778BJF.tsv","ENCFF778BJF.tsv")
download.file("https://www.encodeproject.org/files/ENCFF744NRK/@@download/ENCFF744NRK.tsv","ENCFF744NRK.tsv")
download.file("https://www.encodeproject.org/files/ENCFF054BTH/@@download/ENCFF054BTH.tsv","ENCFF054BTH.tsv")
download.file("https://www.encodeproject.org/files/ENCFF130DRZ/@@download/ENCFF130DRZ.tsv","ENCFF130DRZ.tsv")
download.file("https://www.encodeproject.org/files/ENCFF165EET/@@download/ENCFF165EET.tsv","ENCFF165EET.tsv")

#ctrl
ENCFF082ICE<-read.table("ENCFF082ICE.tsv",header=T)[,c(1,7),drop=F]
ENCFF467SQA<-read.table("ENCFF467SQA.tsv",header=T)[,c(1,7),drop=F]
ENCFF620TAH<-read.table("ENCFF620TAH.tsv",header=T)[,c(1,7),drop=F]
ENCFF778BJF<-read.table("ENCFF778BJF.tsv",header=T)[,c(1,7),drop=F]

#trt
ENCFF744NRK<-read.table("ENCFF744NRK.tsv",header=T)[,c(1,7),drop=F]
ENCFF054BTH<-read.table("ENCFF054BTH.tsv",header=T)[,c(1,7),drop=F]
ENCFF130DRZ<-read.table("ENCFF130DRZ.tsv",header=T)[,c(1,7),drop=F]
ENCFF165EET<-read.table("ENCFF165EET.tsv",header=T)[,c(1,7),drop=F]

x<-join_all(list(ENCFF082ICE,ENCFF467SQA,ENCFF620TAH,ENCFF778BJF,ENCFF744NRK,ENCFF054BTH,ENCFF130DRZ,ENCFF165EET),by="Geneid")
colnames(x)<-c("Geneid","ENCFF082ICE","ENCFF467SQA","ENCFF620TAH","ENCFF778BJF","ENCFF744NRK","ENCFF054BTH","ENCFF130DRZ","ENCFF165EET")
rownames(x)<-x$Geneid
x$Geneid=NULL

#10 read filter
y<-x[which(rowSums(x)/ncol(x)>=(10)),]

#deseq2
dds <- DESeqDataSetFromMatrix(countData = y, colData = samplesheet, design = ~ dex)
res <- DESeq(dds)
z<- results(res)
z<-z[order(z$pvalue),]
write.table(z,"rna_res.tsv",sep="\t")

