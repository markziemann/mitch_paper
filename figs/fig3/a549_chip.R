library("plyr")
library("dplyr")  
library("DESeq2")

pdf("MDS_plots.pdf")
#####################################################
#POL2RA
#####################################################
pol2ra<-as.data.frame(c(1,1,1,1,0,0,0,0))
rownames(pol2ra)<-c("ENCFF380ZHO","ENCFF960AWO","ENCFF471JEL","ENCFF593TPG","ENCFF984VJJ","ENCFF697JDL","ENCFF223EHQ","ENCFF565TWZ")
colnames(pol2ra)="dex"

ENCFF380ZHO<-read.table("ENCFF380ZHO.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF960AWO<-read.table("ENCFF960AWO.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF471JEL<-read.table("ENCFF471JEL.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF593TPG<-read.table("ENCFF593TPG.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF984VJJ<-read.table("ENCFF984VJJ.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF697JDL<-read.table("ENCFF697JDL.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF223EHQ<-read.table("ENCFF223EHQ.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF565TWZ<-read.table("ENCFF565TWZ.bam.cnt",header=T)[,c(1,7),drop=F]

x<-join_all(list(ENCFF380ZHO,ENCFF960AWO,ENCFF471JEL,ENCFF593TPG,ENCFF984VJJ,ENCFF697JDL,ENCFF223EHQ,ENCFF565TWZ),by="Geneid")

rownames(x)<-x$Geneid
x$Geneid=NULL
colnames(x)<-c("ENCFF380ZHO","ENCFF960AWO","ENCFF471JEL","ENCFF593TPG","ENCFF984VJJ","ENCFF697JDL","ENCFF223EHQ","ENCFF565TWZ")

y<-x[which(rowSums(x)/ncol(x)>=(10)),]

#deseq2
dds <- DESeqDataSetFromMatrix(countData = y, colData = pol2ra, design = ~ dex)
res <- DESeq(dds)
z<- results(res)
z<-z[order(z$pvalue),]
write.table(z,"pol2ra_res.tsv",sep="\t")

mds<-scale(y)
plot(cmdscale(dist(t(mds))),xlab="Coordinate 1", ylab="Coordinate 2", type = "p",pch=19,main="MDS POL2RA")
text(cmdscale(dist(t(mds)))+1, labels=colnames(mds))

#####################################################
# H3K4me3
#####################################################
h3k4me3<-as.data.frame(c(1,1,1,0,0,0))
rownames(h3k4me3)<-c("ENCFF626NXS","ENCFF978RET","ENCFF945WGW","ENCFF973TUQ","ENCFF428UWO","ENCFF643FMK")
colnames(h3k4me3)="dex"

ENCFF626NXS<-read.table("ENCFF626NXS.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF978RET<-read.table("ENCFF978RET.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF945WGW<-read.table("ENCFF945WGW.bam.cnt",header=T)[,c(1,7),drop=F]

ENCFF973TUQ<-read.table("ENCFF973TUQ.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF428UWO<-read.table("ENCFF428UWO.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF643FMK<-read.table("ENCFF643FMK.bam.cnt",header=T)[,c(1,7),drop=F]

x<-join_all(list(ENCFF626NXS,ENCFF978RET,ENCFF945WGW,ENCFF973TUQ,ENCFF428UWO,ENCFF643FMK),by="Geneid")
rownames(x)<-x$Geneid
x$Geneid=NULL
colnames(x)<-c("ENCFF626NXS","ENCFF978RET","ENCFF945WGW","ENCFF973TUQ","ENCFF428UWO","ENCFF643FMK")

y<-x[which(rowSums(x)/ncol(x)>=(10)),]

dds <- DESeqDataSetFromMatrix(countData = y, colData = h3k4me3, design = ~ dex)
res <- DESeq(dds)
z<- results(res)
z<-z[order(z$pvalue),]
write.table(z,"h3k4me3_res.tsv",sep="\t")

mds<-scale(y)
plot(cmdscale(dist(t(mds))),xlab="Coordinate 1", ylab="Coordinate 2", type = "p",pch=19,main="MDS H3K4me3")
text(cmdscale(dist(t(mds)))+1, labels=colnames(mds))

#####################################################
# NR3C1
#####################################################
nr3c1<-as.data.frame(c(1,1,1,0,0,0,0,0))
#check 2016 and 2017 (ENCFF668EHX,ENCFF496BTD,ENCFF681JBZ) on MDS plot
rownames(nr3c1)<-c("ENCFF807YIG","ENCFF331QXR","ENCFF038DBH","ENCFF668EHX","ENCFF496BTD","ENCFF681JBZ","ENCFF181HLP","ENCFF870WJP")
colnames(nr3c1)="dex"

ENCFF807YIG<-read.table("ENCFF807YIG.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF331QXR<-read.table("ENCFF331QXR.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF038DBH<-read.table("ENCFF038DBH.bam.cnt",header=T)[,c(1,7),drop=F]

ENCFF668EHX<-read.table("ENCFF668EHX.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF496BTD<-read.table("ENCFF496BTD.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF681JBZ<-read.table("ENCFF681JBZ.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF181HLP<-read.table("ENCFF181HLP.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF870WJP<-read.table("ENCFF870WJP.bam.cnt",header=T)[,c(1,7),drop=F]

x<-join_all(list(ENCFF807YIG,ENCFF331QXR,ENCFF038DBH,ENCFF668EHX,ENCFF496BTD,ENCFF681JBZ,ENCFF181HLP,ENCFF870WJP),by="Geneid")
rownames(x)<-x$Geneid
x$Geneid=NULL
colnames(x)<-c("ENCFF807YIG","ENCFF331QXR","ENCFF038DBH","ENCFF668EHX","ENCFF496BTD","ENCFF681JBZ","ENCFF181HLP","ENCFF870WJP")

y<-x[which(rowSums(x)/ncol(x)>=(10)),]

dds <- DESeqDataSetFromMatrix(countData = y, colData = nr3c1, design = ~ dex)
res <- DESeq(dds)
z<- results(res)
z<-z[order(z$pvalue),]
write.table(z,"nr3c1_res.tsv",sep="\t")

mds<-scale(y)
plot(cmdscale(dist(t(mds))),xlab="Coordinate 1", ylab="Coordinate 2", type = "p",pch=19,main="MDS NR3C1")
text(cmdscale(dist(t(mds)))+1, labels=colnames(mds))

#####################################################
# CTCF
#####################################################
ctcf<-as.data.frame(c(1,1,1,0,0,0))
#check 2016 and 2017 (ENCFF668EHX,ENCFF496BTD,ENCFF681JBZ) on MDS plot
rownames(ctcf)<-c("ENCFF356NSD","ENCFF951KQC","ENCFF180EPT","ENCFF774IBY","ENCFF713UMA","ENCFF810IXF")
colnames(ctcf)="dex"

ENCFF356NSD<-read.table("ENCFF356NSD.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF951KQC<-read.table("ENCFF951KQC.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF180EPT<-read.table("ENCFF180EPT.bam.cnt",header=T)[,c(1,7),drop=F]

ENCFF774IBY<-read.table("ENCFF774IBY.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF713UMA<-read.table("ENCFF713UMA.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF810IXF<-read.table("ENCFF810IXF.bam.cnt",header=T)[,c(1,7),drop=F]

x<-join_all(list(ENCFF356NSD,ENCFF951KQC,ENCFF180EPT,ENCFF774IBY,ENCFF713UMA,ENCFF810IXF),by="Geneid")
rownames(x)<-x$Geneid
x$Geneid=NULL
colnames(x)<-c("ENCFF356NSD","ENCFF951KQC","ENCFF180EPT","ENCFF774IBY","ENCFF713UMA","ENCFF810IXF")

y<-x[which(rowSums(x)/ncol(x)>=(10)),]

dds <- DESeqDataSetFromMatrix(countData = y, colData = ctcf, design = ~ dex)
res <- DESeq(dds)
z<- results(res)
z<-z[order(z$pvalue),]
write.table(z,"ctcf_res.tsv",sep="\t")

mds<-scale(y)
plot(cmdscale(dist(t(mds))),xlab="Coordinate 1", ylab="Coordinate 2", type = "p",pch=19,main="MDS CTCF")
text(cmdscale(dist(t(mds)))+1, labels=colnames(mds))


#####################################################
# ATAC
#####################################################
atac<-as.data.frame(c(1,1,1,0,0,0))
#check 2016 and 2017 (ENCFF668EHX,ENCFF496BTD,ENCFF681JBZ) on MDS plot
rownames(atac)<-c("ENCFF020COS","ENCFF758ORC","ENCFF809EKV","ENCFF597SLV","ENCFF248HDS","ENCFF978DQZ")
colnames(atac)="dex"

ENCFF020COS<-read.table("ENCFF020COS.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF758ORC<-read.table("ENCFF758ORC.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF809EKV<-read.table("ENCFF809EKV.bam.cnt",header=T)[,c(1,7),drop=F]

ENCFF597SLV<-read.table("ENCFF597SLV.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF248HDS<-read.table("ENCFF248HDS.bam.cnt",header=T)[,c(1,7),drop=F]
ENCFF978DQZ<-read.table("ENCFF978DQZ.bam.cnt",header=T)[,c(1,7),drop=F]

x<-join_all(list(ENCFF020COS,ENCFF758ORC,ENCFF809EKV,ENCFF597SLV,ENCFF248HDS,ENCFF978DQZ),by="Geneid")
rownames(x)<-x$Geneid
x$Geneid=NULL
colnames(x)<-c("ENCFF020COS","ENCFF758ORC","ENCFF809EKV","ENCFF597SLV","ENCFF248HDS","ENCFF978DQZ")

y<-x[which(rowSums(x)/ncol(x)>=(10)),]

dds <- DESeqDataSetFromMatrix(countData = y, colData = atac, design = ~ dex)
res <- DESeq(dds)
z<- results(res)
z<-z[order(z$pvalue),]
write.table(z,"atac_res.tsv",sep="\t")

mds<-scale(y)
plot(cmdscale(dist(t(mds))),xlab="Coordinate 1", ylab="Coordinate 2", type = "p",pch=19,main="MDS ATAC")
text(cmdscale(dist(t(mds)))+1, labels=colnames(mds))

dev.off()
