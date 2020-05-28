#install.packages("devtools")
#devtools::install_github("markziemann/dee2/getDEE2")
#devtools::install_github("markziemann/Mitch")
#devtools::install_github("hrbrmstr/taucharts")

library("devtools")
library("getDEE2") # install getDEE2_0.0.4 with devtools::install_github("markziemann/dee2/getDEE2")
library("mitch")
library("DESeq2")


##################################################
# Obtain gene expression counts and run DESeq2
##################################################
mdat<-getDee2Metadata("hsapiens")

#SRP128998
samplesheet<-mdat[which(mdat$GEO_series=="GSE109140"),]
samplesheet<-samplesheet[order(samplesheet$SRR_accession),]
samplesheet$HG<-as.factor(c(1,1,1,1,1,1,0,0,0,0,0,0))
samplesheet$VPA<-as.factor(c(0,0,0,1,1,1,0,0,0,1,1,1))

# samplesheets for the 2 contrasts
s1 <- subset(samplesheet,VPA==0)
s2 <- subset(samplesheet,HG==1)


w<-getDEE2("hsapiens",samplesheet$SRR_accession,metadata=mdat)
x<-Tx2Gene(w)
x<-x$Tx2Gene

# save the genetable for later
gt<-w$GeneInfo[,1,drop=FALSE]
gt$accession<-rownames(gt)

# filter out lowly expressed genes
x<-x[which(rowSums(x)/ncol(x)>=(10)),]

# counts for the 2 contrasts
x1<-x[,which(colnames(x) %in% s1$SRR_accession)]
x2<-x[,which(colnames(x) %in% s2$SRR_accession)]

#run DESeq2 for LG vs HG
y <- DESeqDataSetFromMatrix(countData = round(x1), colData = s1, design = ~ HG)
y <- DESeq(y)
LGvHG <- results(y)
LGvHG<-as.data.frame(LGvHG[order(LGvHG$pvalue),])
rownames(LGvHG)<-sapply(strsplit(rownames(LGvHG),"\\."),"[[",1)

# run DESeq2 for HG vs HGVPA
y <- DESeqDataSetFromMatrix(countData = round(x2), colData = s2, design = ~ VPA)
y <- DESeq(y)
HGvHGVPA <- results(y)
HGvHGVPA<-as.data.frame(HGvHGVPA[order(HGvHGVPA$pvalue),])
rownames(HGvHGVPA)<-sapply(strsplit(rownames(HGvHGVPA),"\\."),"[[",1)

# import
d<-list("LGvHG"=LGvHG,"HGvHGVPA"=HGvHGVPA)
d<-mitch_import(d,DEtype="deseq2",geneTable=gt)

##################################################
# get the gene sets and import
##################################################
download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip", 
    destfile="ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip")
genesets<-gmt_import("ReactomePathways.gmt")

##################################################
# run the analysis significance vs effect
##################################################
res<-mitch_calc(d,genesets,resrows=50,priority="effect")
mitch_plots(res,outfile="HGVPA_eff.pdf")
mitch_report(res,"HGVPA_eff.html")

res<-mitch_calc(d,genesets,resrows=50,priority="significance")
mitch_plots(res,outfile="HGVPA_sig.pdf")
mitch_report(res,"HGVPA_sig.html")

