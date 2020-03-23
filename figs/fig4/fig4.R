library("mitch")

# if there is an error here, then run the a549_chip.sh shell script
gt<-unique(read.table("gencode.v29.annotation.gtf.tss.bed",stringsAsFactors=F)[,c(6,7)])
gt$V5<-sapply(strsplit(gt$V6,"\\."),"[[",1)
gt$V6=NULL
gt<-unique(gt)
rownames(gt)<-gt$V5

# get Reactome gene sets
download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip","ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip")
genesets<-gmt_import("ReactomePathways.gmt")

# get cistrome gene sets
#download.file("https://raw.githubusercontent.com/markziemann/cistro_gmt/master/HumanTfPeaks.gmt","HumanTfPeaks.gmt")
#genesets<-gmt_import("HumanTfPeaks.gmt")


atac<-read.table("atac_res.tsv")
rownames(atac)<-sapply(strsplit(rownames(atac),"\\."),"[[",1)
ctcf<-read.table("ctcf_res.tsv")
rownames(ctcf)<-sapply(strsplit(rownames(ctcf),"\\."),"[[",1)
h3k4me3<-read.table("h3k4me3_res.tsv")
rownames(h3k4me3)<-sapply(strsplit(rownames(h3k4me3),"\\."),"[[",1)
nr3c1<-read.table("nr3c1_res.tsv")
rownames(nr3c1)<-sapply(strsplit(rownames(nr3c1),"\\."),"[[",1)
pol2ra<-read.table("pol2ra_res.tsv")
rownames(pol2ra)<-sapply(strsplit(rownames(pol2ra),"\\."),"[[",1)
rna<-read.table("rna_res.tsv")
rownames(rna)<-sapply(strsplit(rownames(rna),"\\."),"[[",1)

x<-list("CTCF"=ctcf,"H3K4me3"=h3k4me3,"NR3C1"=nr3c1,"POL2RA"=pol2ra,"ATAC"=atac,"RNA"=rna)


#x<-list("ctcf"=ctcf,"h3k4me3"=h3k4me3,"nr3c1"=nr3c1,"pol2ra"=pol2ra,"rna"=rna)

#x<-list("pol2ra"=pol2ra,"rna"=rna)

y<-mitch_import(x, DEtype="deseq2" , geneTable = gt)


res<-mitch_calc(y,genesets,resrows=50,bootstraps=1000,priority="confidence")

mitch_plots(res,outfile="a549_res_conf.pdf")

res<-mitch_calc(y,genesets,resrows=50,bootstraps=1000,priority="significance")

mitch_plots(res,outfile="a549_res_sig.pdf")



