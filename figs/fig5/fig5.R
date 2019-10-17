library("Seurat")
library("mitch")
library("gplots")

load("seurat_data.RData")

# the tbl object contains all the DE contrasts for different cell types
str(tbl)

x<-mitch_import(tbl,DEtype="seurat")

# get gene sets from Reactome
download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip","ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip")
genesets<-gmt_import("ReactomePathways.gmt")

# define a scale bar function
image.scale <- function(z, zlim, col = heat.colors(12),
breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
 if(!missing(breaks)){
  if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
 }
 if(missing(breaks) & !missing(zlim)){
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
 }
 if(missing(breaks) & missing(zlim)){
  zlim <- range(z, na.rm=TRUE)
  zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
  zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
 }
 poly <- vector(mode="list", length(col))
 for(i in seq(poly)){
  poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
 }
 xaxt <- ifelse(horiz, "s", "n")
 yaxt <- ifelse(horiz, "n", "s")
 if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
 if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
 if(missing(xlim)) xlim=XLIM
 if(missing(ylim)) ylim=YLIM
 plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
 for(i in seq(poly)){
  if(horiz){
   polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
  }
  if(!horiz){
   polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
  }
 }
}

# generate some plots
pdf("fig5.pdf")

# first look at gene level DE correlation 
xc<-cor(x,method="p",use="pairwise.complete.obs")
warm_palette <- colorRampPalette(c("dark red","red", "orange", "yellow","white"))(n = 25)
heatmap.2(xc,margin=c(20, 20),cexRow=0.8,trace="none",cexCol=0.8,col=warm_palette,scale="none",main="gene level DE Pearson correlation")
xc<-cor(x,method="s",use="pairwise.complete.obs")
heatmap.2(xc,margin=c(20, 20),cexRow=0.8,trace="none",cexCol=0.8,col=warm_palette,scale="none",main="gene level DE Spearman correlation")

max(apply(xc,2,function(x) {  max(x[order(-x)][2:length(x)])      } ) )
min(cx)
#heatmap(xc,scale="none",main="gene level DE Pearson correlation")

pal.1=colorRampPalette(c("red", "yellow","white"), space="rgb")
breaks <- seq(min(xc), max(xc),length.out=25)
image.scale(xc, col=pal.1(length(breaks)-1), breaks=breaks, horiz=TRUE)

#exclude erythrocytes
x<-x[,which(colnames(x)!="eryth")]

# run mitch with significance prioritisation
res<-mitch_calc(x,genesets,resrows=50,priority="significance")

# create a barplot of DP numbers in each cell type
par(mar=c(5,10,4,2))
counts<-apply(res$enrichment_result[,16:27], 2 , function(x) { length(which(p.adjust(x) <0.05) ) } ) 
counts<-counts[order(-counts)]
barplot( counts, horiz=T,las=1,xlab="no gene sets FDR<0.05")

# histogram of gene sets by cell type
par(mai=c(1.02,0.82,0.82,0.42))
sig<-res$enrichment_result[which(res$enrichment_result$p.adjustMANOVA<0.05),4:15]
sig<-sign(t(sig))
neg<-as.vector(apply(sig,2,function(x) { length(which(x==-1)) } ))
neg<-neg[order(neg)]
pos<-as.vector(apply(sig,2,function(x) { length(which(x==1)) } ))
pos<-pos[order(-pos)]
hist(pos,xlab="number of cell clusters that gene set s is positive",ylab="number of gene sets",main="Regulation of FDR<0.05 gene sets")

# heatmap prioritise by significance
res<-mitch_calc(x,genesets,resrows=50,priority="significance")
sig <- res$enrichment_result[1:30,c(1,4:15)]
rownames(sig)<-sig$set
sig$set=NULL
colnames(sig)<-gsub("s.","",colnames(sig))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
heatmap.2(as.matrix(sig),margin=c(10, 22),cexRow=0.65,trace="none",cexCol=0.8,col=my_palette,scale="none",main="top sets by significance")

# 25
sig <- res$enrichment_result[1:25,c(1,4:15)]
rownames(sig)<-sig$set
sig$set=NULL
colnames(sig)<-gsub("s.","",colnames(sig))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
heatmap.2(as.matrix(sig),margin=c(15, 22),cexRow=0.65,trace="none",cexCol=0.8,col=my_palette,scale="none",main="top sets by significance")


# heatmap prioritise by effect
res<-mitch_calc(x,genesets,resrows=50,priority="effect")
sig <- res$enrichment_result[1:30,c(1,4:15)]
rownames(sig)<-sig$set
sig$set=NULL
colnames(sig)<-gsub("s.","",colnames(sig))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
heatmap.2(as.matrix(sig),margin=c(10, 22),cexRow=0.65,trace="none",cexCol=0.8,col=my_palette,scale="none",main="top sets by effect")

# 25
sig <- res$enrichment_result[1:25,c(1,4:15)]
rownames(sig)<-sig$set
sig$set=NULL
colnames(sig)<-gsub("s.","",colnames(sig))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
heatmap.2(as.matrix(sig),margin=c(15, 22),cexRow=0.65,trace="none",cexCol=0.8,col=my_palette,scale="none",main="top sets by effect")


# heatmap prioritise by SD
res<-mitch_calc(x,genesets,resrows=50,priority="SD")
sig <- res$enrichment_result[1:30,c(1,4:15)]
rownames(sig)<-sig$set
sig$set=NULL
colnames(sig)<-gsub("s.","",colnames(sig))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
heatmap.2(as.matrix(sig),margin=c(10, 22),cexRow=0.65,trace="none",cexCol=0.8,col=my_palette,scale="none",main="top sets by SD")

# 25
sig <- res$enrichment_result[1:25,c(1,4:15)]
rownames(sig)<-sig$set
sig$set=NULL
colnames(sig)<-gsub("s.","",colnames(sig))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
heatmap.2(as.matrix(sig),margin=c(15, 22),cexRow=0.65,trace="none",cexCol=0.8,col=my_palette,scale="none",main="top sets by SD")

dev.off()


mitch_plots(res,outfile="scrnaseq.pdf")
mitch_report(res,"scrnaseq.html")

