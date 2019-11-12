library("Seurat")
library("mitch")

load("seurat_data.RData")

#unique(immune.combined@meta.data$celltype.stim)

cd14.mono<-FindMarkers(immune.combined, ident.1 = "CD14 Mono_STIM", ident.2 = "CD14 Mono_CTRL", print.bar = FALSE)
pdc<-FindMarkers(immune.combined, ident.1 = "pDC_STIM", ident.2 = "pDC_CTRL", print.bar = FALSE)
cd4.mem <- FindMarkers(immune.combined, ident.1 = "CD4 Memory T_STIM", ident.2 = "CD4 Memory T_CTRL", print.bar = FALSE)
t.act <- FindMarkers(immune.combined, ident.1 = "T activated_STIM", ident.2 = "T activated_CTRL", print.bar = FALSE)
cd4.nat <- FindMarkers(immune.combined, ident.1 = "CD4 Naive T_STIM", ident.2 = "CD4 Naive T_CTRL", print.bar = FALSE)
cd8.t <- FindMarkers(immune.combined, ident.1 = "CD8 T_STIM", ident.2 = "CD8 T_CTRL", print.bar = FALSE)
mk <- FindMarkers(immune.combined, ident.1 = "Mk_STIM", ident.2 = "Mk_CTRL", print.bar = FALSE)
b.act <- FindMarkers(immune.combined, ident.1 = "B activated_STIM", ident.2 = "B activated_CTRL", print.bar = FALSE)
b <- FindMarkers(immune.combined, ident.1 = "B_STIM", ident.2 = "B_CTRL", print.bar = FALSE)
dc <- FindMarkers(immune.combined, ident.1 = "DC_STIM", ident.2 = "DC_CTRL", print.bar = FALSE)
cd16.mono <- FindMarkers(immune.combined, ident.1 = "CD16 Mono_STIM", ident.2 = "CD16 Mono_CTRL", print.bar = FALSE)
nk <- FindMarkers(immune.combined, ident.1 = "NK_STIM", ident.2 = "NK_CTRL", print.bar = FALSE)
eryth <- FindMarkers(immune.combined, ident.1 = "Eryth_STIM", ident.2 = "Eryth_CTRL", print.bar = FALSE)

cell_type_de<-list("cd14.mono"=cd14.mono,"pdc"=pdc,"cd4.mem"=cd4.mem,"t.act"=t.act,
  "cd4.nat"=cd4.nat,"cd8.t"=cd8.t,"mk.b"=mk,"b.act"=b.act,"b"=b,"dc"=dc,
  "cd16.mono"=cd16.mono,"nk"=nk,"eryth"=eryth)

x<-mitch_import(cell_type_de,DEtype="seurat")

# get gene sets from Reactome
download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip","ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip")
genesets<-gmt_import("ReactomePathways.gmt")

# run mitch
res<-mitch_calc(x,genesets,resrows=200,bootstraps=1000,priority="significance")

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
xc<-cor(x,use="pairwise.complete.obs")
heatmap(xc,scale="none")
pal.1=colorRampPalette(c("red", "yellow","white"), space="rgb")
breaks <- seq(min(xc), max(xc),length.out=25)
image.scale(xc, col=pal.1(length(breaks)-1), breaks=breaks, horiz=TRUE)

par(mar=c(5,10,4,2))
counts<-apply(res$manova_result[,17:29], 2 , function(x) { length(which(p.adjust(x) <0.05) ) } ) 
counts<-counts[order(-counts)]

barplot( counts, horiz=T,las=1,xlab="no gene sets FDR<0.05")
par(mai=c(1.02,0.82,0.82,0.42))

sig<-res$manova_result[which(res$manova_result$p.adjustMANOVA<0.05),4:16]
sig<-sign(t(sig))
neg<-as.vector(apply(sig,2,function(x) { length(which(x==-1)) } ))
neg<-neg[order(neg)]
pos<-as.vector(apply(sig,2,function(x) { length(which(x==1)) } ))
pos<-pos[order(-pos)]
hist(pos,xlab="number of cell clusters that gene set s is positive",ylab="number of gene sets",main="Regulation of FDR<0.05 gene sets")

sig<-res$manova_result[which(res$manova_result$p.adjustMANOVA<0.05),c(1,4:16)]
rownames(sig)<-sig$set
sig$set=NULL
colnames(sig)<-gsub("s.","",colnames(sig))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
heatmap.2(as.matrix(sig),margin=c(5, 22),cexRow=0.65,trace="none",cexCol=0.8,col=my_palette,scale="none")


### now prioritise by effect size
res<-mitch_calc(x,genesets,resrows=200,bootstraps=1000,priority="significance")

sig<-res$manova_result[which(res$manova_result$p.adjustMANOVA<0.05),c(1,4:16)]
rownames(sig)<-sig$set
sig$set=NULL
colnames(sig)<-gsub("s.","",colnames(sig))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
heatmap.2(as.matrix(sig),margin=c(5, 22),cexRow=0.65,trace="none",cexCol=0.8,col=my_palette,scale="none")






dev.off()


mitch_plots(res,outfile="scrnaseq.pdf")
mitch_report(res,"scrnaseq.html")

