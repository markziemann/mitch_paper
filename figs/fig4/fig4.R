library("getDEE2")
library("DESeq2")
library("limma")
library("edgeR")
library("ABSSeq")
library("plyr")
library("msigdbr")
library("mitch")
library("gplots")
library("PerformanceAnalytics")
library("reshape2")
library("UpSetR")


set_plot_dimensions <- function(width_choice, height_choice) {
        options(repr.plot.width=width_choice, repr.plot.height=height_choice)
        }

mdat<-getDee2Metadata("hsapiens")
m1<-mdat[which(mdat$SRP_accession %in% "SRP096177"),]
SRRlist<-as.vector(m1$SRR_accession)
SRRlist<-SRRlist[order(SRRlist)]
x<-getDEE2("hsapiens",SRRlist)
x<-Tx2Gene(x)
tx<-x$Tx2Gene

#ctrl=c57 trt=DBA
s<-c(1,1,1,0,0,0)
s<-as.data.frame(s)
rownames(s)<-colnames(tx)

# filter counts
tx1<-tx[which(rowSums(tx)/ncol(tx)>=(10)),]

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(tx1), colData = s, design = ~ s)
dres <- DESeq(dds)
res_deseq2<- DESeq2::results(dres)
res_deseq2$pvalue[is.na(res_deseq2$pvalue)] <- 1
res_deseq2$log2FoldChange[is.na(res_deseq2$log2FoldChange)] <- 0
rnk_deseq2<-as.data.frame( sign(res_deseq2$log2FoldChange) * -log10(res_deseq2$pvalue))
rownames(rnk_deseq2)<-rownames(res_deseq2)
colnames(rnk_deseq2)<-"DESeq2"
rnk_deseq2$geneID<-rownames(rnk_deseq2)

# edgeR GLMRT
design<-model.matrix(~s$s)
rownames(design)<-s$s
y<-DGEList(counts=tx1)
y<-calcNormFactors(y)
y <- estimateDisp(y, design,robust=TRUE)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
res_edger_glmrt<-as.data.frame(topTags(lrt,n=Inf))
rnk_edger_glmrt<-as.data.frame( sign(res_edger_glmrt$logFC) * -log10(res_edger_glmrt$PValue))
rownames(rnk_edger_glmrt)<-rownames(res_edger_glmrt)
colnames(rnk_edger_glmrt)<-"edgeR_glmLRT"
rnk_edger_glmrt$geneID<-rownames(rnk_edger_glmrt)

# edgeR QL
design<-model.matrix(~s$s)
rownames(design)<-s$s
y<-DGEList(counts=tx1)
y<-calcNormFactors(y)
y <- estimateDisp(y, design,robust=TRUE)
fit<-glmQLFit(y, design)
lrt<-glmQLFTest(fit)
res_edger_ql<-as.data.frame(topTags(lrt,n=Inf))
rnk_edger_ql<-as.data.frame( sign(res_edger_ql$logFC) * -log10(res_edger_ql$PValue))
rownames(rnk_edger_ql)<-rownames(res_edger_ql)
colnames(rnk_edger_ql)<-"edgeR_QL"
rnk_edger_ql$geneID<-rownames(rnk_edger_ql)

# Voom-Limma
z<-DGEList(counts=tx1)
z <- calcNormFactors(z)
v <- voom(z,design,plot=F)
fit <- lmFit(v, design)
fit.de <- eBayes(fit, robust=TRUE)
res_voom_limma<-topTable(fit.de,n=Inf)
rnk_voom_limma<-as.data.frame( sign(res_voom_limma$logFC) * -log10(res_voom_limma$P.Value))
rownames(rnk_voom_limma)<-rownames(res_voom_limma)
colnames(rnk_voom_limma)<-"voom_limma"
rnk_voom_limma$geneID<-rownames(rnk_voom_limma)

# ABSseq
obj<-ABSDataSet(tx1, factor(s$s))  #default normalisation is qtotal
obj<-ABSSeq(obj)
dge<- as.data.frame(cbind(obj$Amean,obj$Bmean,-obj$foldChange,obj$pvalue,obj$adj.pvalue))
colnames(dge)=c("Amean","Bmean","logFC","PValue","FDR")
res_absseq<-dge[order(dge$PValue),]
rnk_absseq<-as.data.frame( sign(res_absseq$logFC) * -log10(res_absseq$PValue))
rownames(rnk_absseq)<-rownames(res_absseq)
colnames(rnk_absseq)<-"ABSSeq"
rnk_absseq$geneID<-rownames(rnk_absseq)

# genetable
gt<-unique(x$TxInfo[,1:2])
colnames(gt)<-c("geneID","symbol")
gt<-gt[which(gt$geneID %in% rownames(tx1)),]


# join the tables together
xx<-join_all(list("symbol"=gt,"DESeq2"=rnk_deseq2,
  "edgeR_glmLRT"=rnk_edger_glmrt,
  "edgeR_QL"=rnk_edger_ql,"voom-limma"=rnk_voom_limma,
  "ABSSeq"=rnk_absseq),by="geneID")


# aggregate
xx$geneID=NULL
xx<-aggregate(. ~ symbol,xx,sum)
rownames(xx)<-xx$symbol
xx$symbol=NULL


# gene sets from my website but this needs to be updated
#genesets<-gmt_import("reactome.v5.2.symbols_mouse.gmt")
download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip","ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip")
genesets<-gmt_import("ReactomePathways.gmt")

# run mitch analysis
res<-mitch_calc(xx,genesets,resrows=20,bootstraps=1000,priority="effect")

mitch_plots(res,outfile="fig4_plots.pdf")
mitch_report(res,"fig4_report.html")


pdf("fig4.pdf")

sets_deseq2<-res$manova_res[which(p.adjust(res$manova_res$p.DESeq2)<0.05),]$set
sets_edger_glmrt<-res$manova_res[which(p.adjust(res$manova_res$p.edgeR_GLMRT)<0.05),]$set
sets_edger_ql<-res$manova_res[which(p.adjust(res$manova_res$p.edgeR_QL)<0.05),]$set
sets_absseq<-res$manova_res[which(p.adjust(res$manova_res$p.ABSSeq)<0.05),]$set
sets_voom_limma<-res$manova_res[which(p.adjust(res$manova_res$p.voom_limma)<0.05),]$set


# barplt of number of DE sets
my_lengths<-unlist(lapply(list("ABSSeq"=sets_absseq,"voom-limma"=sets_voom_limma,
  "edgeR glmLRT"=sets_edger_glmrt,"edgeR QL"=sets_edger_ql,"DESeq2"=sets_deseq2),length))

par(mar=c(5,10,4,2))
bp<-barplot(my_lengths,xlab="no. sets FDR<0.05",las=1,horiz=T,xlim=c(0,200))
text(x=my_lengths+10,y=bp,label=my_lengths)

# revert default
par(mai=c(1.02,0.82,0.82,0.42))

svals<-res$manova_res[,4:8]

#MDS plot
plot(cmdscale(dist(t(svals))), xlab="Coordinate 1", ylab="Coordinate 2", type = "p", xlim=c(-1,2))
text(cmdscale(dist(t(svals)))+0.01, labels=gsub("s.","",colnames(svals)))

# correlation  plot
chart.Correlation(svals, histogram=F,method = c("pearson"),main="Pearson correlation of s scores")

# correlation heatmaps
heatmap.2(cor(svals,method="s"),scale="none",margin=c(15, 15),trace="none",
  main="Spearman correlation",cellnote=signif(cor(svals,method="s"),3))

heatmap.2(cor(svals,method="p"),scale="none",margin=c(15, 15),trace="none",
  main="Pearson correlation",cellnote=signif(cor(svals,method="s"),3))


# identify high variance gene set
res$manova_res$sd<-apply(res$manova_res[,4:8],1,sd)

res$manova_res$label<-paste(res$manova_res$set,": set size =",res$manova_res$setSize,", SD =",signif(res$manova_res$sd,2))

res_subset<-head( res$manova_res[order(-res$manova_res$sd),] ,20)
rownames(res_subset)<-res_subset$label

# heatmap of gene sets with high variance
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 25)
heatmap.2(as.matrix(res_subset[1:10,4:8]),scale="row",margin=c(15, 15),cexRow=0.8,trace="none",cexCol=0.8,col=my_palette)


res_list<-list("ABSSeq"=sets_absseq,"voom-limma"=sets_voom_limma,"edgeR glmLRT"=sets_edger_glmrt,"edgeR QL"=sets_edger_ql,"DESeq2"=sets_deseq2)
res_list<-lapply(res_list, function(x){
 data.frame(
  list=x,
  val=1)
})
res_df<-ldply(res_list)

colnames(res_df)<-c("list","set","val")
res_df2<-acast(res_df,set~list,value.var="val",fill=0)
res_df2<-as.data.frame(res_df2)
res_df2$name<-rownames(res_df2)
res_df3<-res_df2[,c(ncol(res_df2),1:(ncol(res_df2)-1))]
upset(res_df3,order.by = "freq",text.scale=2)

subset_genesets<-genesets[which(names(genesets) %in% res_subset$set )]
res2<-mitch_calc(xx,subset_genesets,resrows=20,bootstraps=1000,priority="effect")

heatmap.2(res2$detailed_sets$`Peptide chain elongation`,scale="row",col=my_palette,margin=c(12, 12),trace="none",main='Peptide chain elongation')
par(mar=c(5,10,4,2))
barplot(-log10(as.numeric(as.vector(res2$manova_result[2,9:13]))),horiz=T,xlab="-log10(p-value)",las=1,names.arg=names(res2$manova_result[2,9:13]),main=res2$manova_result[2,1])
par(mai=c(1.02,0.82,0.82,0.42))
dev.off()

pdf("fig4b.pdf",width=10,height=7)
heatmap.2(as.matrix(res2$manova_result[1:10,4:8]),scale="row",margin=c(15, 30),cexRow=0.8,trace="none",cexCol=0.8)
dev.off()


mitch_plots(res,outfile="fig4_subset_plots.pdf")
mitch_report(res,"fig4_subset_report.html")

