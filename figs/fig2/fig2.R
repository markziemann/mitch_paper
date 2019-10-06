library("mitch")
load("fig1.RData")
# profile data is "d"
# genesets are "genesets"
x<-d

#run the analysis
res<-mitch_calc(x,genesets,resrows=25,priority="significance")
mitch_report(res,"unrand.html")

#randomise rownames test
xx<-x
rownames(xx)<-sample(rownames(x))
head(xx)
res<-mitch_calc(xx,genesets,priority="significance")
mitch_report(res,"randres_names.html")

#randomise all data
xx<-x
xx[,1]<-sample(x[,1])
xx[,2]<-sample(x[,2])
head(xx)
res<-mitch_calc(xx,genesets,priority="significance")
mitch_report(res,"randres_xy.html")


#run some permutes
PERMUTES=1000

##################################################
# Randomise gene names
##################################################
numsig1_05=NULL
numsig1_01=NULL
print(paste("Running",PERMUTES,"permutes"))
for (i in 1:PERMUTES) {
set.seed(i)
xx<-x
rownames(xx)<-sample(rownames(x))
head(xx)
res<-mitch_calc(xx,genesets)
n_05<-length(which(res$manova_result$p.adjustMANOVA<0.05))
n_01<-length(which(res$manova_result$p.adjustMANOVA<0.01))
numsig1_05=c(numsig1_05,n_05)
numsig1_01=c(numsig1_01,n_01)
}
numsig1_05<-as.data.frame(numsig1_05)
numsig1_01<-as.data.frame(numsig1_01)

##################################################
# Randomise profile data
##################################################
numsig2_05=NULL
numsig2_01=NULL
print(paste("Running",PERMUTES,"permutes"))
for (i in 1:PERMUTES) {
xx<-x
set.seed(i)
xx[,1]<-sample(x[,1])
set.seed(i+i)
xx[,2]<-sample(x[,2])
res<-mitch_calc(xx,genesets)
n_05<-length(which(res$manova_result$p.adjustMANOVA<0.05))
n_01<-length(which(res$manova_result$p.adjustMANOVA<0.01))
numsig2_05=c(numsig2_05,n_05)
numsig2_01=c(numsig2_01,n_01)
}
numsig2_05<-as.data.frame(numsig2_05)
numsig2_01<-as.data.frame(numsig2_01)

##################################################
# Random gene sets
##################################################
numsig3_05=NULL
numsig3_01=NULL
print(paste("Running",PERMUTES,"permutes"))
for (i in 1:PERMUTES) {
xx<-x
set.seed(i)
randgenesets <- sapply(  unname(unlist(lapply(genesets,length))), function(x) {list(as.character(sample(gt$GeneSymbol,x))) } )
names(randgenesets)<-names(genesets)
res<-mitch_calc(xx,randgenesets)
n_05<-length(which(res$manova_result$p.adjustMANOVA<0.05))
n_01<-length(which(res$manova_result$p.adjustMANOVA<0.01))
numsig3_05=c(numsig3_05,n_05)
numsig3_01=c(numsig3_01,n_01)
}
numsig3_05<-as.data.frame(numsig3_05)
numsig3_01<-as.data.frame(numsig3_01)

save.image("randres.RData")

MAX=max(c(max(numsig1),max(numsig2)))

pdf("randres.pdf")
par(mfrow=c(3,1))
plot(numsig1_05$numsig1_05,main="Gene name randomisation",ylab="No. FDR MANOVA<0.05 sets",xlab="Run",pch=19)
plot(numsig2_05$numsig2_05,main="Profile randomisation",ylab="No. FDR MANOVA<0.05 sets",xlab="Run",pch=19)
plot(numsig3_05$numsig3_05,main="Gene set randomisation",ylab="No. FDR MANOVA<0.05 sets",xlab="Run",pch=19)

par(mfrow=c(3,1))
BP=barplot(table(numsig1_05$numsig1_05),ylab="no. runs", xlab="no. false positives",ylim=c(0,1100))
text(BP, as.vector(table(numsig1_05$numsig1_05)), labels= as.vector(table(numsig1_05$numsig1_05)), pos=3)
BP=barplot(table(numsig2_05$numsig2_05),ylab="no. runs", xlab="no. false positives",ylim=c(0,1100))
text(BP, as.vector(table(numsig2_05$numsig2_05)), labels= as.vector(table(numsig2_05$numsig2_05)), pos=3)
BP=barplot(table(numsig3_05$numsig3_05),ylab="no. runs", xlab="no. false positives",ylim=c(0,1100))
text(BP, as.vector(table(numsig3_05$numsig3_05)), labels= as.vector(table(numsig3_05$numsig3_05)), pos=3)

par(mfrow=c(3,1))
plot(numsig1_01$numsig1_01,main="Gene name randomisation",ylab="No. FDR MANOVA<0.01 sets",xlab="Run",pch=19)
plot(numsig2_01$numsig2_01,main="Profile randomisation",ylab="No. FDR MANOVA<0.01 sets",xlab="Run",pch=19)
plot(numsig3_01$numsig3_01,main="Gene set randomisation",ylab="No. FDR MANOVA<0.01 sets",xlab="Run",pch=19)

par(mfrow=c(3,1))
BP=barplot(table(numsig1_01$numsig1_01),ylab="no. runs", xlab="no. false positives",ylim=c(0,1100))
text(BP, as.vector(table(numsig1_01$numsig1_01)), labels= as.vector(table(numsig1_01$numsig1_01)), pos=3)
BP=barplot(table(numsig2_01$numsig2_01),ylab="no. runs", xlab="no. false positives",ylim=c(0,1100))
text(BP, as.vector(table(numsig2_01$numsig2_01)), labels= as.vector(table(numsig2_01$numsig2_01)), pos=3)
BP=barplot(table(numsig3_01$numsig3_01),ylab="no. runs", xlab="no. false positives",ylim=c(0,1100))
text(BP, as.vector(table(numsig3_01$numsig3_01)), labels= as.vector(table(numsig3_01$numsig3_01)), pos=3)

dev.off()


save.image("randres.RData")
