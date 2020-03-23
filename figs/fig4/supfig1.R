# suplementary figure to determine whether the dip in mitch precision
# at 0.05 is important or not.

myres<-res$enrichment_result
all<-apply(myres[,4:9],1,function(x) {mean(abs(x)) })

myres<-subset(res$enrichment_result,p.adjustMANOVA<0.05 )
sig<-apply(myres[,4:9],1,function(x) {mean(abs(x)) })

myres<-subset(res$enrichment_result,p.adjustMANOVA>0.05 )
nsig<-apply(myres[,4:9],1,function(x) {mean(abs(x)) })

sig_gt=length(which(sig>0.08))
sig_lt=length(which(sig<0.08))

nsig_gt=length(which(nsig>0.08))
nsig_lt=length(which(nsig<0.08))

slices <- c(sig_gt, sig_lt, nsig_gt, nsig_lt)
lbls <- c("FDR<0.05, mean abs s>0.08", 
"FDR<0.05, mean abs s<0.08", 
"FDR>0.05, mean abs s>0.08",
"FDR>0.05, mean abs s<0.08")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels

pdf("hist.pdf")

hist( all ,xlab="mean abs s",main="all sets")
hist( sig,xlab="mean abs s",main="FDR<0.05 sets")
hist( nsig ,xlab="mean abs s", main="FDR>0.05 sets")
boxplot(all,sig,nsig,ylim=c(0,0.4),names=c("all","FDR<0.05","FDR>0.05"),
ylab="mean absolute s",main="mean abs s scores in case study 2")
grid()
pie(slices, labels = lbls, main="mean abs s values in case study 2",
radius=0.6,cex=0.6)
dev.off()


