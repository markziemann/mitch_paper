#source("../fig1/fig1.R")
#save.image("f1.Rdata")
load("f1.Rdata")
library("mitch")
# d is the imported object
r<-res$ranked_profile

########################################
# functions
########################################

# generate some gene sets
randomGeneSets<-function(a){
library("stringi")
gsets<-sapply( rep(50,1000) , function(x) {list(as.character(sample(rownames(a),x))) } )
names(gsets)<-stri_rand_strings(length(gsets), 15, pattern = "[A-Za-z]")
gsets
}

# ranking function
mitch_rank <- function(x) {

    for (i in seq_len(ncol(x))) {
        LEN = length(x[, i])
        UNIQLEN = length(unique(x[, i]))
        if (UNIQLEN/LEN < 0.4) {
            warning("Warning: >60% of genes have the same score.")
        }
    }

    rank_adj <- function(x) {
        xx <- rank(x,na.last = "keep")
        num_neg = length(which(x < 0))
        num_zero = length(which(x == 0))
        num_adj = num_neg + (num_zero/2)
        adj <- xx - num_adj
        adj
    }
    adj <- apply(x, 2, rank_adj)
    adj
}


########################################
# analysis
########################################
gs<-randomGeneSets(r)

# create 20 de sets
de<-sample( gs , (length(gs) / 50 ) )

# create some s scores for these gene sets
# run the analysis for 1% PDR
SMRY=NULL
SD<-c(0,0.01,0.025,0.05,0.1,0.15,0.2,0.25)
for ( sd in SD)  {
    s<-lapply(de,function(x) { rnorm(5, mean = 0, sd = sd) } )
    r1<-sample(r[,1])
    r2<-sample(r[,1])
    r3<-sample(r[,1])
    r4<-sample(r[,1])
    r5<-sample(r[,1])
    rr_orig<-data.frame(r1,r2,r3,r4,r5)
    rr<-data.frame(r1,r2,r3,r4,r5)
    n=nrow(rr)

    smry1pcFDR<-NULL
    smry1pcFDR<-NULL
    for ( rep in seq(1000) ) {
        for ( gset in seq(de) ) {
            for ( dim in seq(length(s[[1]])) ) {
                ss<-s[[gset]][dim]
                diff=ss*n/2
                rrr<-rownames(rr) %in% de[[gset]] *diff
                rr[,dim]<-rr[,dim]+rrr
            }
        }

        rr<-mitch_rank(rr)
        res<-mitch_calc(rr,gs,priority="significance")
        myres<-subset(res$enrichment_result,p.adjustMANOVA<0.01)
        myres$set %in% names(de)
        tp<-length(which(myres$set %in% names(de)))
        fp<-length(which(!myres$set %in% names(de)))
        fn<-length(which(!names(de) %in% myres$set))
        tn=length(gs) - tp - fp - fn
        smry1<-cbind(tp,tn,fp,fn)
        smry1pcFDR<-rbind(smry1pcFDR,smry1)

        myres<-subset(res$enrichment_result,p.adjustMANOVA<0.05)
        myres$set %in% names(de)
        tp<-length(which(myres$set %in% names(de)))
        fp<-length(which(!myres$set %in% names(de)))
        fn<-length(which(!names(de) %in% myres$set))
        tn=length(gs) - tp - fp - fn
        smry1<-cbind(tp,tn,fp,fn)
        smry5pcFDR<-rbind(smry5pcFDR,smry1)

    }

    smry1pcFDR<-as.data.frame(smry1pcFDR)
    TP<-mean(smry1pcFDR$tp)
    TN<-mean(smry1pcFDR$tn)
    FP<-mean(smry1pcFDR$fp)
    FN<-mean(smry1pcFDR$fn)
    A=(TP+TN)/(TP+FP+FN+TN)
    Prec = TP/(TP+FP)
    Rec = TP/(TP+FN)
    F1 = 2*(Rec * Prec) / (Rec + Prec)
    SMRY1<-cbind(TP,TN,FP,FN,A,Prec,Rec,F1)
    SMRY1pcFDR<-rbind(SMRY1pcFDR,SMRY1)

    smry5pcFDR<-as.data.frame(smry5pcFDR)
    TP<-mean(smry5pcFDR$tp)
    TN<-mean(smry5pcFDR$tn)
    FP<-mean(smry5pcFDR$fp)
    FN<-mean(smry5pcFDR$fn)
    A=(TP+TN)/(TP+FP+FN+TN)
    Prec = TP/(TP+FP)
    Rec = TP/(TP+FN)
    F1 = 2*(Rec * Prec) / (Rec + Prec)
    SMRY1<-cbind(TP,TN,FP,FN,A,Prec,Rec,F1)
    SMRY5pcFDR<-rbind(SMRY5pcFDR,SMRY1)

}

SMRY1pcFDR<-as.data.frame(SMRY1pcFDR)
rownames(SMRY1pcFDR)<-SD
write.table(SMRY1pcFDR,file="md_1pcFDR.tsv",sep="\t",quote=FALSE)

SMRY5pcFDR<-as.data.frame(SMRY5pcFDR)
rownames(SMRY5pcFDR)<-SD
write.table(SMRY5pcFDR,file="md_5pcFDR.tsv",sep="\t",quote=FALSE)


# plot results
# md_1pcFDR.tsv  md_5pcFDR.tsv  md.tsv


md1pc<-read.table("md_1pcFDR.tsv",header=TRUE,row.names=1)
md1pc$SD<-rownames(md1pc)
md1pc<-md1pc[2:nrow(md1pc),]

md5pc<-read.table("md_5pcFDR.tsv",header=TRUE,row.names=1)
md5pc$SD<-rownames(md5pc)
md5pc<-md5pc[2:nrow(md5pc),]



pdf("res.pdf")

plot(md1pc$SD,md1pc$Rec,type="b", ylim=c(0.94,1),xlab="SD",col="blue",
  ylab="recall / precision / F1", main="1% FDR")
points(md1pc$SD,md1pc$Prec,type="b",col="red")
points(md1pc$SD,md1pc$F1,type="b",col="black") 

legend("bottomright", legend=c("recall", "precision", "F1 score"),
       col=c("blue", "red", "black"), lty=1, cex=1)

plot(md5pc$SD,md5pc$Rec,type="b", ylim=c(0.85,1),xlab="SD",col="blue",
  ylab="recall / precision / F1", main="5% FDR")
points(md5pc$SD,md5pc$Prec,type="b",col="red")
points(md5pc$SD,md5pc$F1,type="b",col="black") 

legend("bottomright", legend=c("recall", "precision", "F1 score"),
       col=c("blue", "red", "black"), lty=1, cex=1)

dev.off()

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

