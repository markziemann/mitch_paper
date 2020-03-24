
library("tidyverse")
library("parallel")
library("topconfects")
library("edgeR")
library("DESeq2")
library("limma")
library("ABSSeq")
library("stringi")
library("mdgsa")


source("simpw2d_func.R")

# obtain count data
a<-countData()

# generate some gene sets
gsets<-randomGeneSets(a)

test_simpw2d<-function(){
 # create some random data with two contraste
 N_REPS=3 ; SUM_COUNT=40000000 ; VARIANCE=0.3 ; FRAC_DE=0.05 ; FC=1 ; DGE_FUNC="deseq2"

 x<-simrna2d(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,gsets)

 # run the DESeq2 DE analysis
 x<-deseq2(x)

 # run mitch
 x<-run_mitch(x,DGE_FUNC,gsets, N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC)

 # run fgsea
 x<-run_fgsea(x,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC)

 # run mdgsa
# x<-run_mdgsa(x,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC)

 # look at the results
 str(x)

 # run the full shbang
 SIMS=10
 test1<-agg_dge(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,gsets)
}

###############################################
# run simulations over a range of parameters
###############################################
SIMS=10
unlink("simpw2d_res_running.tsv")
res=NULL
for ( FRAC_DE in c(0.05)) {
  for (FC in c(1)) {
    for (N_REPS in c(3)) {
      for (DGE_FUNC in c("deseq2")) {
        for ( SUM_COUNT in c( 40000000 )) {
          for  ( VARIANCE in c(0.3)) {
            x<-agg_dge(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,gsets)
            write.table(x,file="simpw2d_res_running.tsv",quote=F,sep='\t',append=T)
            res=c(res,x)
          }
        }
      }
    }
  }
}
