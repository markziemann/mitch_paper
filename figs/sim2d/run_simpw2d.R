
library("tidyverse")
library("parallel")
library("edgeR")
library("DESeq2")
library("limma")
library("stringi")
library("mdgsa")
library("fgsea")

source("simpw2d_func.R")

# obtain count data
a<-countData()

# generate some gene sets
gsets<-randomGeneSets(a)

test_simpw2d<-function(){
 # create some random data with two contraste
 N_REPS=3 ; SUM_COUNT=10000000 ; VARIANCE=0.3 ; FRAC_DE=0.05 ; FC=1 ; DGE_FUNC="deseq2" ; SIMS=10
 x<-simrna2d(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,gsets)
 # run the DESeq2 DE analysis
 x<-deseq2(x)
 # run mitch
 system.time(x<-run_mitch(x,DGE_FUNC,gsets, N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC))
 # run fgsea
 system.time(x<-run_fgsea(x,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC))
 # run mdgsa
 system.time(x<-run_mdgsa(x,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC))
 # run mavtgsa
 system.time(x<-run_mavtgsa(x,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC) ) #SSLLOOOOWWWW
 # look at the results
 str(x)
 # run the full shbang
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
        for ( SUM_COUNT in c(10000000, 40000000 , 100000000)) {
          for  ( VARIANCE in c(0, 0.3, 0.6, 0.9)) {
            x<-agg_dge(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,gsets)
            write.table(x,file="simpw2d_res_running.tsv",quote=F,sep='\t',append=T)
            res=c(res,x)
          }
        }
      }
    }
  }
}
