library("tidyverse")
library("parallel")
library("DESeq2")
library("stringi")
library("mitch")
library("fgsea")
library("mdgsa")
library("MAVTgsa")
library("edgeR")

########################################
# get some counts
########################################
countData<-function() {
# a is orig expression data
a<-read.table("https://raw.githubusercontent.com/markziemann/simde/master/start_data/ERR2539161/ERR2539161.se.tsv")

# GI is the gene information. It has gene name information for mapping accessions to gene symbols in GMT files
gi<-read.table("https://raw.githubusercontent.com/markziemann/simde/master/start_data/ERR2539161/GeneInfo.tsv",header=T,row.names=1)

# merge gene names
aa<-merge(a,gi,by=0)
aa<-aa[,-c(4:7)]
aa<-aggregate(. ~ GeneSymbol,aa,function(x) sum(as.numeric(as.character(x))))
aa$Row.names=NULL
rownames(aa)<-aa$GeneSymbol
aa$GeneSymbol=NULL
a<-aa[which(aa$ERR2539161>=10),,drop=F]
a
}

########################################
# generate some gene sets
########################################
randomGeneSets<-function(a){
gsets<-sapply( rep(50,1000) , function(x) {list(as.character(sample(rownames(a),x))) } )
names(gsets) <- stri_rand_strings(length(gsets), 15, pattern = "[A-Za-z]")
gsets
}

########################################
# simulate some gene expression data
########################################
simrna2d<-function(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,gsets) {
# x<-simrna2d(a,3,10000000,0.2,0.2,1,gsets)
# N_REPS=3 ; SUM_COUNT=10000000 ; VARIANCE=0.3 ; FRAC_DE=0.1 ; FC=1 
NUM_SAMPLE_GROUPS=3
df = NULL
for (k in paste0("data",1:(N_REPS*3)))  {
        b <- thinCounts(a,target.size = SUM_COUNT)
        colnames(b) = k
        df = cbind(df,b)
     }

# now need to only include gsets with 10 members in the 
gsets_sub <- which(unlist( lapply(gsets,function(x) {
  length(which(rownames(a) %in% as.character(unlist(x)))) >10 
}  ) ) )
gsets <- gsets[which(names(gsets) %in% names(gsets_sub))]
#Number of differential genes
NDIF = round(length(gsets)*FRAC_DE)
# introduce some variance here
if (VARIANCE>0) {
  #create some random values centred around 1 with some% error
  rand <- matrix(log2(rnorm(nrow(a)*N_REPS * NUM_SAMPLE_GROUPS , 2, VARIANCE)),
    ncol = N_REPS * NUM_SAMPLE_GROUPS)
  #incorporate the noise
  df <- round(df*rand)
  #set any negative counts to zero
  df <- apply(df, 2, function(x) {ifelse(x < 0, 0, x)})
}

# prepare some fold changes
if (NDIF>0) {
  message("prep fold changes")
  # Make NDIF a multiple of 4
  if ( NDIF %% 2 == 1 ) { print("odd") ; NDIF = NDIF-1 }
  if ( NDIF %% 4 == 2 ) { print("not multiple of 4") ; NDIF = NDIF-2 }
  # select some sets to be DE
  DE_LIST <- sample(gsets , NDIF)
  # we have gene sets can be in segments 1 to 4 or 6 to 9.
  MYSAMPLE = sample(c(1,2,3,4,6,7,8,9), NDIF, replace = TRUE)
  names(MYSAMPLE) <- names(DE_LIST)
  # define upreg genes in contrast 1
  UP_LIST1 <- DE_LIST[which(MYSAMPLE%%3==0)]
  # now find a list of genes inside the pathways
  UP_DE1 <- unique(unlist(unname(UP_LIST1)))
  # select the ones that are also in the profile
  UP_DE1 <- UP_DE1[which(UP_DE1 %in% row.names(df))]
  # same for down genes
  DN_LIST1 <- DE_LIST[which(MYSAMPLE%%3==1)]
  DN_DE1 <- unique(unlist(unname(DN_LIST1)))
  DN_DE1 <- DN_DE1[which(DN_DE1 %in% row.names(df))]
  ITX <- intersect(UP_DE1,DN_DE1)
  # need to eliminate the overlapping ones for simplicity
  UP_DE1 <- setdiff(UP_DE1,ITX)
  DN_DE1 <- setdiff(DN_DE1,ITX)
  #reformat as df and add fold change
  UP_DE1 <- as.data.frame(UP_DE1)
  UP_DE1$V1 <- 2^FC
  colnames(UP_DE1) = c("Gene","FC")
  DN_DE1 <- as.data.frame(DN_DE1)
  DN_DE1$V1 <- 2^-FC
  colnames(DN_DE1) = c("Gene","FC")
  ALL_DE1 <- rbind(DN_DE1,UP_DE1)
  #Go back to list for downstream work
  UP_DE1 <- UP_DE1$Gene
  DN_DE1 <- DN_DE1$Gene
  NON_DE1 <- as.data.frame(setdiff(rownames(df),ALL_DE1$Gene))
  colnames(NON_DE1) = "Gene"
  NON_DE1$FC = 1
  ALL_DE1 <- rbind(ALL_DE1,NON_DE1)
  ALL_DE1 <- data.frame(ALL_DE1[ order(as.vector(ALL_DE1$Gene)) , ])
  # now DE2
  UP_LIST2 <- DE_LIST[which(MYSAMPLE<4)]
  UP_DE2 <- unique(unlist(unname(UP_LIST2)))
  UP_DE2 <- UP_DE2[which(UP_DE2 %in% row.names(df))]
  DN_LIST2 <- DE_LIST[which(MYSAMPLE>6)]
  DN_DE2 <- unique(unlist(unname(DN_LIST2)))
  DN_DE2 <- DN_DE2[which(DN_DE2 %in% row.names(df))]
  ITX<-intersect(UP_DE2,DN_DE2)
  UP_DE2 <- setdiff(UP_DE2,ITX)
  DN_DE2 <- setdiff(DN_DE2,ITX)
  UP_DE2 <- as.data.frame(UP_DE2)
  UP_DE2$V1 <- 2^FC
  colnames(UP_DE2) = c("Gene","FC")
  DN_DE2 <- as.data.frame(DN_DE2)
  DN_DE2$V1 <- 2^-FC
  colnames(DN_DE2) = c("Gene","FC")
  ALL_DE2 <- rbind(DN_DE2,UP_DE2)
  UP_DE2 <- UP_DE2$Gene
  DN_DE2 <- DN_DE2$Gene
  NON_DE2 <- as.data.frame(setdiff(rownames(df),ALL_DE2$Gene))
  colnames(NON_DE2) = "Gene"
  NON_DE2$FC = 1
  ALL_DE2 <- rbind(ALL_DE2,NON_DE2)
  ALL_DE2 <- data.frame(ALL_DE2[ order(as.vector(ALL_DE2$Gene)) , ])
  message("incorporate changes")
  df <- df[ order(row.names(df)), ]
} else {
  ALL_DE1 <- df[,1,drop=F] ; ALL_DE1[,1] <- 1 ; colnames(ALL_DE1) = "FC"
  ALL_DE2 <- df[,1,drop=F] ; ALL_DE2[,1] <- 1 ; colnames(ALL_DE2) = "FC"
  ALL_DE1 <- data.frame(ALL_DE1)
  ALL_DE2 <- data.frame(ALL_DE2)
  UP_DE1 = NULL ; DN_DE1 = NULL ; UP_LIST1 = NULL ; DN_LIST1 = NULL
  UP_DE2 = NULL ; DN_DE2 = NULL ; UP_LIST2 = NULL ; DN_LIST2 = NULL
}
ONE_IDX_COLS = (1:ncol(df))[c(TRUE,FALSE,FALSE)]
TWO_IDX_COLS = (1:ncol(df))[c(FALSE,TRUE,FALSE)]
THREE_IDX_COLS = (1:ncol(df))[c(FALSE,FALSE,TRUE)]
controls1 <- df[,ONE_IDX_COLS]
colnames(controls1) = paste0( "ctrl1_" ,1:ncol(controls1) )
treatments1 <- round(df[,TWO_IDX_COLS]*ALL_DE1$FC)
colnames(treatments1) = paste0( "trt1_" ,1:ncol(treatments1) )
treatments2 <- round(df[,c(THREE_IDX_COLS)]*ALL_DE2$FC)
colnames(treatments2) = paste0( "trt2_" ,1:ncol(treatments2) )
x <- cbind(controls1,treatments1,treatments2)
rownames(x) = rownames(df)
x <- x[which(rowSums(x)/ncol(x)>10),]
x <- as.data.frame(x)
UP_DE1 <- intersect(UP_DE1,rownames(x)) ; DN_DE1 <- intersect(DN_DE1,rownames(x))
UP_DE2 <- intersect(UP_DE2,rownames(x)) ; DN_DE2 <- intersect(DN_DE2,rownames(x))
xx <- list("x" = x, "truth" = MYSAMPLE)
xx
}

#################################################
# A parallel repeat function
#################################################
#Thanks Gray Calhoun gcalhoun@iastate.edu for the following function
RepParallel <- function(n, expr, simplify = "array",...) {
      answer <-
        mclapply(integer(n), eval.parent(substitute(function(...) expr)),...)
      if (!identical(simplify, FALSE) && length(answer)) 
        return(simplify2array(answer, higher = (simplify == "array")))
      else return(answer)
    }
# RepParallel usage
#xxx<-RepParallel(10,simrna(a,5,10000000,0.2,20), simplify=F, mc.cores = detectCores() )

#################################################
# define DESeq2 2D function
##################################################
deseq2 <- function(x) {
library("DESeq2")
label = "simulate"
y <- x[[1]]
samplesheet <- as.data.frame(colnames(y))
colnames(samplesheet) = "sample"

y1 <- y[,c(grep("ctrl1",samplesheet$sample), grep("trt1",samplesheet$sample)) ]
samplesheet1 <- as.data.frame(colnames(y1))
colnames(samplesheet1) = "sample"
samplesheet1$trt <- as.numeric(grepl("trt",colnames(y1)))
dds <- DESeqDataSetFromMatrix(countData = y1, colData = samplesheet1, design = ~ trt )
res <- DESeq(dds)
z<- DESeq2::results(res)
vsd <- vst(dds, blind=FALSE)
zz <- cbind(z,assay(vsd))
x[[6]] <- as.data.frame(zz[order(zz$padj),])
sig <- subset(zz,padj<0.05)
x[[7]] <- rownames(sig[which(sig$log2FoldChange>0),])
x[[8]] <- rownames(sig[which(sig$log2FoldChange<0),])

y2 <- y[,c(grep("ctrl1",samplesheet$sample), grep("trt2",samplesheet$sample)) ]
samplesheet2 <- as.data.frame(colnames(y2))
colnames(samplesheet2) = "sample"
samplesheet2$trt <- as.numeric(grepl("trt",colnames(y2)))
dds <- DESeqDataSetFromMatrix(countData = y2, colData = samplesheet2, design = ~ trt )
res <- DESeq(dds)
z <- DESeq2::results(res)
vsd <- vst(dds, blind=FALSE)
zz <- cbind(z,assay(vsd))
x[[9]] <- as.data.frame(zz[order(zz$padj),])
sig <- subset(zz,padj<0.05)
x[[10]] <- rownames(sig[which(sig$log2FoldChange>0),])
x[[11]] <- rownames(sig[which(sig$log2FoldChange<0),])
names(x)[6:11] <- c("DGE1","DGE1_up","DGE1_dn","DGE2","DGE2_up","DGE2_dn")
x
}

#################################################
# define mitch function
##################################################
run_mitch <- function(x,DGE_FUNC,gsets, N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC) {
library("mitch")
dge1 <- x[[6]]
dge2 <- x[[9]]
dge <- list("dge1"=dge1,"dge2"=dge2)
w <- mitch_import(dge, DGE_FUNC )
res <- mitch_calc(w,gsets,priority="significance",cores=1)
res <- res$enrichment_result
mitch_sig <- res[which(res$p.adjustMANOVA<0.05),1]
gt <- names(x$truth)
true_pos=length(intersect(mitch_sig , gt ))
false_pos=length(setdiff( mitch_sig, gt  ))
false_neg=length(setdiff(gt , mitch_sig ))
true_neg=length(gsets)-sum(true_pos,false_pos,false_neg)
p<-true_pos/(true_pos+false_pos)
r<-true_pos/(true_pos+false_neg)
f<-2*p*r/(p+r)
attr(x,'mitch_res') <-data.frame(N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,DGE_FUNC,
  true_pos,false_pos,true_neg,false_neg,p,r,f)
x
}

##################################
# FGSEA function
##################################
run_fgsea <- function(x,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC){
dge1 <- x[[6]]
dge2 <- x[[9]]
s1 <- dge1$stat
names(s1) <- rownames(dge1)
s1[is.na(s1)] <- 0
p1 <- as.data.frame(fgsea(pathways=gsets, stats=s1 ,nperm=2000 ))
obs1 <- subset(p1,padj<0.05 )[,1]
s2 <- dge2$stat
names(s2) <- rownames(dge2)
s2[is.na(s2)] <- 0
p2 <- as.data.frame(fgsea(pathways=gsets, stats=s2 ,nperm=2000))
obs2 <- subset(p2,padj<0.05 )[,1]
obs <- union(obs1,obs2)
gt <- names(x$truth)
true_pos = length(intersect( obs , gt ))
false_pos = length(setdiff(obs , gt ))
false_neg = length(setdiff(gt , obs ))
true_neg=length(gsets)-sum(true_pos,false_pos,false_neg)
p<-true_pos/(true_pos+false_pos)
r<-true_pos/(true_pos+false_neg)
f<-2*p*r/(p+r)
attr(x,'fgsea_res') <-data.frame(N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,DGE_FUNC,
  true_pos,false_pos,true_neg,false_neg,p,r,f)
x
}

#################################################
# define mdgsa function
##################################################
run_mdgsa<-function(x,DGE_FUNC,gsets, N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC) {
library("mdgsa")
dge1<-x[[6]]
dge2<-x[[9]]
dge<-list("dge1"=dge1,"dge2"=dge2)
w<-mitch_import(dge, DGE_FUNC )
w[is.na(w)] <- 0
res <- mdGsa (w, gsets)
sig <- rownames(subset(res,padj.dge1 < 0.05 | padj.dge2 < 0.05 | padj.I < 0.05))
gt <- names(x$truth)
true_pos=length(intersect( sig , gt ))
false_pos=length(setdiff( sig, gt  ))
false_neg=length(setdiff(gt , sig ))
true_neg=length(gsets)-sum(true_pos,false_pos,false_neg)
p<-true_pos/(true_pos+false_pos)
r<-true_pos/(true_pos+false_neg)
f<-2*p*r/(p+r)
attr(x,'mdgsa_res') <-data.frame(N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,DGE_FUNC,
  true_pos,false_pos,true_neg,false_neg,p,r,f)
x
}

#################################################
# define MAVTgsa function
##################################################
run_mavtgsa<-function(x,DGE_FUNC,gsets, N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC) {
mx <- x[[1]]
gx <- lapply(gsets , function(x) { as.numeric( rownames (mx ) %in% x ) })
gx <- data.frame(gx)
samplegroups <- as.numeric(factor(sapply(strsplit(colnames(mx),"_"),"[[",1)))
mx <- rbind(t(data.frame(samplegroups)),as.matrix(mx))
res <- MAVTn(mx, gx, alpha = 0.05, nbPerm = 1000 )
sig<- names(gsets)[which(res[[1]]$`MANOVA p-value`<0.05)]
gt <- names(x$truth)
true_pos=length(intersect( sig , gt ))
false_pos=length(setdiff( sig, gt  ))
false_neg=length(setdiff(gt , sig ))
true_neg=length(gsets)-sum(true_pos,false_pos,false_neg)
p<-true_pos/(true_pos+false_pos)
r<-true_pos/(true_pos+false_neg)
f<-2*p*r/(p+r)

attr(x,'mdgsa_res') <-data.frame(N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,DGE_FUNC,
  true_pos,false_pos,true_neg,false_neg,p,r,f)
x
}


##################################
# aggregate function
##################################
agg_dge<-function(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,gsets) {
#SIMS=10

myagg<-function(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,gsets) {
 x<-simrna2d(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,gsets)
 x<-deseq2(x)
# x<-run_mitch(x,DGE_FUNC,gsets, N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC)
 x<-run_fgsea(x,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC)
# x<-run_mdgsa(x,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC)
# x<-run_mavtgsa(x,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC) #SSLLOOOOWWWW

 g=list()
 for (f in 2:length(attributes(x))) {
  PWAY_FUNC<-names(attributes(x)[f])
  PWAY_FUNC<-as.data.frame(PWAY_FUNC)
  g[[f-1]]<-cbind(unname(attributes(x)[f]),PWAY_FUNC)
 }
 g<-as.data.frame(do.call(rbind, g))
 g
}

g<-RepParallel(SIMS,myagg(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,gsets), simplify=F, mc.cores = 10 )
g<-as.data.frame(do.call(rbind, g))
g
}
# N_REPS=3 ; SUM_COUNT=40000000 ; VARIANCE=0 ; FRAC_DE=0.05 ; FC=1 ; SIMS=10 ; DGE_FUNC="deseq2" ; gsets=gsets
# res<-agg_dge(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,gsets)
# res<-agg_dge(a,10,40000000,0.4,0.2,1,10,"deseq2",gsets) 

