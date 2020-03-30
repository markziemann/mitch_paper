library("mitch")
library("fgsea")
library("mdgsa")
library("reshape2")
library("stringi")
library("mdgsa")
library("MAVTgsa")

########################################
# randy makes random gene names
########################################
randy <- function(n = 5000) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

########################################
# generate some gene sets # don't use
########################################
randomGeneSets<-function(a){
gsets<-sapply( rep(50,1000) , function(x) {list(as.character(sample(rownames(a),x))) } )
names(gsets)<-stri_rand_strings(length(gsets), 15, pattern = "[A-Za-z]")
gsets
}

# reverse the row order of a dataframe
rev_df_row<-function(df){
    df[seq(dim(df)[1],1),]
}

#################################################
# first example - compare mitch and fgsea
#   cores = 1 2 4 8 16
#   genes = 20000
#   1 dimensions
#   1000 sets of genes
#   500 genes per set
#################################################
genes<-randy(20000)
k1<-as.data.frame(genes,stringsAsFactors=F)

k1$a<-rnorm(20000, 0, 3)

rownames(k1)<-k1[,1]
k1[,1]=NULL

# create random gene sets
s_500_1000<-sapply(rep(500,1000), function(x) {list(as.character(sample(genes,x))) } )
names(s_500_1000)<-paste("set",1:1000)

# run mitch with different cores
m_c1<-system.time( res<-mitch_calc(k1,s_500_1000,priority="significance",cores=1) )
m_c2<-system.time( res<-mitch_calc(k1,s_500_1000,priority="significance",cores=2) )
m_c4<-system.time( res<-mitch_calc(k1,s_500_1000,priority="significance",cores=4) )
m_c8<-system.time( res<-mitch_calc(k1,s_500_1000,priority="significance",cores=8) )
m_c16<-system.time( res<-mitch_calc(k1,s_500_1000,priority="significance",cores=16) )

times_m<-as.data.frame(rbind(
  m_c16,
  m_c8,
  m_c4,
  m_c2,
  m_c1))

cores<-c(16,8,4,2,1)

times_m$cores<-cores
times_m<-rev_df_row(times_m)

# run fgsea with different cores
s<-k1$a
names(s)<-rownames(k1)
f_c16<-system.time(fgsea(pathways=s_500_1000, stats=s, nperm=1000, nproc=16))
f_c8<-system.time(fgsea(pathways=s_500_1000, stats=s, nperm=1000, nproc=8))
f_c4<-system.time(fgsea(pathways=s_500_1000, stats=s, nperm=1000, nproc=4))
f_c2<-system.time(fgsea(pathways=s_500_1000, stats=s, nperm=1000, nproc=2))
f_c1<-system.time(fgsea(pathways=s_500_1000, stats=s, nperm=1000, nproc=1))

times_f<-as.data.frame(rbind(
  f_c16,
  f_c8,
  f_c4,
  f_c2,
  f_c1))

times_f$cores<-cores
times_f<-rev_df_row(times_f)

# data
write.table(times_f,file="fig7a_fgsea.tsv",sep="\t",quote=F,row.names=TRUE)
write.table(times_m,file="fig7a_mitch.tsv",sep="\t",quote=F,row.names=TRUE)

times_f<-read.table("fig7a_fgsea.tsv")
times_m<-read.table("fig7a_mitch.tsv")

pdf("fig7a.pdf")
matplot(times_m$elapsed , pch=1,type = c("b"),xlab="parallel CPU threads",ylab="elapsed time (s)" , 
  main="Execution time",axes=F , ylim=c(0,10))
points(times_f$elapsed, type="b",col="red")
grid()
axis(2)
legend("topright",inset=0.1, legend = c("Mitch","FGSEA"), col=c("black","red"), pch=1) # optional legend
mtext("1 dimension, 20000 genes, 1000 sets, 50 genes per set")
axis(side=1,at=1:nrow(times_m),labels=times_m$cores)
dev.off()

#################################################
# second example compare 2xfgsea mitch mdgsa and MAVtsa
#   genes = 20000
#   dimensions = 2,3,5
#################################################

# create the large dataset
g50k <- randy(50000)
k <- as.data.frame(g50k,stringsAsFactors=F)
k$a1 <- rnorm(50000, 0, 3)
k$a2 <- rnorm(50000, 0, 3)
k$a3 <- rnorm(50000, 0, 3)
k$a4 <- rnorm(50000, 0, 3)
k$a5 <- rnorm(50000, 0, 3)

rownames(k) <- k[,1]
k[,1] = NULL

k2_g20k <- k[1:20000,1:2]
k3_g20k <- k[1:20000,1:3]
k5_g20k <- k[1:20000,1:5]

s <- sapply(rep(50,1000), function(x) {list(as.character(sample(rownames(k2_g20k),x))) } )
names(s)<-paste("set",1:1000)

# mitch
mi_k2 <- system.time( res<-mitch_calc(k2_g20k,s,priority="significance",cores=8) )
mi_k3 <- system.time( res<-mitch_calc(k3_g20k,s,priority="significance",cores=8) )
mi_k5 <- system.time( res<-mitch_calc(k5_g20k,s,priority="significance",cores=8) )

times_mi <- as.data.frame(rbind(mi_k2,
mi_k3,
mi_k5))

times_mi$dims<-c(2,3,5)

# fgsea
s1 <- k5_g20k[,1] ; names(s1)<-rownames(k2_g20k)
s2 <- k5_g20k[,2] ; names(s2)<-rownames(k2_g20k)
s3 <- k5_g20k[,3] ; names(s3)<-rownames(k2_g20k)
s4 <- k5_g20k[,4] ; names(s4)<-rownames(k2_g20k)
s5 <- k5_g20k[,5] ; names(s5)<-rownames(k2_g20k)

f_k2_a <- system.time( fgsea(pathways=s, stats=s1, nperm=1000, nproc=8) )
f_k2_b <- system.time(fgsea(pathways=s, stats=s2, nperm=1000, nproc=8) )
f_k2 = f_k2_a + f_k2_b

f_k3_a <- system.time( fgsea(pathways=s, stats=s1, nperm=1000, nproc=8) )
f_k3_b <- system.time( fgsea(pathways=s, stats=s2, nperm=1000, nproc=8) )
f_k3_c <- system.time( fgsea(pathways=s, stats=s3, nperm=1000, nproc=8) )
f_k3 = f_k3_a + f_k3_b + f_k3_c

f_k5_a <- system.time( fgsea(pathways=s, stats=s1, nperm=1000, nproc=8) )
f_k5_b <- system.time( fgsea(pathways=s, stats=s2, nperm=1000, nproc=8) )
f_k5_c <- system.time( fgsea(pathways=s, stats=s3, nperm=1000, nproc=8) )
f_k5_d <- system.time( fgsea(pathways=s, stats=s4, nperm=1000, nproc=8) )
f_k5_e <- system.time( fgsea(pathways=s, stats=s5, nperm=1000, nproc=8) )
f_k5 <- f_k5_a + f_k5_b + f_k5_c + f_k5_d + f_k5_e

times_f<-as.data.frame(rbind(
  f_k2,
  f_k3,
  f_k5))

times_f$dims<-c(2,3,5)

# mdgsa
md_k2 <- system.time( res <- mdGsa (k2_g20k, s) )
md_k3 <- system.time( res <- mdGsa (k3_g20k, s) )
md_k5 <- system.time( res <- mdGsa (k5_g20k, s) )

times_md <- as.data.frame(rbind(md_k2,
md_k3,
md_k5))

times_md$dims<-c(2,3,5)

# MAVTgsa
mx = NULL
for (k in paste0("data", 1:18 ))  {
        b <- as.data.frame( sample(1:nrow(k1),nrow(k1),replace=TRUE) )
        #b <- thinCounts(a,target.size = SUM_COUNT)
        colnames(b) = k
        if ( is.null(mx) ) { mx <- b }
        else { mx = cbind(mx,b) }
     }

rownames(mx) <- rownames(k1)
s <- sapply(rep(50,1000), function(x) {list(as.character(sample(rownames(mx),x))) } )
names(s)<-paste("set",1:1000)
gx <- lapply(s , function(x) { as.numeric( rownames (mx ) %in% x ) })
gx <- data.frame(gx)
samplegroups <- as.vector(sapply(1:6 ,function(x) { rep(x,3) } ))
mx <- rbind(t(data.frame(samplegroups)),as.matrix(mx))

mx2 <- mx[,1:9]
mx3 <- mx[,1:12]
mx5 <- mx
ma_k2 <- system.time(  res <- MAVTn(mx2, gx, alpha = 0.05, nbPerm = 1000 ) )
ma_k3 <- system.time(  res <- MAVTn(mx3, gx, alpha = 0.05, nbPerm = 1000 ) )
ma_k5 <- system.time(  res <- MAVTn(mx5, gx, alpha = 0.05, nbPerm = 1000 ) )

times_ma <- as.data.frame(rbind(ma_k5,
ma_k3,
ma_k2))

times_ma$dims<-c(2,3,5)

# data
write.table(times_f,file="fig7x_fgsea.tsv",sep="\t",quote=F,row.names=TRUE)
write.table(times_mi,file="fig7x_mitch.tsv",sep="\t",quote=F,row.names=TRUE)
write.table(times_md,file="fig7x_mdgsa.tsv",sep="\t",quote=F,row.names=TRUE)
write.table(times_ma,file="fig7x_MAVTgsa.tsv",sep="\t",quote=F,row.names=TRUE)

times_f <- read.table("fig7x_fgsea.tsv")
times_mi <- read.table("fig7x_mitch.tsv")
times_md <- read.table("fig7x_mdgsa.tsv")
times_ma <- read.table("fig7x_MAVTgsa.tsv")



pdf("fig7x.pdf")
matplot(times_md$elapsed , pch=19 ,type = c("b"),xlab="no. contrasts",ylab="elapsed time (s)" ,
  main="Execution time",axes=F , ylim=c(0.1,11000), col="blue" , log="y")
points(times_ma$elapsed, type="b",col="darkgreen",pch=19)
points(times_f$elapsed, type="b",col="red", pch=19)
points(times_mi$elapsed, type="b",col="black",pch=19)
grid()
axis(2)
legend("topleft",inset=0.1, legend = c("MAVTgsa","mdgsa","mitch","FGSEA"),
 col=c("darkgreen","blue","black","red"), pch=1) # optional legend
mtext("20000 genes, 1000 sets, 50 genes per set")
axis(side=1,at=1:nrow(times_mi),labels=times_mi$dims)
dev.off()

#################################################
# third example - vary the no genes in profile and number of dimensions
#   genes = 1000 5000 20000 100000 
#   dimensions = 3,5,10,15
#   cores = 8
#   bootstraps = 500
#   no genesets = 1000 sets of genes
#   genes per set = 500
#################################################

# create the large dataset
g50k<-randy(50000)

k<-as.data.frame(g50k,stringsAsFactors=F)

k$a<-rnorm(50000, 0, 3)
k$a2<-rnorm(50000, 0, 3)
k$a3<-rnorm(50000, 0, 3)
k$a4<-rnorm(50000, 0, 3)
k$a5<-rnorm(50000, 0, 3)
k$a6<-rnorm(50000, 0, 3)
k$a7<-rnorm(50000, 0, 3)
k$a8<-rnorm(50000, 0, 3)
k$a9<-rnorm(50000, 0, 3)
k$a10<-rnorm(50000, 0, 3)
k$a11<-rnorm(50000, 0, 3)
k$a12<-rnorm(50000, 0, 3)
k$a13<-rnorm(50000, 0, 3)
k$a14<-rnorm(50000, 0, 3)
k$a15<-rnorm(50000, 0, 3)

rownames(k)<-k[,1]
k[,1]=NULL

k15_g50k<-k
k15_g20k<-k[1:20000,]
k15_g5k<-k[1:5000,]
k15_g1k<-k[1:1000,]
k10_g50k<-k[,1:10]
k10_g20k<-k[1:20000,1:10]
k10_g5k<-k[1:5000,1:10]
k10_g1k<-k[1:1000,1:10]
k5_g50k<-k[,1:5]
k5_g20k<-k[1:20000,1:5]
k5_g5k<-k[1:5000,1:5]
k5_g1k<-k[1:1000,1:5]
k3_g50k<-k[,1:3]
k3_g20k<-k[1:20000,1:3]
k3_g5k<-k[1:5000,1:3]
k3_g1k<-k[1:1000,1:3]
k1_g50k<-k[,1:3]
k1_g20k<-k[1:20000,1:3]
k1_g5k<-k[1:5000,1:3]
k1_g1k<-k[1:1000,1:3]


# k15 
s<-sapply(rep(50,1000), function(x) {list(as.character(sample(rownames(k15_g50k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k15_g50k<-system.time( res<-mitch_calc(k15_g50k,s,priority="significance",cores=8) )

s<-sapply(rep(50,1000), function(x) {list(as.character(sample(rownames(k15_g20k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k15_g20k<-system.time( res<-mitch_calc(k15_g20k,s,priority="significance",cores=8) )

s<-sapply(rep(50,1000), function(x) {list(as.character(sample(rownames(k15_g5k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k15_g5k<-system.time( res<-mitch_calc(k15_g5k,s,priority="significance",cores=8) )

s<-sapply(rep(50,1000), function(x) {list(as.character(sample(rownames(k15_g1k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k15_g1k<-system.time( res<-mitch_calc(k15_g1k,s,priority="significance",cores=8) )

# k10
s<-sapply(rep(50,1000), function(x) {list(as.character(sample(rownames(k10_g50k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k10_g50k<-system.time( res<-mitch_calc(k10_g50k,s,priority="significance",cores=8) )

s<-sapply(rep(50,1000), function(x) {list(as.character(sample(rownames(k10_g20k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k10_g20k<-system.time( res<-mitch_calc(k10_g20k,s,priority="significance",cores=8) )

s<-sapply(rep(50,1000), function(x) {list(as.character(sample(rownames(k10_g5k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k10_g5k<-system.time( res<-mitch_calc(k10_g5k,s,priority="significance",cores=8) )

s<-sapply(rep(50,1000), function(x) {list(as.character(sample(rownames(k10_g1k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k10_g1k<-system.time( res<-mitch_calc(k10_g1k,s,priority="significance",cores=8) )

# k5
s<-sapply(rep(50,1000), function(x) {list(as.character(sample(rownames(k5_g50k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k5_g50k<-system.time( res<-mitch_calc(k5_g50k,s,priority="significance",cores=8) )

s<-sapply(rep(50,1000), function(x) {list(as.character(sample(rownames(k5_g20k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k5_g20k<-system.time( res<-mitch_calc(k5_g20k,s,priority="significance",cores=8) )

s<-sapply(rep(50,1000), function(x) {list(as.character(sample(rownames(k5_g5k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k5_g5k<-system.time( res<-mitch_calc(k5_g5k,s,priority="significance",cores=8) )

s<-sapply(rep(50,1000), function(x) {list(as.character(sample(rownames(k5_g1k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k5_g1k<-system.time( res<-mitch_calc(k5_g1k,s,priority="significance",cores=8) )

# k3
s<-sapply(rep(50,1000), function(x) {list(as.character(sample(rownames(k3_g50k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k3_g50k<-system.time( res<-mitch_calc(k3_g50k,s,priority="significance",cores=8) )

s<-sapply(rep(50,1000), function(x) {list(as.character(sample(rownames(k3_g20k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k3_g20k<-system.time( res<-mitch_calc(k3_g20k,s,priority="significance",cores=8) )

s<-sapply(rep(50,1000), function(x) {list(as.character(sample(rownames(k3_g5k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k3_g5k<-system.time( res<-mitch_calc(k3_g5k,s,priority="significance",cores=8) )

s<-sapply(rep(50,1000), function(x) {list(as.character(sample(rownames(k3_g1k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k3_g1k<-system.time( res<-mitch_calc(k3_g1k,s,priority="significance",cores=8) )

times_b<-as.data.frame(rbind(t_k15_g50k,
t_k15_g20k,
t_k15_g5k,
t_k15_g1k,
t_k10_g50k,
t_k10_g20k,
t_k10_g5k,
t_k10_g1k,
t_k5_g50k,
t_k5_g20k,
t_k5_g5k,
t_k5_g1k,
t_k3_g50k,
t_k3_g20k,
t_k3_g5k,
t_k3_g1k))

dims<-c(15,15,15,15,10,10,10,10,5,5,5,5,3,3,3,3)

ngenes<-rep(c(100000,20000,5000,1000),4)

times_b$dims<-dims
times_b$ngenes<-ngenes

times_b_wide<-t(acast(times_b,ngenes ~ dims, value.var="elapsed"))

write.table(times_b_wide,file="fig7b.tsv")
times_b_wide<-read.table("fig7b.tsv")

pdf("fig7b.pdf")
matplot(times_b_wide,pch=1,type = c("b"),xlab="dimensions",ylab="elapsed time (s)" ,
  main="variation of n genes and dimensions", axes=F, ylim=c(0,20))
grid()
axis(2)
axis(side=1,at=1:nrow(times_b_wide),labels=rownames(times_b_wide))
legend("topleft", inset=0.1, legend = c("50k","20k","5k","1k"), title="no. genes in profile", col=4:1, pch=1)
mtext("1000 sets, 50 genes per set, 8 parallel CPU threads")
dev.off()



#################################################
# next, see vary the number and size of gene sets
#   genesets = 100 500 5000 20000
#   geneset size = 20 100 500 2000
#   cores = 8
#   bootstraps = 500 
#   genes = 20000
#   5 dimensions
#################################################

genes<-randy(20000)
k1<-as.data.frame(genes,stringsAsFactors=F)

k1$a<-rnorm(20000, 0, 3)
k1$a2<-rnorm(20000, 0, 3)
k1$a3<-rnorm(20000, 0, 3)
k1$a4<-rnorm(20000, 0, 3)
k1$a5<-rnorm(20000, 0, 3)

rownames(k1)<-k1[,1]
k1[,1]=NULL

# no genesets=100 now vay the set size
s_20_100<-sapply(rep(20,100), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_20_100)<-paste("set",1:100)
t_20g_100l<-system.time( res<-mitch_calc(k1,s_20_100,priority="significance",cores=8) )

s_100_100<-sapply(rep(100,100), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_100_100)<-paste("set",1:100)
t_100g_100l<-system.time( res<-mitch_calc(k1,s_100_100,priority="significance",cores=8) )

s_500_100<-sapply(rep(500,100), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_500_100)<-paste("set",1:100)
t_500g_100l<-system.time( res<-mitch_calc(k1,s_500_100,priority="significance",cores=8) )

s_2000_100<-sapply(rep(2000,100), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_2000_100)<-paste("set",1:100)
t_2000g_100l<-system.time( res<-mitch_calc(k1,s_2000_100,priority="significance",cores=8) )


# no genesets=500 now vay the set size
s_20_500<-sapply(rep(20,500), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_20_500)<-paste("set",1:500)
t_20g_500l<-system.time( res<-mitch_calc(k1,s_20_500,priority="significance",cores=8) )

s_100_500<-sapply(rep(100,500), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_100_500)<-paste("set",1:500)
t_100g_500l<-system.time( res<-mitch_calc(k1,s_100_500,priority="significance",cores=8) )

s_500_500<-sapply(rep(500,500), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_500_500)<-paste("set",1:500)
t_500g_500l<-system.time( res<-mitch_calc(k1,s_500_500,priority="significance",cores=8) )

s_2000_500<-sapply(rep(2000,500), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_2000_500)<-paste("set",1:500)
t_2000g_500l<-system.time( res<-mitch_calc(k1,s_2000_500,priority="significance",cores=8) )


# no genesets=5000 now vay the set size
s_20_5000<-sapply(rep(20,5000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_20_5000)<-paste("set",1:5000)
t_20g_5000l<-system.time( res<-mitch_calc(k1,s_20_5000,priority="significance",cores=8) )

s_100_5000<-sapply(rep(100,5000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_100_5000)<-paste("set",1:5000)
t_100g_5000l<-system.time( res<-mitch_calc(k1,s_100_5000,priority="significance",cores=8) )

s_500_5000<-sapply(rep(500,5000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_500_5000)<-paste("set",1:5000)
t_500g_5000l<-system.time( res<-mitch_calc(k1,s_500_5000,priority="significance",cores=8) )

s_2000_5000<-sapply(rep(2000,5000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_2000_5000)<-paste("set",1:5000)
t_2000g_5000l<-system.time( res<-mitch_calc(k1,s_2000_5000,priority="significance",cores=8) )


# no genesets = 20000 now vay the set size
s_20_20000<-sapply(rep(20,20000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_20_20000)<-paste("set",1:20000)
t_20g_20000l<-system.time( res<-mitch_calc(k1,s_20_20000,priority="significance",cores=8) )

s_100_20000<-sapply(rep(100,20000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_100_20000)<-paste("set",1:20000)
t_100g_20000l<-system.time( res<-mitch_calc(k1,s_100_20000,priority="significance",cores=8) )

s_500_20000<-sapply(rep(500,20000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_500_20000)<-paste("set",1:20000)
t_500g_20000l<-system.time( res<-mitch_calc(k1,s_500_20000,priority="significance",cores=8) )

s_2000_20000<-sapply(rep(2000,20000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_2000_20000)<-paste("set",1:20000)
t_2000g_20000l<-system.time( res<-mitch_calc(k1,s_2000_20000,priority="significance",cores=8) )

times_c<-as.data.frame(rbind(t_20g_100l,
t_100g_100l,
t_500g_100l,
t_2000g_100l,
t_20g_500l,
t_100g_500l,
t_500g_500l,
t_2000g_500l,
t_20g_5000l,
t_100g_5000l,
t_500g_5000l,
t_2000g_5000l,
t_20g_20000l,
t_100g_20000l,
t_500g_20000l,
t_2000g_20000l))

nsets<-c(100,100,100,100,500,500,500,500,5000,5000,5000,5000,20000,20000,20000,20000)

setsize<-rep(c(20,100,500,2000),4)

times_c$nsets<-nsets
times_c$ngenes<-setsize

times_c_wide<-t(acast(times_c,nsets ~ setsize, value.var="elapsed"))
write.table(times_c_wide,file="fig7c.tsv")
times_c_wide<-read.table("fig7c.tsv")

pdf("fig7c.pdf")
matplot(times_c_wide,pch=1,type = c("b"),xlab="gene set size",ylab="elapsed time (s)" , ylim=c(0,120) , 
    main="variation of number and size of gene sets", axes=F)
grid()
axis(2)
axis(side=1,at=1:nrow(times_c_wide),labels=rownames(times_c_wide))
legend("topleft", inset=0.1 ,legend = c("20000","5000","500","100"), title="no. gene sets",col=4:1, pch=1)
mtext("5 dimensions, 20000 genes, 8 parallel CPU threads")
dev.off()

save.image("fig7.RData")
