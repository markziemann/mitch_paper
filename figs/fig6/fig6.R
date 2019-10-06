library("mitch")
library("reshape2")
# randy makes random gene names
randy <- function(n = 5000) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

#################################################
# first example - vary the cores and bootstraps with the following prototypical example
#   cores = 1 2 4 8 16
#   bootstraps = 0 100 500 1000
#   genes = 20000
#   5 dimensions
#   1000 sets of genes
#   500 genes per set
#################################################

figa<-function(){
genes<-randy(20000)
k1<-as.data.frame(genes,stringsAsFactors=F)

k1$a<-rnorm(20000, 0, 3)
k1$a2<-rnorm(20000, 0, 3)
k1$a3<-rnorm(20000, 0, 3)
k1$a4<-rnorm(20000, 0, 3)
k1$a5<-rnorm(20000, 0, 3)

rownames(k1)<-k1[,1]
k1[,1]=NULL

# create random gene sets
s_500_1000<-sapply(rep(500,1000), function(x) {list(as.character(sample(genes,x))) } )
names(s_500_1000)<-paste("set",1:1000)

t_c16_b1000<-system.time( res<-mitch_calc(k1,s_500_1000,bootstraps=1000,priority="significance",cores=16) )
t_c16_b500<-system.time( res<-mitch_calc(k1,s_500_1000,bootstraps=500,priority="significance",cores=16) )
t_c16_b100<-system.time( res<-mitch_calc(k1,s_500_1000,bootstraps=100,priority="significance",cores=16) )
t_c16_b0<-system.time( res<-mitch_calc(k1,s_500_1000,priority="significance",cores=16) )

t_c8_b1000<-system.time( res<-mitch_calc(k1,s_500_1000,bootstraps=1000,priority="significance",cores=8) )
t_c8_b500<-system.time( res<-mitch_calc(k1,s_500_1000,bootstraps=500,priority="significance",cores=8) )
t_c8_b100<-system.time( res<-mitch_calc(k1,s_500_1000,bootstraps=100,priority="significance",cores=8) )
t_c8_b0<-system.time( res<-mitch_calc(k1,s_500_1000,priority="significance",cores=8) )

t_c4_b1000<-system.time( res<-mitch_calc(k1,s_500_1000,bootstraps=1000,priority="significance",cores=4) )
t_c4_b500<-system.time( res<-mitch_calc(k1,s_500_1000,bootstraps=500,priority="significance",cores=4) )
t_c4_b100<-system.time( res<-mitch_calc(k1,s_500_1000,bootstraps=100,priority="significance",cores=4) )
t_c4_b0<-system.time( res<-mitch_calc(k1,s_500_1000,priority="significance",cores=4) )

t_c2_b1000<-system.time( res<-mitch_calc(k1,s_500_1000,bootstraps=1000,priority="significance",cores=2) )
t_c2_b500<-system.time( res<-mitch_calc(k1,s_500_1000,bootstraps=500,priority="significance",cores=2) )
t_c2_b100<-system.time( res<-mitch_calc(k1,s_500_1000,bootstraps=100,priority="significance",cores=2) )
t_c2_b0<-system.time( res<-mitch_calc(k1,s_500_1000,priority="significance",cores=2) )

t_c1_b1000<-system.time( res<-mitch_calc(k1,s_500_1000,bootstraps=1000,priority="significance",cores=1) )
t_c1_b500<-system.time( res<-mitch_calc(k1,s_500_1000,bootstraps=500,priority="significance",cores=1) )
t_c1_b100<-system.time( res<-mitch_calc(k1,s_500_1000,bootstraps=100,priority="significance",cores=1) )
t_c1_b0<-system.time( res<-mitch_calc(k1,s_500_1000,priority="significance",cores=1) )

times_a<-as.data.frame(rbind(t_c16_b1000,
  t_c16_b500,
  t_c16_b100,
  t_c16_b0,
  t_c8_b1000,
  t_c8_b500,
  t_c8_b100,
  t_c8_b0,
  t_c4_b1000,
  t_c4_b500,
  t_c4_b100,
  t_c4_b0,
  t_c2_b1000,
  t_c2_b500,
  t_c2_b100,
  t_c2_b0,
  t_c1_b1000,
  t_c1_b500,
  t_c1_b100,
  t_c1_b0))

cores<-c(16,16,16,16,8,8,8,8,4,4,4,4,2,2,2,2,1,1,1,1)

bootstraps<-rep(c(1000,500,100,0),5)

times_a$cores<-cores
times_a$bootstraps<-bootstraps

times_a_wide<-t(acast(times_a,bootstraps ~ cores, value.var="elapsed"))

pdf("fig6a.pdf")
matplot(times_a_wide,pch=1,type = c("b"),xlab="Cores",ylab="elapsed time (s)" , main="variation of cores and bootstraps",axes=F,ylim=c(0,900))
axis(2)
legend("topright", legend = c(1000,500,100,0), col=4:1, pch=1) # optional legend
mtext("5 dimensions, 20000 genes, 1000 sets, 500 genes per set")
axis(side=1,at=1:nrow(times_a_wide),labels=rownames(times_a_wide))

matplot(times_a_wide[,1],pch=1,type = c("b"),xlab="Cores",ylab="elapsed time (s)" , main="variation of cores without bootstraps",axes=F,ylim=c(0,18))
axis(2)
axis(side=1,at=1:nrow(times_a_wide),labels=rownames(times_a_wide))
mtext("5 dimensions, 20000 genes, 1000 sets, 500 genes per set")
dev.off()

write.table(times_a_wide,files="fig6a.tsv")

times_a_wide
}

#################################################
# second example - vary the no genes in profile and number of dimensions
#   genes = 1000 5000 20000 100000 
#   dimensions = 3,5,10,15
#   cores = 8
#   bootstraps = 500
#   no genesets = 1000 sets of genes
#   genes per set = 500
#################################################

# create the large dataset
figb<-function(){
g100k<-randy(100000)

k<-as.data.frame(g100k,stringsAsFactors=F)

k$a<-rnorm(100000, 0, 3)
k$a2<-rnorm(100000, 0, 3)
k$a3<-rnorm(100000, 0, 3)
k$a4<-rnorm(100000, 0, 3)
k$a5<-rnorm(100000, 0, 3)
k$a6<-rnorm(100000, 0, 3)
k$a7<-rnorm(100000, 0, 3)
k$a8<-rnorm(100000, 0, 3)
k$a9<-rnorm(100000, 0, 3)
k$a10<-rnorm(100000, 0, 3)
k$a11<-rnorm(100000, 0, 3)
k$a12<-rnorm(100000, 0, 3)
k$a13<-rnorm(100000, 0, 3)
k$a14<-rnorm(100000, 0, 3)
k$a15<-rnorm(100000, 0, 3)

rownames(k)<-k[,1]
k[,1]=NULL

k15_g100k<-k
k15_g20k<-k[1:20000,]
k15_g5k<-k[1:5000,]
k15_g1k<-k[1:1000,]
k10_g100k<-k[,1:10]
k10_g20k<-k[1:20000,1:10]
k10_g5k<-k[1:5000,1:10]
k10_g1k<-k[1:1000,1:10]
k5_g100k<-k[,1:5]
k5_g20k<-k[1:20000,1:5]
k5_g5k<-k[1:5000,1:5]
k5_g1k<-k[1:1000,1:5]
k3_g100k<-k[,1:3]
k3_g20k<-k[1:20000,1:3]
k3_g5k<-k[1:5000,1:3]
k3_g1k<-k[1:5000,1:3]


# k15 
s<-sapply(rep(500,1000), function(x) {list(as.character(sample(rownames(k15_g100k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k15_g100k<-system.time( res<-mitch_calc(k15_g100k,s,bootstraps=500,priority="significance",cores=8) )

s<-sapply(rep(500,1000), function(x) {list(as.character(sample(rownames(k15_g20k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k15_g20k<-system.time( res<-mitch_calc(k15_g20k,s,bootstraps=500,priority="significance",cores=8) )

s<-sapply(rep(500,1000), function(x) {list(as.character(sample(rownames(k15_g5k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k15_g5k<-system.time( res<-mitch_calc(k15_g5k,s,bootstraps=500,priority="significance",cores=8) )

s<-sapply(rep(500,1000), function(x) {list(as.character(sample(rownames(k15_g1k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k15_g1k<-system.time( res<-mitch_calc(k15_g1k,s,bootstraps=500,priority="significance",cores=8) )

# k10
s<-sapply(rep(500,1000), function(x) {list(as.character(sample(rownames(k10_g100k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k10_g100k<-system.time( res<-mitch_calc(k10_g100k,s,bootstraps=500,priority="significance",cores=8) )

s<-sapply(rep(500,1000), function(x) {list(as.character(sample(rownames(k10_g20k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k10_g20k<-system.time( res<-mitch_calc(k10_g20k,s,bootstraps=500,priority="significance",cores=8) )

s<-sapply(rep(500,1000), function(x) {list(as.character(sample(rownames(k10_g5k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k10_g5k<-system.time( res<-mitch_calc(k10_g5k,s,bootstraps=500,priority="significance",cores=8) )

s<-sapply(rep(500,1000), function(x) {list(as.character(sample(rownames(k10_g1k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k10_g1k<-system.time( res<-mitch_calc(k10_g1k,s,bootstraps=500,priority="significance",cores=8) )

# k5
s<-sapply(rep(500,1000), function(x) {list(as.character(sample(rownames(k5_g100k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k5_g100k<-system.time( res<-mitch_calc(k5_g100k,s,bootstraps=500,priority="significance",cores=8) )

s<-sapply(rep(500,1000), function(x) {list(as.character(sample(rownames(k5_g20k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k5_g20k<-system.time( res<-mitch_calc(k5_g20k,s,bootstraps=500,priority="significance",cores=8) )

s<-sapply(rep(500,1000), function(x) {list(as.character(sample(rownames(k5_g5k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k5_g5k<-system.time( res<-mitch_calc(k5_g5k,s,bootstraps=500,priority="significance",cores=8) )

s<-sapply(rep(500,1000), function(x) {list(as.character(sample(rownames(k5_g1k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k5_g1k<-system.time( res<-mitch_calc(k5_g1k,s,bootstraps=500,priority="significance",cores=8) )

# k3
s<-sapply(rep(500,1000), function(x) {list(as.character(sample(rownames(k3_g100k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k3_g100k<-system.time( res<-mitch_calc(k3_g100k,s,bootstraps=500,priority="significance",cores=8) )

s<-sapply(rep(500,1000), function(x) {list(as.character(sample(rownames(k3_g20k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k3_g20k<-system.time( res<-mitch_calc(k3_g20k,s,bootstraps=500,priority="significance",cores=8) )

s<-sapply(rep(500,1000), function(x) {list(as.character(sample(rownames(k3_g5k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k3_g5k<-system.time( res<-mitch_calc(k3_g5k,s,bootstraps=500,priority="significance",cores=8) )

s<-sapply(rep(500,1000), function(x) {list(as.character(sample(rownames(k3_g1k),x))) } ) ; names(s)<-paste("set",1:1000)
t_k3_g1k<-system.time( res<-mitch_calc(k3_g1k,s,bootstraps=500,priority="significance",cores=8) )


times_b<-as.data.frame(rbind(t_k15_g100k,
t_k15_g20k,
t_k15_g5k,
t_k15_g1k,
t_k10_g100k,
t_k10_g20k,
t_k10_g5k,
t_k10_g1k,
t_k5_g100k,
t_k5_g20k,
t_k5_g5k,
t_k5_g1k,
t_k3_g100k,
t_k3_g20k,
t_k3_g5k,
t_k3_g1k))

dims<-c(15,15,15,15,10,10,10,10,5,5,5,5,3,3,3,3)

ngenes<-rep(c(100000,20000,5000,1000),4)

times_b$dims<-dims
times_b$ngenes<-ngenes

times_b_wide<-t(acast(times_b,ngenes ~ dims, value.var="elapsed"))

pdf("fig6b.pdf")
matplot(times_b_wide,pch=1,type = c("b"),xlab="dimensions",ylab="elapsed time (s)" , main="variation of n genes and dimensions", axes=F, ylim=c(0,500))
grid()
axis(2)
axis(side=1,at=1:nrow(times_b_wide),labels=rownames(times_b_wide))
legend("right", legend = c("100k","20k","5k","1k"), col=4:1, pch=1)
mtext("8 cores, 500 bootstraps, 1000 sets, 500 genes per set")
dev.off()

write.table(times_b_wide,files="fig6b.tsv")

times_b_wide
}

#################################################
# next, see vary the number and size of gene sets
#   genesets = 100 500 5000 20000
#   geneset size = 20 100 500 2000
#   cores = 8
#   bootstraps = 500 
#   genes = 20000
#   5 dimensions
#################################################

figc<-function(){
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
t_20g_100l<-system.time( res<-mitch_calc(k1,s_20_100,bootstraps=0,priority="significance",cores=8) )

s_100_100<-sapply(rep(100,100), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_100_100)<-paste("set",1:100)
t_100g_100l<-system.time( res<-mitch_calc(k1,s_100_100,bootstraps=0,priority="significance",cores=8) )

s_500_100<-sapply(rep(500,100), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_500_100)<-paste("set",1:100)
t_500g_100l<-system.time( res<-mitch_calc(k1,s_500_100,bootstraps=0,priority="significance",cores=8) )

s_2000_100<-sapply(rep(2000,100), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_2000_100)<-paste("set",1:100)
t_2000g_100l<-system.time( res<-mitch_calc(k1,s_2000_100,bootstraps=0,priority="significance",cores=8) )


# no genesets=500 now vay the set size
s_20_500<-sapply(rep(20,500), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_20_500)<-paste("set",1:500)
t_20g_500l<-system.time( res<-mitch_calc(k1,s_20_500,bootstraps=0,priority="significance",cores=8) )

s_100_500<-sapply(rep(100,500), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_100_500)<-paste("set",1:500)
t_100g_500l<-system.time( res<-mitch_calc(k1,s_100_500,bootstraps=0,priority="significance",cores=8) )

s_500_500<-sapply(rep(500,500), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_500_500)<-paste("set",1:500)
t_500g_500l<-system.time( res<-mitch_calc(k1,s_500_500,bootstraps=0,priority="significance",cores=8) )

s_2000_500<-sapply(rep(2000,500), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_2000_500)<-paste("set",1:500)
t_2000g_500l<-system.time( res<-mitch_calc(k1,s_2000_500,bootstraps=0,priority="significance",cores=8) )


# no genesets=5000 now vay the set size
s_20_5000<-sapply(rep(20,5000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_20_5000)<-paste("set",1:5000)
t_20g_5000l<-system.time( res<-mitch_calc(k1,s_20_5000,bootstraps=0,priority="significance",cores=8) )

s_100_5000<-sapply(rep(100,5000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_100_5000)<-paste("set",1:5000)
t_100g_5000l<-system.time( res<-mitch_calc(k1,s_100_5000,bootstraps=0,priority="significance",cores=8) )

s_500_5000<-sapply(rep(500,5000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_500_5000)<-paste("set",1:5000)
t_500g_5000l<-system.time( res<-mitch_calc(k1,s_500_5000,bootstraps=0,priority="significance",cores=8) )

s_2000_5000<-sapply(rep(2000,5000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_2000_5000)<-paste("set",1:5000)
t_2000g_5000l<-system.time( res<-mitch_calc(k1,s_2000_5000,bootstraps=0,priority="significance",cores=8) )


# no genesets = 20000 now vay the set size
s_20_20000<-sapply(rep(20,20000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_20_20000)<-paste("set",1:20000)
t_20g_20000l<-system.time( res<-mitch_calc(k1,s_20_20000,bootstraps=0,priority="significance",cores=8) )

s_100_20000<-sapply(rep(100,20000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_100_20000)<-paste("set",1:20000)
t_100g_20000l<-system.time( res<-mitch_calc(k1,s_100_20000,bootstraps=0,priority="significance",cores=8) )

s_500_20000<-sapply(rep(500,20000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_500_20000)<-paste("set",1:20000)
t_500g_20000l<-system.time( res<-mitch_calc(k1,s_500_20000,bootstraps=0,priority="significance",cores=8) )

s_2000_20000<-sapply(rep(2000,20000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_2000_20000)<-paste("set",1:20000)
t_2000g_20000l<-system.time( res<-mitch_calc(k1,s_2000_20000,bootstraps=0,priority="significance",cores=8) )

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

pdf("fig6c.pdf")
matplot(times_c_wide,pch=1,type = c("b"),xlab="gene set size",ylab="elapsed time (s)" , main="variation of number and size of gene sets", axes=F)
grid()
axis(2)
axis(side=1,at=1:nrow(times_c_wide),labels=rownames(times_c_wide))
legend("topleft", legend = c("20000","5000","500","100"), col=4:1, pch=1)
mtext("5 dimensions, 20000 genes, 8 cores, 0 bootstraps")
dev.off()

write.table(times_c_wide,files="fig6c.tsv")

times_c_wide
}

#################################################
# next, see vary the number and size of gene sets
#   genesets = 100 500 5000 20000
#   geneset size = 20 100 500 2000
#   cores = 8
#   bootstraps = 500 
#   genes = 20000
#   5 dimensions
#################################################

figd<-function(){
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
t_20g_100l<-system.time( res<-mitch_calc(k1,s_20_100,bootstraps=500,priority="significance",cores=8) )

s_100_100<-sapply(rep(100,100), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_100_100)<-paste("set",1:100)
t_100g_100l<-system.time( res<-mitch_calc(k1,s_100_100,bootstraps=500,priority="significance",cores=8) )

s_500_100<-sapply(rep(500,100), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_500_100)<-paste("set",1:100)
t_500g_100l<-system.time( res<-mitch_calc(k1,s_500_100,bootstraps=500,priority="significance",cores=8) )

s_2000_100<-sapply(rep(2000,100), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_2000_100)<-paste("set",1:100)
t_2000g_100l<-system.time( res<-mitch_calc(k1,s_2000_100,bootstraps=500,priority="significance",cores=8) )


# no genesets=500 now vay the set size
s_20_500<-sapply(rep(20,500), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_20_500)<-paste("set",1:500)
t_20g_500l<-system.time( res<-mitch_calc(k1,s_20_500,bootstraps=500,priority="significance",cores=8) )

s_100_500<-sapply(rep(100,500), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_100_500)<-paste("set",1:500)
t_100g_500l<-system.time( res<-mitch_calc(k1,s_100_500,bootstraps=500,priority="significance",cores=8) )

s_500_500<-sapply(rep(500,500), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_500_500)<-paste("set",1:500)
t_500g_500l<-system.time( res<-mitch_calc(k1,s_500_500,bootstraps=500,priority="significance",cores=8) )

s_2000_500<-sapply(rep(2000,500), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_2000_500)<-paste("set",1:500)
t_2000g_500l<-system.time( res<-mitch_calc(k1,s_2000_500,bootstraps=500,priority="significance",cores=8) )


# no genesets=5000 now vay the set size
s_20_5000<-sapply(rep(20,5000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_20_5000)<-paste("set",1:5000)
t_20g_5000l<-system.time( res<-mitch_calc(k1,s_20_5000,bootstraps=500,priority="significance",cores=8) )

s_100_5000<-sapply(rep(100,5000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_100_5000)<-paste("set",1:5000)
t_100g_5000l<-system.time( res<-mitch_calc(k1,s_100_5000,bootstraps=500,priority="significance",cores=8) )

s_500_5000<-sapply(rep(500,5000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_500_5000)<-paste("set",1:5000)
t_500g_5000l<-system.time( res<-mitch_calc(k1,s_500_5000,bootstraps=500,priority="significance",cores=8) )

s_2000_5000<-sapply(rep(2000,5000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_2000_5000)<-paste("set",1:5000)
t_2000g_5000l<-system.time( res<-mitch_calc(k1,s_2000_5000,bootstraps=500,priority="significance",cores=8) )


# no genesets = 20000 now vay the set size
s_20_20000<-sapply(rep(20,20000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_20_20000)<-paste("set",1:20000)
t_20g_20000l<-system.time( res<-mitch_calc(k1,s_20_20000,bootstraps=500,priority="significance",cores=8) )

s_100_20000<-sapply(rep(100,20000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_100_20000)<-paste("set",1:20000)
t_100g_20000l<-system.time( res<-mitch_calc(k1,s_100_20000,bootstraps=500,priority="significance",cores=8) )

s_500_20000<-sapply(rep(500,20000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_500_20000)<-paste("set",1:20000)
t_500g_20000l<-system.time( res<-mitch_calc(k1,s_500_20000,bootstraps=500,priority="significance",cores=8) )

s_2000_20000<-sapply(rep(2000,20000), function(x) {list(as.character(sample(genes,x))) } ) ; names(s_2000_20000)<-paste("set",1:20000)
t_2000g_20000l<-system.time( res<-mitch_calc(k1,s_2000_20000,bootstraps=500,priority="significance",cores=8) )

times_d<-as.data.frame(rbind(t_20g_100l,
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

times_d$nsets<-nsets
times_d$ngenes<-setsize

times_d_wide<-t(acast(times_d,nsets ~ setsize, value.var="elapsed"))

pdf("fig6d.pdf")
matplot(times_d_wide,pch=1,type = c("b"),xlab="gene set size",ylab="elapsed time (s)" , main="variation of number and size of gene sets", axes=F)
grid()
axis(2)
axis(side=1,at=1:nrow(times_d_wide),labels=rownames(times_d_wide))
legend("topleft", legend = c("20000","5000","500","100"), col=4:1, pch=1)
mtext("5 dimensions, 20000 genes, 8 cores, 500 bootstraps")
dev.off()

write.table(times_d_wide,files="fig6d.tsv")

times_d_wide
}


times_a_wide<-figa()
times_b_wide<-figb()
times_c_wide<-figc()
times_d_wide<-figd()

save.image("fig6.RData")
