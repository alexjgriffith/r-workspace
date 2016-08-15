library(CCCA)
library(functional)
library(parallel)
library(permute)
library(abind)


cats<-readCategories("~/Dropbox/UTX-Alex/jan/catagories")
test<-read.table("~/thesis-november/22x22-pvalue=0_heights.matrix",header=T)
colnames(test)<-cats
covBase<-cor(test)

y<-dim(test)[1]
x<-dim(test)[2]

ind<-repR(shuffle,10,len)

n=1
n2=2
cs<-makeForkCluster(n,renice=0)

permuteMatrix<-function(matrix)
    apply(matrix,2,function(x)x[shuffle(length(x))])

ret<-parallel::parLapply(cs,seq(n2),function(p,mat) cov(permuteMatrix,), shuffleMatrix,mat,x,y)

stopCluster(cs)

shuffleMatrix<-function(x,y)
    t(repR(shuffle,x,y))

function(x) head())


sm<-shuffleMatrix(x,y)

t<-outer(sm,seq(22),Vectorize(function(x,y,test){test[x[,y],y]}) ,test)


mat<-matrix(seq(100),ncol=10)

subset<-shuffleMatrix(10,10)
