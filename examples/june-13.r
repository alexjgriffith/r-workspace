str(array(a,dim=c(sapply(seq(3),rep,10))))

library(scatterplot3d)
a<-rnorm(100,-1,1)
#d3p<-array(a,dim=c(10,10,10))
d2p<-array(a,dim=c(10,10))



scatterplot3d(d2p[,1],d2p[,2])

library(CCCA)

source("~/r-workspace/project.r")

env<-getPRC20(2)


png("~/Dropbox/UTX-Alex/Paper/3d-scatterplot.png")
scatterplot3d(pca2Matr(env$prc)[,c(2,4,1)],pch=16,color=swapFunC(swapFunD(env$categories)),type="h",angle=45,box=FALSE,x.ticklabs="",y.ticklabs="",z.ticklabs="")
dev.off()

png("~/Dropbox/UTX-Alex/Paper/3d-scatterplot-alt.png")
scatterplot3d(pca2Matr(env$prc)[,c(2,4,1)],pch=16,color=swapFunC(swapFunD(env$categories)),type="p",angle=45,box=FALSE,x.ticklabs="",y.ticklabs="",z.ticklabs="")
dev.off()

n<-prcomp(env$prc$normData)

(n$sdev)^2/sum(n$sdev^2)

# 40.4%
# 17.5%
# 6.71%
# 5.78%
# 3.85%


# pca2matr is from CCC.r in r-workspace

library(ggplot2)

plotPCMat2D(pca2Matr(env$prc),c(1,3),env$categories,swapFun,swapFunD,swapFunC)

weights<-pca2Matr(env$prc)

png("~/Dropbox/UTX-Alex/Paper/PC1-PC3.png")
plotPCs(env$prc$eigenVectors,c("PC1","PC3"),env$prc$normData,env$categories)
dev.off()
