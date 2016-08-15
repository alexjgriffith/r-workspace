source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/jan-functions.r")
source("~/r-workspace/feb-functions.r")
source("~/r-workspace/feb-variables.r")


env=getPRC20(2)

#png("~/Dropbox/Data/new-figs/ecfc-selection-hist-0.5.png")
width=0.5
xlim=c(-2,8)
v="ECFC"
par(mfrow=c(2,2),mar=c(2,2,3,2),oma=c(0, 0, 2, 0))
x<-log(addColnames(do.call(cbind,lapply(unique(swapFunD(env$categories)),function(x) apply(cbind(env$prc$normData[,x== swapFunD(env$categories)],0),1,sum))),unique(swapFunD(env$categories)))[,v])
breaks<-diff(range(Filter(function(x)x!=-Inf,x)))/width
hist(x,breaks,xlim=xlim,main="No Filter",ylab="",xlab="",col="green")
##
x<-log(addColnames(do.call(cbind,lapply(unique(swapFunD(env$categories)),function(x) apply(cbind(env$prc$normData[,x== swapFunD(env$categories)],0),1,sum))),unique(swapFunD(env$categories)))[env$ECFCB,v])
breaks<-diff(range(Filter(function(x)x!=-Inf,x)))/width
hist(x,breaks,xlim=xlim,main="PC2",ylab="",xlab="",col="green")
##
x<-log(addColnames(do.call(cbind,lapply(unique(swapFunD(env$categories)),function(x) apply(cbind(env$prc$normData[,x== swapFunD(env$categories)],0),1,sum))),unique(swapFunD(env$categories)))[env$ECFCA,v])
breaks<-diff(range(Filter(function(x)x!=-Inf,x)))/width
hist(x,breaks,xlim=xlim,main="PC4",ylab="",xlab="",col="green")
##
x<-log(addColnames(do.call(cbind,lapply(unique(swapFunD(env$categories)),function(x) apply(cbind(env$prc$normData[,x== swapFunD(env$categories)],0),1,sum))),unique(swapFunD(env$categories)))[env$ECFC,v])
breaks<-diff(range(Filter(function(x)x!=-Inf,x)))/width
hist(x,breaks,xlim=xlim,main="PC1 & PC2",ylab="",xlab="",col="green")
title("Selecting ECFC Specifying Peaks",outer=TRUE)
xlab("11")
#dev.off()


png("~/Dropbox/Data/new-figs/selection-hist-0.5.png")
width=0.25
xlim=c(3,8)
par(mfrow=c(2,2),mar=c(2,2,3,2),oma=c(0, 0, 2, 0))
x<-log(addColnames(do.call(cbind,lapply(unique(swapFunD(env$categories)),function(x) apply(cbind(env$prc$normData[,x== swapFunD(env$categories)],0),1,sum))),unique(swapFunD(env$categories)))[env$ECFC,"ECFC"])
breaks<-diff(range(Filter(function(x)x!=-Inf,x)))/width
hist(x,breaks,xlim=xlim,main="ECFC",ylab="",xlab="",col="green")
##
x<-log(addColnames(do.call(cbind,lapply(unique(swapFunD(env$categories)),function(x) apply(cbind(env$prc$normData[,x== swapFunD(env$categories)],0),1,sum))),unique(swapFunD(env$categories)))[env$Erythroid,"Erythroid"])
breaks<-diff(range(Filter(function(x)x!=-Inf,x)))/width
hist(x,breaks,xlim=xlim,main="Erythroid",ylab="",xlab="",col="red")
##
x<-log(addColnames(do.call(cbind,lapply(unique(swapFunD(env$categories)),function(x) apply(cbind(env$prc$normData[,x== swapFunD(env$categories)],0),1,sum))),unique(swapFunD(env$categories)))[env$HSC,"HSC"])
breaks<-diff(range(Filter(function(x)x!=-Inf,x)))/width
hist(x,breaks,xlim=xlim,main="HSC",ylab="",xlab="",col="orange")
##
x<-log(addColnames(do.call(cbind,lapply(unique(swapFunD(env$categories)),function(x) apply(cbind(env$prc$normData[,x== swapFunD(env$categories)],0),1,sum))),unique(swapFunD(env$categories)))[env$Leukemia,"Leukemia"])
breaks<-diff(range(Filter(function(x)x!=-Inf,x)))/width
hist(x,breaks,xlim=xlim,main="Leukemia",ylab="",xlab="",col="blue")
##
dev.off()



env=getPRC20(0.25)

width=0.25
xlim=c(0,8)
par(mfrow=c(2,2),mar=c(2,2,3,2),oma=c(0, 0, 2, 0))
x<-log(addColnames(do.call(cbind,lapply(unique(swapFunD(env$categories)),function(x) apply(cbind(env$prc$normData[,x== swapFunD(env$categories)],0),1,sum))),unique(swapFunD(env$categories)))[env$ECFC,"ECFC"])
breaks<-diff(range(Filter(function(x)x!=-Inf,x)))/width
hist(x,breaks,xlim=xlim,main="ECFC",ylab="",xlab="",col="green")
##
x<-log(addColnames(do.call(cbind,lapply(unique(swapFunD(env$categories)),function(x) apply(cbind(env$prc$normData[,x== swapFunD(env$categories)],0),1,sum))),unique(swapFunD(env$categories)))[env$Erythroid,"Erythroid"])
breaks<-diff(range(Filter(function(x)x!=-Inf,x)))/width
hist(x,breaks,xlim=xlim,main="Erythroid",ylab="",xlab="",col="red")
##
x<-log(addColnames(do.call(cbind,lapply(unique(swapFunD(env$categories)),function(x) apply(cbind(env$prc$normData[,x== swapFunD(env$categories)],0),1,sum))),unique(swapFunD(env$categories)))[env$HSC,"HSC"])
breaks<-diff(range(Filter(function(x)x!=-Inf,x)))/width
hist(x,breaks,xlim=xlim,main="HSC",ylab="",xlab="",col="orange")
##
x<-log(addColnames(do.call(cbind,lapply(unique(swapFunD(env$categories)),function(x) apply(cbind(env$prc$normData[,x== swapFunD(env$categories)],0),1,sum))),unique(swapFunD(env$categories)))[env$Leukemia,"Leukemia"])
breaks<-diff(range(Filter(function(x)x!=-Inf,x)))/width
hist(x,breaks,xlim=xlim,main="Leukemia",ylab="",xlab="",col="blue")


hist()




source("~/Dropbox/thesis/R/contrib_new.R")

penv<-function(env,pc)with(env,{    
    stackedContrib(prc$eigenVectors[,pc],"contrib2",mergeFun(over[4:dim(over)[2]],swapFunD),swapFun=swapFunD,colourOveride =swapFunC,n=6)    
})

penv(env,2)

env<-getPRC20(0.5)

