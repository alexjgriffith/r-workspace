source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/jan-functions.r")
source("~/r-workspace/feb-functions.r")
source("~/r-workspace/feb-variables.r")

source("~/Dropbox/thesis/R/contrib_new.R")

env<-getPRC(20)

penv<-function(env)with(env,{    
    stackedContrib(prc$eigenVectors[,4],"contrib2",mergeFun(over[4:dim(over)[2]],swapFunD),swapFun=swapFunD,colourOveride =swapFunC,n=6)    
})

penv<-function(env)with(env,{    
    plotPCMat2D(pca2Matr(prc),c(1,2),categories,swapFun,swapFunD,swapFunC) 
})


mergeFun(with(env,{over[normalize(prc$eigenVectors[,1])>1,4:dim(over)[2]]}),swapFunD)

providing<-function(env,name,pc,fun=">",x=seq(6),swapFun=swapFunD){
    sapply(x,function(sd){
        x<-mergeFun(with(env,{over[subsetPRC(prc,pc,sd,fun),4:dim(over)[2]]}),swapFun)
        #x<-mergeFun(with(env,{over[do.call(fun,list(normalize(prc$eigenVectors[,pc]),sd)),4:dim(over)[2]]}),swapFunD)
        length(which(x[, name]==1))/dim(x)[1]})}


p2l<-function(env,name,pc,fun=">",x=seq(6),swapFun=swapFunD){
    r=providing(env,name,pc,fun,x=x)
    plot(r,ylab="Proportion",xlab="SDs",main=name,ylim=c(0,1))
    z<-sapply(x,function(sd)length(which(with(env,{subsetPRC(prc,pc,sd,fun)}))))
    points(z/z[1],col="red")
    list(z/z[1],r)
}



par(mfrow=c(2,2))
p2l(env,"HSC",4,">",seq(6))
p2l(env,"ECFC",4,"<",-seq(6))
p2l(env,"Erythroid",1,">",seq(6))
p2l(env,"Leukemia",1,"<",-seq(6))



with(getPRC(20,"single"),{
    stem=normalize(prc$eigenVectors[,4])<(-2)
    x<-mergeFun(over[which(apply(cbind(heights[stem,swapFunD(categories)=="ECFC"],0),1,sum)<20),4:dim(over)[2]],swapFunD)
    apply(x,2,mean)
})


hist(with(getPRC(20,"single"),{
    stem=normalize(prc$eigenVectors[,4])<(-3)
    x<-apply(cbind(heights[stem,swapFunD(categories)=="ECFC"],0),1,sum)
    #l<-floor(length(x))
    #sort(x)[floor(l*0.125):floor(l-l*0.125)]
    x
}),breaks=500,xlim=c(0,700),xlab="",ylab="",main="")


hist(with(env,{
    stem=normalize(prc$eigenVectors[,4])<(-2)
    #x<-apply(cbind(heights[stem,swapFunD(categories)=="ECFC"],0),1,sum)
    x<-apply(cbind(qn(heights)[stem,swapFunD(categories)=="ECFC"],0),1,sum)
    x
}))


with(getPRC(20,"combined",forward="22"),{
    n=1
    par(mfrow=c(2,2))
    stem=normalize(prc$eigenVectors[,1])>(n)
    x<-apply(cbind(qn(heights)[stem,swapFunD(categories)!="Erythroid"],0),1,max)
    hist(x,xlab="",ylab="",main="Erythroid",xlim=c(0,200),breaks=500)
    stem=normalize(prc$eigenVectors[,1])<(-n)
    x<-apply(cbind(qn(heights)[stem,swapFunD(categories)=="Leukemia"],0),1,max)
    hist(x,xlab="",ylab="",main="Leukemia",xlim=c(0,200),breaks=500)
    stem=normalize(prc$eigenVectors[,4])>(n)
    x<-apply(cbind(qn(heights)[stem,swapFunD(categories)=="HSC"],0),1,max)
    hist(x,xlab="",ylab="",main="HSC",xlim=c(0,200),breaks=500)
    stem=normalize(prc$eigenVectors[,4])<(-(n^2/2))
    stem2=normalize(prc$eigenVectors[,2])<(-(n^2/2))
    x<-apply(cbind(qn(heights)[stem&stem2,swapFunD(categories)=="ECFC"],0),1,max)
    hist(x,xlab="",ylab="",main="ECFC",xlim=c(0,200),breaks=500)     
})


penv(getPRC(20))


with(env,{
eryt=bed[normalize(prc$eigenVectors[,1])>2 & bed$chr!="chrM",]
jurk=bed[normalize(prc$eigenVectors[,1])<(-2) & bed$chr!="chrM",]
ecfc=bed[normalize(prc$eigenVectors[,4])<(-1.3) &normalize(prc$eigenVectors[,2])<(-1.3) & bed$chr!="chrM",]
stem=bed[normalize(prc$eigenVectors[,4])>2 & bed$chr!="chrM",]
list(dim(eryt)[1],dim(jurk)[1],dim(ecfc)[1],dim(stem)[1])
#printBB(eryt,"~/Dropbox/eryt_sd_2_new.bb")
#printBB(jurk,"~/Dropbox/jurk_sd_2_new.bb")
printBB(ecfc,"~/Dropbox/ecfc_sd_2_new.bb")
#printBB(stem,"~/Dropbox/hsc_sd_2_new.bb")
tts<-function(x){
    paste0(x[,1],":",x[,2]-2000,"-",x[,3]+2000)
}
write.table(tts(eryt),"~/Dropbox/UTX-Alex/br-data/eryt_sd_2.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(tts(jurk),"~/Dropbox/UTX-Alex/br-data/jurk_sd_2.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(tts(ecfc),"~/Dropbox/UTX-Alex/br-data/ecfc_sd_2.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(tts(stem),"~/Dropbox/UTX-Alex/br-data/stem_sd_2.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
})



plot(providing(env,"HSC",4,">",seq(6)),ylab="Proportion",main="HSC")
plot(providing(env,"ECFC",4,"<",-seq(6)),ylab="Proportion",main="ECFC")
plot(providing(env,"Erythroid",1,">",seq(6)),ylab="Proportion",main="Erythroid")
plot(providing(env,"Leukemia",1,"<",-seq(6)),ylab="Proportion",main="Leukemia")

plot(sapply(seq(6),function(sd){
x<-mergeFun(with(env,{over[subsetPRC(prc,1,sd),4:dim(over)[2]]}),swapFunD)
length(which(x[,"Erythroid"]==1))/dim(x)[1]}),ylab="Proportion",main="Erythroid")
plot(sapply(seq(6),function(sd){
x<-mergeFun(with(env,{over[subsetPRC(prc,1,-1*sd,"<"),4:dim(over)[2]]}),swapFunD)
length(which(x[,"Leukemia"]==1))/dim(x)[1]}),ylab="Proportion",main="Leukemia")
plot(sapply(seq(6),function(sd){
x<-mergeFun(with(env,{over[subsetPRC(prc,4,-1*sd,"<"),4:dim(over)[2]]}),swapFunD)
length(which(x[,"ECFC"]==1))/dim(x)[1]}),ylab="Proportion",main="ECFC")



plotPCMat2D(pca2Matr(prc),c(1,2),categories,swapFun=swapFunM,swapFunB=swapFunB,swapFunC=swapFunC)

plotPCs(prc$eigenVectors,c(1,6),prc$normData,categories)

plot(hclust(as.dist(1-cor(heights))),hang=-1)



env=getPRC(20,forward=23)
plotPCs(env$prc$eigenVectors,c(1,4),env$prc$normData,swapFun(env$categories))
penv(env)

env<-getPRC(20,"combined",forward=22)
#penv(env)

sort(with(env,{heights[normalize(prc$eigenVectors[,4])<(-1.3)&normalize(prc$eigenVectors[,2])<(-1.3),"ecfc.tsa"]}))

env<-getPRC20(2)

