### Visualize Correlatoins Using Heatmaps
## correl and overl are from nov-11.r
correl<-getCor(pvalues[1])
for(i in seq(2,n))
    correl<-abind(correl,getCor(pvalues[i]),along=3)

overl<-getOve(pvalues[1])
cats<-colnames(overl)
for(i in seq(2,n))
    overl<-abind(overl,getOve(pvalues[i],cats),along=3)
correl<-renameColumn(correl,"ecfc.tsa","ecfc-tsa")


oi<-addNames(outer(seq(22),seq(22),Vectorize(function(i,j) lm(y~x,data.frame(y=correl[i,j,],x=pvalues))$coef[2])),swapFun(categories),swapFun(categories))

ol<-addNames(outer(seq(22),seq(22),Vectorize(function(i,j) lm(y~x,data.frame(y=overl[i,j,],x=pvalues))$coef[2])),swapFun(categories),swapFun(categories))


png("~/Dropbox/pvalue-trend-scale.png")
y<-seq(range(oi)[1],range(oi)[2],length=12)
image(y,0,matrix(seq(12)),ylab="",xlab="")
dev.off()
png("~/Dropbox/overlap-pvalue-trend.png")
heatmapWrapper(ol,range(oi))
dev.off()
png("~/Dropbox/correlation-pvalue-trend.png")
heatmapWrapper(oi,range(oi))
dev.off()

### Coloured Dendrograms

x<-addColnames(read.table("~/Dropbox/Data/22x22-pvalue=0_heights.matrix",T),swapFun(categories))

png("~/Dropbox/22x22-dend-cap.png")
dend(x,norm=pass,n=3,linkage="average",colours=c("green","red","blue"))
dev.off()
png("~/Dropbox/22x22-dend-ccp.png")
dend(x,norm=pass,n=3,linkage="complete",colours=c("green","red","blue"))
dev.off()
png("~/Dropbox/22x22-dend-caq.png")
dend(x,n=5,linkage="average",colour=c("green","red","blue","brown","orange"))
dev.off()
png("~/Dropbox/22x22-dend-ccq.png")
dend(x,n=4,linkage="complete",colour=c("green","blue","red","orange"))
dev.off()
png("~/Dropbox/22x22-dend-map.png")
dend(x,n=2,norm=pass,dist=function(x)dist(t(x),"manh"),linkage="average",distanceName="Manhatten",colour=c("green","blue","orange","brown","red"))
dev.off()
png("~/Dropbox/22x22-dend-mcp.png")
dend(x,n=2,norm=pass,dist=function(x)dist(t(x),"manh"),linkage="complete",distanceName="Manhatten",colour=c("green","blue","orange","brown","red"))
dev.off()
png("~/Dropbox/22x22-dend-maq.png")
dend(x,n=5,norm=qn,dist=function(x)dist(t(x),"manh"),linkage="average",distanceName="Manhatten",colour=c("green","blue","orange","brown","red"))
dev.off()
png("~/Dropbox/22x22-dend-mcq.png")
dend(x,n=5,norm=qn,dist=function(x)dist(t(x),"manh"),linkage="complete",distanceName="Manhatten",colour=c("blue","green","red","brown","orange"))
dev.off()


#plot(hclust(dist(t(qn(x))),"average"),hang=-1)

#plot(hclust(as.dist(1-cor(x)),"average"),hang=-1)

### pvlaues and pcs
q<-lapply(c(3,5,7,9,11),function(i) stackedContribWrapper(mfn22(pvalues[i],TRUE),mfn22(pvalues[i]),categories,seq(4),swapFunB,swapFunC))


png("~/Dropbox/22x22-pc1-4-pvalue=5-20-r.png",height=1250,width=1000)
multiPlot(lapply(createZip(X2(5,4,FALSE)),function(i)q[[i[1]]][[i[2]]]),5,4)
dev.off()

p<-stackedContribWrapper(mfn22(pvalues[1],TRUE),mfn22(pvalues[1]),categories,c(1,3,5,7),swapFunB,swapFunC)

png("~/Dropbox/22x22-pc1,3,5,7-pvalue=5-r.png",height=(4*480),width=480)
multiPlot(p,4,1)
dev.off()

### Correlatoin Data For the following Sections

matr<-do.call(abind, append(lapply( pvalues[1:30],function(i){
    heights<-read.table(mfn22(i,TRUE),header=T)
    colnames(heights)<-categories
    prc<-pca(heights)
    t(prc$normData) %*% prc$eigenVectors
}), list(along=3)))


matrC<-do.call(abind, append(lapply(pvalues[1:30], function(i) pca2Matr(pca(addColnames(read.table(paste("~/thesis-november/cont-pvalue=",i,"_heights.matrix",sep=""),header = T),conts)))), list(along=3)))
                         


### Treatment
png("~/Dropbox/22x22-summary.png")
plots<-lapply(list(c(1,2),c(1,3),c(3,5),c(3,7)),function(i){
    plotPCMat(matr,c(pcsO[i[1]],pcsO[i[2]]),categories,swapFun,swapFunB,n=1,sdf=NULL,text=NULL,legend=NULL)})
multiPlot(plots,2,2)
dev.off()

### Permutate
png("~/Dropbox/22x22-pca-perm.png")
heights<-addColnames(read.table(mfn22(pvalues[1],height=TRUE),header=TRUE),categories)
rems<-list(c("k562_1","k562_2"),
           c("jurk_sandar","jurk","jurk_sandar_1"),           
           c("ecfc-tsa"),
           c("cd34","cd34_new"),
           c("rpmi_1","rpmi_2"),
           c("tall_p1"))
matrR<-lapply(rems,Compose(function(x)removeColumn(heights,x),pca,pca2Matr))
p<-lapply(seql(rems), function(i) plotPCMat2D(matrR[[i]],c(pcsO[1],pcsO[3]),removeColumn(categories,rems[[i]]),swapFun,swapFunB,swapFunC))
multiPlot(p[c(1,6,3,4)],2,2)
dev.off()

### control
png("~/Dropbox/22x22-pca-control.png")
plots<-lapply(list(c(1,2),c(1,3),c(3,5),c(3,7)),function(i){
    plotPCMat(matrC,c(pcsO[i[1]],pcsO[i[2]]),conts,Compose(mockToNorm,swapFun),Compose(mockToNorm,swapFunB),n=1,sdf=NULL,text=NULL,legend=NULL)})
multiPlot(plots,2,2)
dev.off()

### Pvalue variations

#png("~/Dropbox/22x22-pca-pvalue=5-20.png",width=480,height=480*5/2)
plots<-lapply(list(c(1,2),c(1,3),c(1,4),c(3,5),c(3,7)),function(i){
    plotPCMat(matr,c(pcsO[i[1]],pcsO[i[2]]),categories,swapFun,swapFunB,n=1,sdf=NULL,text=NULL,legend=NULL)})
plots<-c(rbind(plots,lapply(list(c(1,2),c(1,3),c(1,4),c(3,5),c(3,7)),function(i){
    plotPCMat(matr,c(pcsO[i[1]],pcsO[i[2]]),categories,swapFun,swapFunB,n=9,sdf=NULL,text=NULL,legend=NULL)})))
multiPlot(plots,5,2)
#dev.off()


### Recrate down mean square distance
png("~/Dropbox/22x22-stability-1-4.png")
par(mfrow=c(2,2),mar=c(2,2,2,1),oma = c(4, 4, 0, 0))
for(i in seq(4)){
t<-mapply(function(i,j,pc) MSD(matr[,pc,i],matr[,pc,j]),3:29,4:30,MoreArgs = list(i))
plot(pvalues[3:29],log(t,10),main=paste("PC",i,sep=""),cex=1.5,cex.main=2,cex.axis=1.5,pch=19,xlab="",ylab="",ylim=c(2,5.3))
}
mtext('P-Values', side = 1, outer = TRUE, line = 2,cex=2)
mtext('Log Mean Squared Difference', side = 2, outer = TRUE, line = 2,cex=2)
dev.off()



### 22x22 AFS showing 3 points
nafs<-do.call(rbind,lapply( pvalues,function(i){  
    bedData<-readTagMatrix(mfn22(i),categories)
    dim(bedData)
   }
))

png("~/Dropbox/22x22-afs-size.png")
par(mar=c(5,5,5,2))
plot(NULL,xlim=range(pvalues[3:29]),ylim=range(nafs[3:29]/nafs[3]),ylab="Peak Proportions",xlab="P-Values",main="AFS Peak Proportions",cex=1.5,type="l",lwd=4,cex.axis=1.5,cex.main=2,cex.lab =2 )
lines(pvalues[3:29],(nafs[3:29]/nafs[3]),lwd=4)
lines(x=c(10,10),y=c(-100,nafs[which(pvalues==10)]/nafs[3]),col="red",lwd=4)
lines(x=c(-100,10),y=rep(nafs[which(pvalues==10)]/nafs[3],2),col="red",lwd=4)
points(x=10,y=nafs[which(pvalues==10)]/nafs[3],col="red",cex=1.5,pch=19)
lines(x=c(20,20),y=c(-100,nafs[which(pvalues==20)]/nafs[3]),col="blue",lwd=4)
lines(x=c(-100,20),y=rep(nafs[which(pvalues==20)]/nafs[3],2),col="blue",lwd=4)
points(x=20,y=nafs[which(pvalues==20)]/nafs[3],col="blue",cex=1.5,pch=19)
lines(x=c(30,30),y=c(-100,nafs[which(pvalues==30)]/nafs[3]),col="green",lwd=4)
lines(x=c(-100,30),y=rep(nafs[which(pvalues==30)]/nafs[3],2),col="green",lwd=4)
points(x=30,y=nafs[which(pvalues==30)]/nafs[3],col="green",cex=1.5,pch=19)
text(c(17,27,37),sapply(c(10,20,30),function(x)nafs[which(pvalues==x)]/nafs[3])+c(.025,.025,.025), labels=c("(10 , 0.71)","(20 , 0.46)","(30 , 0.36)"),cex=1.5)
dev.off()


### PCA pvalue=0 (9) 1,3
sdf<-shiftdf("ecfc.tsa",0,-0.05)
sdf<-shiftdf("cd34",-0.12,0,sdf)
sdf<-shiftdf("cd133",-0.12,0,sdf)
sdf<-shiftdf("cd34_new",0,-.05,sdf)
sdf<-shiftdf("meka",0,-.05,sdf)
sdf<-shiftdf("eryt_f",0,.05,sdf)
sdf<-shiftdf("eryt",0,.05,sdf)
sdf<-shiftdf("eryt_a",-.17,0,sdf)
sdf<-shiftdf("k562_1",-.15,0,sdf)
sdf<-shiftdf("k562_2",-.05,-.05,sdf)
sdf<-shiftdf("jurk_sandar_1",0,.05,sdf)
sdf<-shiftdf("rpmi_1",.15,0.04,sdf)
sdf<-shiftdf("rpmi_2",.15,0,sdf)
sdf<-shiftdf("tall_p1",.01,.07,sdf)
sdf<-shiftdf("tall_p2_1",.15,.02,sdf)
sdf<-shiftdf("tall_p2_2",.2,0,sdf)
sdf<-shiftdf("jurk",-.02,.04,sdf)
sdf<-shiftdf("tall_p3_1",-0.1,-0.06,sdf)
sdf<-shiftdf("tall_p3_2",0.15,-0.06,sdf)
sdf<-shiftdf("cem_1",-0.11,.05,sdf)
sdf<-shiftdf("cem_2",.15,-.02,sdf)
sdf<-shiftdf("jurk_sandar",-.03,-.04,sdf)
p<-plotPCMat2D(matr[,,1] ,c(pcsO[1],pcsO[3]), categories=categories, sdf=sdf,swapFun=swapFun,swapFunB=swapFunB,swapFunC=swapFunC,text=TRUE)
ggsave(p,filename="~/Dropbox/22x22-pca-1,3-pvalue=0.png")

### PCA pvalue=20 (9) 1,2 1,4
sdf<-shiftdf("ecfc.tsa",0,+0.05)
sdf<-shiftdf("cd34",-0.12,0,sdf)
sdf<-shiftdf("cd133",-0.12,0,sdf)
sdf<-shiftdf("cd34_new",0,-.05,sdf)
sdf<-shiftdf("meka",0,-.05,sdf)
sdf<-shiftdf("eryt_f",0,-.05,sdf)
sdf<-shiftdf("eryt",-.16,-.02,sdf)
sdf<-shiftdf("eryt_a",-.16,.03,sdf)
sdf<-shiftdf("k562_1",.02,-.05,sdf)
sdf<-shiftdf("k562_2",.01,.05,sdf)
sdf<-shiftdf("jurk_sandar_1",.2,0,sdf)
sdf<-shiftdf("rpmi_1",0.1,-.07,sdf)
sdf<-shiftdf("rpmi_2",-.07,-0.08,sdf)
sdf<-shiftdf("tall_p1",0.01,.05,sdf)
sdf<-shiftdf("tall_p2_1",0,-.07,sdf)
sdf<-shiftdf("tall_p2_2",.15,-0.04,sdf)
sdf<-shiftdf("jurk",.03,0.05,sdf)
sdf<-shiftdf("tall_p3_1",0.15,0,sdf)
sdf<-shiftdf("tall_p3_2",0.15,0,sdf)
sdf<-shiftdf("cem_1",-0.15,0,sdf)
sdf<-shiftdf("cem_2",0.15,0,sdf)
sdf<-shiftdf("jurk_sandar",-.03,0.05,sdf)
p<-plotPCMat2D(matr[,,9] ,c(pcsO[1],pcsO[2]), categories=categories, sdf=sdf,swapFun=swapFun,swapFunB=swapFunB,swapFunC=swapFunC,text=TRUE)
ggsave(p,file="~/Dropbox/22x22-pca-1,2-pvalue=20.png")

sdf<-shiftdf("ecfc.tsa",0,+0.05)
sdf<-shiftdf("cd34",-0.12,0,sdf)
sdf<-shiftdf("cd133",-0.12,0,sdf)
sdf<-shiftdf("cd34_new",0,-.05,sdf)
sdf<-shiftdf("meka",0,-.05,sdf)
sdf<-shiftdf("eryt_f",-0.16,-.01,sdf)
sdf<-shiftdf("eryt",-.16,0,sdf)
sdf<-shiftdf("eryt_a",-.16,0.01,sdf)
sdf<-shiftdf("k562_1",.02,-.05,sdf)
sdf<-shiftdf("k562_2",.01,.05,sdf)
sdf<-shiftdf("jurk_sandar_1",.12,0,sdf)
sdf<-shiftdf("rpmi_1",0,-.05,sdf)
sdf<-shiftdf("rpmi_2",0,+0.05,sdf)
sdf<-shiftdf("tall_p1",0.01,.05,sdf)
sdf<-shiftdf("tall_p2_1",0.14,-.02,sdf)
sdf<-shiftdf("tall_p2_2",0.16,.02,sdf)
sdf<-shiftdf("jurk",.04,0.06,sdf)
sdf<-shiftdf("tall_p3_1",0.15,.03,sdf)
sdf<-shiftdf("tall_p3_2",0.15,-.03,sdf)
sdf<-shiftdf("cem_1",.05,-.05,sdf)
sdf<-shiftdf("cem_2",0.17,-.02,sdf)
sdf<-shiftdf("jurk_sandar",-.03,-0.05,sdf)
p<-plotPCMat2D(matr[,,9] ,c(pcsO[1],pcsO[4]), categories=categories, sdf=sdf,swapFun=swapFun,swapFunB=swapFunB,swapFunC=swapFunC,text=TRUE)
ggsave(p,file="~/Dropbox/22x22-pca1,4-pvalue=20.png")

ggsave(q[[1]][[1]],file="~/Dropbox/22x22-contrib-1-pvalue=5-nl.png")

ggsave(q[[1]][[1]],file="~/Dropbox/test.png")

ggsave(q[[1]][[3]],file="~/Dropbox/22x22-contrib-3-pvalue=5-nl.png")
ggsave(q[[4]][[1]],file="~/Dropbox/22x22-contrib-1-pvalue=20-nl.png")
ggsave(q[[4]][[2]],file="~/Dropbox/22x22-contrib-2-pvalue=20-nl.png")
ggsave(q[[4]][[4]],file="~/Dropbox/22x22-contrib-4-pvalue=20-nl.png")

p2<-stackedContribWrapper(mfn22(pvalues[1],TRUE),mfn22(pvalues[1]),categories,c(1,3,5),swapFunB,swapFunC,sum=FALSE)

ggsave(p2[[1]],filename="~/Dropbox/22x22-contrib-1-ns-pvalue=5-nl.png")
ggsave(p2[[2]],filename="~/Dropbox/22x22-contrib-3-ns-pvalue=5-nl.png")
ggsave(p2[[3]],filename="~/Dropbox/22x22-contrib-5-ns-pvalue=5-nl.png")

source("~/r-workspace/nov-functions.r")

q[[1]]
  

p[1]

q[1]
