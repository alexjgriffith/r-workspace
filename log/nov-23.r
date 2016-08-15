# make file name
mfnMake<-function(fun)
    function(pvalue=0,heights=FALSE){    
        if(heights)
            heights="_heights"
        else
            heights=""
        fun(pvalue,heights)
}

mfn4<-mfnMake(function(pvalue,heights)
    paste("~/Dropbox/Data/matrix/4x4-pvalue=",pvalue,heights,".matrix",sep=""))

mfn22<-mfnMake(function(pvalue,heights)
    paste("~/Dropbox/Data/22x22-pvalue=",pvalue,heights,".matrix",sep=""))

mfn


categories<-colnames(read.table(mfn22(0),header=T)[,4:(22+3)])

bedData<-read.table(mfn22(0),header=T)[,1:3]

heights<-read.table(mfn22(0,TRUE),header=T)

colnames(heights)<-categories

pvalues<-seq(from=0,by=2.5,length=30)



npeaks<-do.call(rbind,lapply( pvalues,function(i){
    bedData<-read.table(mfn4(i),header=T)[,4:7]
    colnames(bedData)<-swapFun(colnames(bedData))
    apply(bedData,2,sum)
   }
))

nafs<-do.call(rbind,lapply( pvalues,function(i){
    bedData<-read.table(mfn4(i),header=T)[,4:7]
    dim(bedData)
   }
))

nafs<-do.call(rbind,lapply( pvalues,function(i){  
    bedData<-readTagMatrix(mfn22(i),categories)
    dim(bedData)
   }
))

png("~/Dropbox/22x22-AFS-pvalue.png")
plot(pvalues[4:30],nafs[4:30,1],type="l",xlab="MACS Score", ylab="Peaks in AFS",lwd =3)
dev.off()

plotV<-function(x,y,ylab,title,sf=1){
    df<-data.frame(pvalues=x,peaks=y)
    ggplot(df,aes(x=pvalues,y=peaks))+geom_line(size=2*sf)+xlab("Cut Off (10^(-pvalue))")+ylab(ylab)+ggtitle(title)+
        theme_bw() +
            theme(axis.line = element_line(colour = "black"),
                  axis.text=element_text(size=12*sf),
        axis.title=element_text(size=14*sf,face="bold"),
        axis.line.x=element_line(size=1*sf),
        axis.line.y=element_line(size=1*sf),
        axis.ticks.y=element_line(size=1*sf),
        axis.ticks.x=element_line(size=1*sf),
        title=element_text(size=14*sf,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

}

pvWrap <-function(n)
    plotV(pvalues[3:30],npeaks[3:30,n]/npeaks[3,n],"Proportion of Peaks",colnames(npeaks)[n])

p1<-pvWrap(1)
p2<-pvWrap(2)
p3<-pvWrap(3)
p4<-pvWrap(4)
library(grid)

png("~/Dropbox/4x4-peaks-pvalue.png")
pushViewport(viewport(layout=grid.layout(2,2)))    
print(p1,vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(p2,vp=viewport(layout.pos.row=1,layout.pos.col=2))
print(p3,vp=viewport(layout.pos.row=2,layout.pos.col=1))
print(p4,vp=viewport(layout.pos.row=2,layout.pos.col=2))
dev.off()


png("~/Dropbox/4x4-afs-size.png")
plotV(pvalues[3:30],nafs[3:30,1]/nafs[3,1],"Proportion of Peaks","AFS",sf=2)
dev.off()

## get overlap info for venn diagram
## functions from nov-11



## bunch of histograms
dend(heights[,categoriesB],n=1,colours=c("red","orange","blue"),swapFun=swapFun)

p2<-dend(heights[,categoriesB],"euclidian",n=1,colours=c("red","orange","blue"),swapFun=swapFun)


sf=4
png("~/Dropbox/4x4-hclust-afs-combined-new-big.png",height=480*sf,width=480*sf)
par(mfrow=c(3,4),mar=c(0.5,1,3,1),cex=sf,lwd=sf)
data<-applySwapFun(heights[,c(categoriesB)],swapFun)
dlist<-list("Dissim"=function(x) as.dist(1-cor(x)),"Manhat"=function(x) dist(t(x),"man"),"Euclid"=function(x)dist(t(x),"euc"),"Pnorm"=function(x)dist(t(x),"min",100))
for(linkage in c("average","complete","single"))
for(distance in c("Dissim","Manhat","Euclid","Pnorm"))
    {
        plot(hclust(dlist[[distance]](data),linkage),hang=-1,sub="",xlab="",ylab=NULL,main=paste(distance,linkage),yaxt='n')
    }
dev.off()

        #dend(heights[,categoriesB],distance,linkage=linkage,n=1,swapFun=

pushViewport(viewport(layout=grid.layout(2,2)))    
print(p1,vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(p2,vp=viewport(layout.pos.row=1,layout.pos.col=2))



########## overlap stuff
dataSets<-lapply(categoriesB,function(i)
    read.table(paste("~/Dropbox/Data/august_peaks/",i,"~combined_mock_peaks.xls",sep=""),header=T))

afs4<-read.table("~/Dropbox/Data/matrix/4x4-pvalue=0.matrix",header=T)



d<-lapply(seq(4),function(x) afs4[which(apply(afs4[,4:7],1,sum)==x),])

a<-lapply(d,function(r) unlist(apply(r,1,function(x) lapply(dataSets[which(as.logical(as.integer(x[4:7])))],findClosed,x))))


a2<-unlist(apply(d[[2]],1,function(x) max(unlist(lapply(dataSets[which(as.logical(as.integer(x[4:7])))],findClosed,x)))))
a3<-unlist(apply(d[[3]],1,function(x) max(unlist(lapply(dataSets[which(as.logical(as.integer(x[4:7])))],findClosed,x)))))
a4<-unlist(apply(d[[4]],1,function(x) max(unlist(lapply(dataSets[which(as.logical(as.integer(x[4:7])))],findClosed,x)))))

ag1<-c(a2,a3,a4)

png("~/Dropbox/4-pv-hist.png")
par(mfrow=c(2,2))
hist(log(a[[1]]),breaks=seq(1,6,by=0.5),xlab="",ylab="",main="Unique")
hist(log(a[[2]]),breaks=seq(1,6,by=0.5),xlab="",ylab="",main="One Overlap")
hist(log(a[[3]]),breaks=seq(1,6,by=0.5),xlab="",ylab="",main="Two Overlaps")
hist(log(a[[4]]),breaks=seq(1,6,by=0.5),xlab="",ylab="",main="Three Overlaps")
dev.off()


#png("~/Dropbox/4x4-hist-pvalue.png")
par(mfrow=c(2,2))
for(i in seq(4))
    hist(dataSets[[i]]$X.log10.pvalue.,main=swapFun(categoriesB[i]),xlab="",ylab="")
#dev.off()

findClosed<-function(dataSet,peak){
    ss<-dataSet[as.character(dataSet$chr)==as.character(peak[1]),]
    d<-abs(ss$abs_summit-(as.numeric(peak[2])+as.numeric(peak[3]))/2)
    wi<-which.min(d)    
    ss[wi,"X.log10.pvalue."]#,d[wi])
}


# Find overlaps to make up to date ven diagram

# Make 22x22 plot of the first second third and fourth principle components at 10-5


library(abind)

categories<-colnames(read.table(mfn22(0),header=T)[,4:(22+3)])

bedData<-read.table(mfn22(0),header=T)[,1:3]

heights<-read.table(mfn22(0,TRUE),header=T)
prc<-pca(heights)



matr<-do.call(abind, append(lapply( pvalues[1:30],function(i){
    heights<-read.table(mfn22(i,TRUE),header=T)
    colnames(heights)<-categories
    prc<-pca(heights)
    t(prc$normData) %*% prc$eigenVectors
}), list(along=3)))




plotPCs(prc$eigenVectors,c(1,2),prc$normData,cats=swapFun(categories))

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
sdf<-shiftdf("jurk_sandar_1",.1,0.05,sdf)
sdf<-shiftdf("rpmi_1",.1,0.05,sdf)
sdf<-shiftdf("rpmi_2",.15,.02,sdf)
sdf<-shiftdf("tall_p1",.01,-.03,sdf)
sdf<-shiftdf("tall_p2_1",.15,.02,sdf)
sdf<-shiftdf("tall_p2_2",.2,0,sdf)
sdf<-shiftdf("jurk",-.02,.03,sdf)
sdf<-shiftdf("tall_p3_1",-0.1,-0.05,sdf)
sdf<-shiftdf("tall_p3_2",0.15,-0.05,sdf)
sdf<-shiftdf("cem_1",-0.11,.04,sdf)
sdf<-shiftdf("cem_2",.15,-.01,sdf)
sdf<-shiftdf("jurk_sandar",-.02,-.03,sdf)
df<-data.frame(x=matr[,"PC1",1],y=matr[,"PC3",1],categories=swapFun(categories),Conditions=swapFunB(categories))
postext<-shiftCols(df$x,df$y,categories,sdf)
p<-ggplot(df,aes(x=x,y=y,col=Conditions,label=categories))+ylab("PC3")+xlab("PC1")+geom_point(size=3,shape=1)+theme_bw()+geom_text(x=postext$x,y=postext$y,show_guide=F,size=4)+scale_color_manual(values=c("green","red" , "orange","blue","brown"))


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
sdf<-shiftdf("jurk_sandar_1",.1,0.05,sdf)
sdf<-shiftdf("rpmi_1",.1,0.05,sdf)
sdf<-shiftdf("rpmi_2",.15,.02,sdf)
sdf<-shiftdf("tall_p1",.01,-.03,sdf)
sdf<-shiftdf("tall_p2_1",.15,.02,sdf)
sdf<-shiftdf("tall_p2_2",.2,0,sdf)
sdf<-shiftdf("jurk",-.02,.03,sdf)
sdf<-shiftdf("tall_p3_1",-0.1,-0.05,sdf)
sdf<-shiftdf("tall_p3_2",0.15,-0.05,sdf)
sdf<-shiftdf("cem_1",-0.11,.04,sdf)
sdf<-shiftdf("cem_2",.15,-.01,sdf)
sdf<-shiftdf("jurk_sandar",-.02,-.03,sdf)
n=1
p<-plotPCMat(matr,c("PC1","PC3"),categories,swapFun,swapFunB,n=n,sdf=sdf)
save(p,matr,swapFun,swapFunB,n,categories,sdf,file="~/Dropbox/22x22-pca-1-3-pvalue=0.RData")

plotPCMat<-function(matr,pcs,categories,swapFun,swapFunB,n=1,sdf=NULL,colours=c("green","red" , "orange","blue","brown"),legend=TRUE,text=TRUE){    
    df<-data.frame(x=matr[,pcs[1],n],y=matr[,pcs[2],n],categories=swapFun(categories),Conditions=swapFunB(categories))
    if(is.null(sdf))
       postext<-df
   else
       postext<-shiftCols(df$x,df$y,categories,sdf)
    p<-ggplot(df,aes(x=x,y=y,col=Conditions,label=categories))+ylab(pcs[2])+xlab(pcs[1])+geom_point(size=3,shape=1)+theme_bw()
    if(! is.null(text))
        p<-p+geom_text(x=postext$x,y=postext$y,show_guide=F,size=4)
    if(is.null(legend))
        p<-p+ theme(legend.position="none")
    if(! is.null(colours))
        p<-p+scale_color_manual(values=colours)
    p
}

pass<-function(x) x


sdf<-shiftdf("ecfc.tsa",0,-0.05)
sdf<-shiftdf("cd34",-0.05,-0.05,sdf)
sdf<-shiftdf("cd133",0,-0.05,sdf)
sdf<-shiftdf("cd34_new",-.2,0,sdf)
sdf<-shiftdf("meka",0,-.05,sdf)
sdf<-shiftdf("eryt_f",0,.05,sdf)
sdf<-shiftdf("eryt",0,.05,sdf)
sdf<-shiftdf("eryt_a",-.17,0,sdf)
sdf<-shiftdf("k562_1",-.15,0,sdf)
sdf<-shiftdf("k562_2",-.15,.01,sdf)
sdf<-shiftdf("jurk_sandar_1",.2,0.05,sdf)
sdf<-shiftdf("rpmi_1",0,-0.05,sdf)
sdf<-shiftdf("rpmi_2",0,.05,sdf)
sdf<-shiftdf("tall_p1",.15,0,sdf)
sdf<-shiftdf("tall_p2_1",.17,-.02,sdf)
sdf<-shiftdf("tall_p2_2",.2,0,sdf)
sdf<-shiftdf("jurk",.1,-.02,sdf)
sdf<-shiftdf("tall_p3_1",0.2,0.025,sdf)
sdf<-shiftdf("tall_p3_2",0.2,-0.025,sdf)
sdf<-shiftdf("cem_1",0.15,0,sdf)
sdf<-shiftdf("cem_2",.15,0.02,sdf)
sdf<-shiftdf("jurk_sandar",-.02,-.03,sdf)

n=1

pcsO<-paste("PC",seq(16),sep="")

plots<-lapply(list(c(1,2),c(1,3),c(1,5),c(1,7)),function(i){
    plotPCMat(matr,c(pcsO[i[1]],pcsO[i[2]]),categories,swapFun,swapFunB,n=n,sdf=NULL,text=NULL,legend=NULL)})
multiPlot(plots,2,2)

plots<-lapply(seq(4,12),function(i){
    plotPCMat(matr,c(pcsO[1],pcsO[4]),categories,swapFun,swapFunB,n=i,sdf=NULL,text=NULL,legend=NULL)})


multiPlot<-function(plots,rows,cols){
pushViewport(viewport(layout=grid.layout(rows,cols)))
for(i in seq(length(plots))){
    r<-modulous((i-1),cols)+1
    c<-floor((i-1)/cols)+1
    print(plots[[i]],vp=viewport(layout.pos.row=c,layout.pos.col=r))
    if(i>=rows*cols)
        break
}
}


#save(p,matr,swapFun,swapFunB,n,categories,sdf,file="~/Dropbox/22x22-pca-3-5-pv=0.RData")




plotPCs(prc$eigenVectors,c(1,3),prc$normData,categories)



#categoriesB<-colnames(read.table(mfn4(pvalues[1]),header=T)[4:7])

for(i in c(1,5,10,15,25,30)){
    tags<-read.table(mfn4(pvalues[i]),header=T)
    #tags<-renameColumn(tags,"ecfc.tsa","ecfc-tsa")
    tags<-tags[,4:dim(tags)[2]]    
    tags<-(function(){
        n<-1
        newMatrix<-do.call(cbind,lapply(seq(length(categoriesB)),function(i){ if(categoriesB[i] %in% colnames(tags)) {r<-tags[,n]; n<<-n+1;r} else rep(0,dim(tags)[1])}))
        colnames(newMatrix)<-categoriesB
        newMatrix
    })()
    #tags<-tags[,categoriesB]
    heights<-read.table(mfn4(pvalues[i],TRUE),header=T)[,unlist(sapply(categoriesB,function(x) which(x==categories)))]
    prc<-pca(heights)
    t2<-mergeFun(tags)
    for(pc in c(1,2,3)){
        fname<-paste("~/Dropbox/4x4-contrib-",pc,"-pvalue=",pvalues[i],"-peak.png",sep="")   
        p<-stackedContrib(prc$eigenVectors[,pc],"contrib2",t2,swapFun=swapFunB,colourOveride=swapFunC,sum=FALSE)
        #ggsave(fname,plot=p)
  #      dev.off()
    }
}



for(i in c(3,5,7,9,11,13,15,17,19,21))
    for(j in seq(10))
        ggsave(p[[i]][[j]],filename=paste("~/Dropbox/")

    








