library(CCCA)
library(parallel)
source("~/r-workspace/utils.r")

source("~/r-workspace/project.r")

categories<-readCategories("~/Dropbox/Data/categories/22-categories.txt")

conts<-readCategories("~/Dropbox/UTX-Alex/jan/contcats")

eqiv<-gsub("_mock","",conts)

catfiles<-paste0("~/Dropbox/Data/august_peaks/",categories,"~combined_mock_peaks.xls")

cfiles<-paste0("/mnt/brand01-00/mbrand_analysis/data_sets/", conts,"/",conts,"_unique_nodupes.bed")

makeEnv<-function(categories,eqiv,catfiles,contfiles){
    afs<-makeAFS(catfiles,categories,pValue=20)
    n<-min(c(20,length(cfiles)))
    cl<-makeForkCluster(n)
    score<-addColnames(CCCA::pileUp(afs,contfiles,n=n,clust=cl),eqiv)
    stopCluster(cl)    
    prc<-pca(score)
    list(categories=categories,control=eqiv,over=afs,bed=afs[,1:3],heights=score,prc=prc)
}

cont20<-getPRC20F(makeEnv(categories,eqiv,catfiles,contfiles),2)

save(cont20,file="~/Dropbox/Data/cont20.RData")

load("~/Dropbox/Data/cont20.RData")


source("~/r-workspace/pca.r")

plotPCMat2D(pca2Matr(cont20$prc),c("PC1","PC4"),cont20$control,swapFun,swapFunB,swapFunC)


sumHeights<-function(heights,categories,cats=NULL,swapFun=pass){
    if(is.null(cats))
        cats<-unique(swapFun(categories))
    selSum<-function(x)
        apply(cbind(heights[,x== swapFun(categories)],0),1,sum)
    addColnames(do.call(cbind,lapply(cats,selSum)),cats)
}



## histograms
png("~/Dropbox/UTX-Alex/bio-results-images/bio-results-1-cont.png")
width=0.5
xlim=c(-2,8)
par(mfrow=c(2,2),mar=c(2,2,3,2),oma=c(0, 0, 2, 0))
cats<-c("ECFC","Erythroid","HSC","Leukemia")
sh<-sumHeights(cont20$prc$normData,cont20$control,swapFun=swapFunD)
for( v in lzip(cats,c("green","red","orange","blue"))){
    x<-log(sh[,v[[1]]])
    breaks<-diff(range(Filter(function(x)x!=-Inf,x)))/width
    hist(x,breaks,xlim=xlim,main=v[[1]],ylab="",xlab="",col=v[[2]])
}
dev.off()


env<-getPRC20(2)


#png("~/Dropbox/UTX-Alex/bio-results-images/bio-results-1-cont-reg.png")
width=0.5
xlim=c(-2,8)
par(mfrow=c(2,2),mar=c(2,2,3,2),oma=c(0, 0, 2, 0))
cats<-c("ECFC","Erythroid","HSC","Leukemia")
for( v in lzip(cats,c("green","red","orange","blue"))){
    sh<-sumHeights(cont20$prc$normData[env$reg[,v[[1]]],]
                  ,cont20$control,swapFun=swapFunD)
    x<-log(sh[,v[[1]]])
    breaks<-diff(range(Filter(function(x)x!=-Inf,x)))/width
    hist(x,breaks,xlim=xlim,main=v[[1]],ylab="",xlab="",col=v[[2]])
}
#dev.off()

png("~/Dropbox/UTX-Alex/bio-results-images/bio-results-1-e.png")
width=0.2
xlim=c(3,8)
par(mfrow=c(2,2),mar=c(2,2,3,2),oma=c(0, 0, 2, 0))
cats<-c("ECFC","Erythroid","HSC","Leukemia")
for( v in lzip(cats,c("green","red","orange","blue"))){
    sh<-sumHeights(env$prc$normData[env$reg[,v[[1]]],]
                  ,env$categories,swapFun=swapFunD)
    x<-log(sh[,v[[1]]])
    breaks<-diff(range(Filter(function(x)x!=-Inf,x)))/width
    b<-breaks
    hist(x,b,xlim=xlim,main=v[[1]],ylab="",xlab="",col=v[[2]])
}
dev.off()


png("~/Dropbox/UTX-Alex/bio-results-images/bio-results-1-d.png")
width=0.5
xlim=c(-2,8)
par(mfrow=c(2,2),mar=c(2,2,3,2),oma=c(0, 0, 2, 0))
    x<-log(sumHeights(env$prc$normData
                  ,env$categories,swapFun=swapFunD)[,"ECFC"])
    breaks<-diff(range(Filter(function(x)x!=-Inf,x)))/width
    b<-breaks
    hist(x,b,xlim=xlim,main=v[[1]],ylab="",xlab="",col="green")
    x<-log(sumHeights(env$prc$normData[normalize(env$prc$eigenVectors[,2])<(-2),]
                  ,env$categories,swapFun=swapFunD)[,"ECFC"])
    breaks<-diff(range(Filter(function(x)x!=-Inf,x)))/width
    b<-breaks
    hist(x,b,xlim=xlim,main=v[[1]],ylab="",xlab="",col="green")
    x<-log(sumHeights(env$prc$normData[normalize(env$prc$eigenVectors[,4])<(-2),]
                  ,env$categories,swapFun=swapFunD)[,"ECFC"])
    breaks<-diff(range(Filter(function(x)x!=-Inf,x)))/width
    b<-breaks
    hist(x,b,xlim=xlim,main=v[[1]],ylab="",xlab="",col="green")
    x<-log(sumHeights(env$prc$normData[normalize(env$prc$eigenVectors[,4])<(-2)&normalize(env$prc$eigenVectors[,2])<(-2),]
                  ,env$categories,swapFun=swapFunD)[,"ECFC"])
    breaks<-diff(range(Filter(function(x)x!=-Inf,x)))/width
    b<-breaks
    hist(x,b,xlim=xlim,main=v[[1]],ylab="",xlab="",col="green")
dev.off()


# define env
# env (categories, over, bed, heights,prc)

# do AFS using the treatment peaks and the background data sets

n<-4
hsc<-apply(mergeFun(env$over[normalize(env$prc$eigenVectors[,4])>(n),4:dim(env$over)[2]],swapFunD),2,sum)
nhsc<-apply(mergeFun(env$over[normalize(env$prc$eigenVectors[,4])<(-n),4:dim(env$over)[2]],swapFunD),2,sum)
rbind(nhsc/sum(nhsc),hsc/sum(hsc))

unlist(addNames(lapply(cats,function(x) sum(hsc[x==swapFunB(categories)])/sum(hsc)),cats,list=TRUE))

unlist(addNames(lapply(cats,function(x) sum(nhsc[x==swapFunB(categories)])/sum(nhsc)),cats,list=TRUE))

## count how many specific contexts are in each set

source("~/Masters/CCCA/inst/scipts/eboxFrequency.r")


library(BSgenome.Hsapiens.UCSC.hg19)

env$fasta<-getSeq(BSgenome.Hsapiens.UCSC.hg19,env$bed$chr,start=env$bed$start+150,width=300)

motifs<-sapply(c(genEboxCombs(),"CANNTG","GATAA","ACCACA","AGGCGG","SCACTG","TTGTTM","NNNN"),IUPACtoBase)

cList<-sapply(motifs,compliment)
lM2<-lapply(motifs,grep,env$fasta)
lC2<-lapply(cList,grep,env$fasta)

hE<-motifHist(env$fasta,motifs,cList,lM2,lC2,12,11,env$reg[,"Erythroid"])

# make table for presence in the 6 contexts

occur<-mapply(function(y,yc) apply(env$reg,2,function(x) sum(x & (makeLogic(y,length(x))|makeLogic(yc,length(x))))), lM2,lC2)

write.table(occur,"~/Dropbox/UTX-Alex/br-data/motif-occurance.txt")

do.call(rbind,lapply(seq(4),function(x) apply(getPRC20F(list(prc=env$prc),x)$reg,2,sum)))

o<-addRownames(apply(env$reg,2 ,function(r) do.call(cbind,lapply(1:11, function(x){
    #h<-motifHist(env$fasta,motifs,cList,lM2,lC2,x,12,r)
    #sum(h>-17 & h<(-14))
    sum(makeLogic(union(intersect(lM2[[x]],lM2[[12]]),intersect(lC2[[x]],lC2[[12]])),length(r)) & r)    
    #length(h)
}))),motifs[1:11])


or<-addRownames(apply(env$reg,2 ,function(r) do.call(cbind,lapply(1:11, function(x){
    #h<-motifHist(env$fasta,motifs,cList,lM2,lC2,x,12,r)
    #sum(h>-17 & h<(-14))
    sum(makeLogic(union(intersect(lM2[[x]],lM2[[13]]),intersect(lC2[[x]],lC2[[13]])),length(r)) & r)    
    #length(h)
}))),motifs[1:11])


aw<-addRownames(apply(env$reg,2 ,function(r) do.call(cbind,lapply(1:17, function(x){
    #h<-motifHist(env$fasta,motifs,cList,lM2,lC2,x,12,r)
    #sum(h>-17 & h<(-14))
    sum(makeLogic(union(lM2[[x]],lC2[[x]]),length(r)) & r)    
    #length(h)
}))),sapply(motifs[1:17],consenusIUPAC))



       
