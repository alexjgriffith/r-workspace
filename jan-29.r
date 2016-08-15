## working to redo read density matirix generation for multiple data sets
## the issue is that the bed file filtering failed
## Step 1 create unified peak list for pvalues from 5-75 at steps of 2.5
## Step 2 create density files for each
## Step 3 plot stability as a function of p-value
## Step 4 look at changes in pca plot

source("~/r-workspace/nov-functions.r")
options(max.print=200)

fullAFS<-wrapAFS(makePeakFiles,categories)
pvalues<-seq(from=5,to=75,by=2.5)
n=length(pvalues)
cs<-makeForkCluster(n,renice=0)
ret<-parallel::parLapply(cs,pvalues,fullAFS)
stopCluster(cs)


genCommand<-function(x,name){
    paste0("grep ",x$chr," ",name,"/",name,"_unique_nodupes.bed | awk '$3>",x$start," && $2<",x$end,"{print $0}' | wc -l\n")
}



m=(20-5)/2.5+1
temp<-ret[[m]][1,1:3]

sapply(genCommand(temp,categories),function(x) system(x))

data.frame(heights=as.numeric(strsplit("74 1 1 0 0 2 11 5 6 1 0 0 20 15 7 2 13 6 3 1 67 64"," ")[[1]]),categories=categories)



rawDataFiles<-function(categories)paste("/mnt/brand01-00/mbrand_analysis/data_sets/",categories,"/",categories,"_sorted.bed",sep="")

#score1<-pileUp(ret[[1]],rawDataFiles(categories),20,TRUE)
score7<-pileUp(ret[[7]],rawDataFiles(categories),20,TRUE)

ret[[1]][1:3]

order<-unique(as.character($chr))

#orderBed(ret[[7]])[,1:3]

#score7<-pileUp(addColnames(temp,c("chro","start","end")),rawDataFiles(categories),20,TRUE)



library(CCCA)
source("~/r-workspace/nov-functions.r")

orderBed<-function(ret)
    ret[order(as.character(ret[,1]),ret[,3]),]

fullAFS<-wrapAFS(makePeakFiles,categories)


ret<-fullAFS(20);
bed=orderBed(ret)[,1:3]

file="/mnt/brand01-00/mbrand_analysis/data_sets/cd34/cd34_sorted.bed"



score=as.integer(rep(0,length(bed$chr)))

results<-.C("pileup",file,chrom=as.character(bed$chr)
           ,start=as.integer(bed$start),end=as.integer(bed$end),peaknum=as.integer(length(bed$chr)-1),score=score)


rawDataFiles<-function(categories)paste("/mnt/brand01-00/mbrand_analysis/data_sets/",categories,"/",categories,"_sorted.bed",sep="")


score7<-pileUp(bed,rawDataFiles(categories),20,TRUE)

getPileUp(file,bed,as.character(bed$chr),length(bed$chr))

source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/jan-functions.r")
source("~/r-workspace/jan-variables.r")


bed<-orderBed(bed20)

heights<-addColnames(read.table("~/r-workspace/test.txt"),categories)
prc=pca(heights)

eryt=bed[normalize(prc$eigenVectors[,1])>2,]
jurk=bed[normalize(prc$eigenVectors[,1])<(-2),]
ecfc=bed[normalize(prc$eigenVectors[,4])<(-2),]
stem=bed[normalize(prc$eigenVectors[,4])>2,]

#plotPCMat2D(pca2Matr(prc),c(1,4),categories,swapFun,swapFunB,swapFunC)

printBB(eryt,"~/Dropbox/eryt_sd_2.bb")
printBB(jurk,"~/Dropbox/jurk_sd_2.bb")
printBB(ecfc,"~/Dropbox/ecfc_sd_2.bb")
printBB(stem,"~/Dropbox/jurk_sd_2.bb")

tts<-function(x){
    paste0(x[,1],":",x[,2]-2000,"-",x[,3]+2000)
}

write.table(tts(eryt),"~/Dropbox/UTX-Alex/br-data/eryt_sd_2.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(tts(jurk),"~/Dropbox/UTX-Alex/br-data/jurk_sd_2.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(tts(ecfc),"~/Dropbox/UTX-Alex/br-data/ecfc_sd_2.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(tts(stem),"~/Dropbox/UTX-Alex/br-data/stem_sd_2.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
