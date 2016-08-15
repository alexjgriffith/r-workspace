
library(CCCA)

data<-read.table("~/Dropbox/Data/thesis-november/4x4-0-overlap-matrix.data",header=TRUE)

bedtag<-read.table("~/Dropbox/Data/thesis-november/4x4-0-score-matrix.data",header=TRUE)


tags<-bedtag[,4:dim(bedtag)[2]]

pcs<-pca(data)
plotPCs(pcs$eigenVectors,c(1,2),pcs$normData,names(data))


library(ggplot2)
library(parallel)

splitZip<-function(inList)
    list(sapply(inList,"[",1),sapply(inList,"[",2))

createZip<-function(x)
    Map(function(i,j)cbind(x[i],x[j]),seq(1,length(x),2),seq(2,length(x),2))


getSwapCats<-function(cats,names){    
   swCat<-cats
   names(swCat)<-names
   function(x){as.character(unlist(lapply(x,function(x) swCat[x])))}
}

stringToSwap<-function(x)do.call(getSwapCats,splitZip(createZip(strsplit(x," ")[[1]])))
prepApp<-function(prefix,suffix,...){function(x) paste(prefix,x,suffix,...)}

swapFun<-stringToSwap(paste (rev(strsplit("rpmi Leukemia tall_p3 Leukemia tall_p2 Leukemia tall_p1 Leukemia tall_p2_1 Leukemia tall_p2_2 Leukemia tall_p3_1 Leukemia tall_p3_2 Leukemia jurk_sandar_1 Leukemia jurk_sandar Leukemia jurk Leukemia rpmi_1 Leukemia rpmi_2 Leukemia cem_1 Leukemia cem_2 Leukemia cem Leukemia ecfc-tsa ECFC ecfc ECFC meka MEKA cd133 HSC cd34 HSC cd34_new HSC eryt Erythroid eryt_f Erythroid eryt_a Erythroid k562 Erythroid k562_1 Erythroid k562_2 Erythroid", " ")[[1]]),collapse=" "))


stackedContrib(pcs$eigenVectors[,2],"contrib2",n=3,tags=tags,swapFun=swapFun)



afs<-makeAFS(unlist(lapply(categories,function(x) paste("~/.peaktemp/",x,"~combined_mock_peaks.xls",sep=""))),categories)

rawDataFile<-sapply(categories, function(x) paste("data_sets/",x,"/",x,"_unique_nodupes.bed",sep=""))

changecolname<-function(frame,labels){colnames(frame)<-labels ;frame}

score<-pileUp(hg19Sort(changecolname(afs[,1:3],c("chro","start","end"))),rawDataFile  ,n=20)
