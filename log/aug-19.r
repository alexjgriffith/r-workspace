#source("http://bioconductor.org/biocLite.R")
library(CCCA)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
source(system.file("scripts","eboxFrequency.r",package="CCCA"))
source("~/r-workspace/project/ccca.r")
source("~/r-workspace/project/project.r")
source("~/r-workspace/project/project-variables.r")

peakFiles<-paste0("~/Dropbox/UTX-Alex/Paper/Raw\ Data/peaks/",categories,"~combined_mock_peaks.xls")

env<-getPRC20(2)
contPeaks<-mergeFun(env$over[,4:dim(env$over)[[2]]],swapFunD)
unique<-rowSums(contPeaks)==1


peaksBED<-addNames(lapply(contexts,function(cont){i<-swapFunD(categories)==cont
       makeAFS(peakFiles[i],categories[i],pValue=20)[,c(1,2,3)]}
       ),contexts,list=TRUE)

peaksUNI<-addNames(lapply(contexts,function(cont) env$bed[unique&as.logical(contPeaks[,cont]),]),contexts,list=TRUE)

peaksPCA<-addNames(lapply(contexts,function(cont) env$over[env$reg[,cont],c(1,2,3)]),contexts,list=TRUE)

    
envFromPeaks<-function(peaks){
    contexts<-names(peaks)
    combSum<-function(a,b) c(a,b+a[length(a)])
    reglength<-Reduce(combSum,sapply(peaks, function(x)dim(x)[[1]]))
    BED$reg<-matrix(FALSE,dim(do.call(rbind,peaks))[1],length(peaks))
    k<-1
    for(i in seq(length(peaks))){
        BED$reg[k:reglength[i],i]=TRUE
        k<-reglength[i]+1
    }
    colnames(BED$reg)<-contexts
    BED$bed<-do.call(rbind,peaks)
    BED<-addFasta(BED)
    BED
}


BED<-envFromPeaks(peaksBED)

UNI<-envFromPeaks(peaksUNI)

PCA<-envFromPeaks(peaksPCA)

stem(seq(8,20),sapply(seq(8,20),function(i)length(findMotifDist(PCA$fasta,"GATAA","CANNTG",PCA$reg[,"Erythroid"],i+1,i+2,i+3))))

ops<-list(list("GNN",
               function(fasta){
                   all<-rep(TRUE,length(fasta))
                   locations<-findMotifDist(fasta,"GATAA","CANNTG",all,13,14,15)
          makeLogic(locations,length(fasta)))
},
     list("GG",
          function(fasta) makeLogic(grepMotifs("CAGGTG",fasta),length(fasta))))

addMotifs<-function(env,ops,name="motifs"){
    env[[name]]=
        do.call(cbind,lapply(ops,function(o)
        o[[2]](env$fasta)))
    colnames(env[[name]])<-lapply(ops, function(x) x[[1]])
    env
}

PCA=addMotifs(PCA,ops)
BED=addMotifs(BED,ops)
UNI=addMotifs(UNI,ops)

save(PCA,BED,UNI,file="~/Dropbox/UTX-Alex/Paper/Raw Data/peakLocations.RData")
