library(CCCA)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
#library(org.Hs.eg.db)
#library(RDAVIDWebService)
#library(xlsx)
#library(org.Hs.eg.db)
source(system.file("scripts","eboxFrequency.r",package="CCCA"))
source("~/r-workspace/project/ccca.r")
source("~/r-workspace/project/project.r")
source("~/r-workspace/project/project-variables.r")

options(java.parameters = "-Xmx1024m")

library(xlsx)


load("~/Dropbox/UTX-Alex/Paper/Raw Data/peakLocations.RData")

eboxPref<-function(env){
    motifs<-c(genEboxCombs(),"CANNTG","NNNNNN")
    eboxs<-do.call(rbind,lapply(contexts,function(reg,env) sapply(motifs,function(x)length(grepMotifs(x,env$fasta[env$reg[,reg],]))),env))
    eboxs<-cbind(eboxs,noEbox=eboxs[,"NNNNNN"]-eboxs[,"CANNTG"])

    eboxsRank<-rbind(apply(eboxs[,1:10],1,function(x) abs(rank(x)-11)),0,0,0)

    eboxsRatio<-apply(eboxs,1,function(x) x/x[length(x)-1])
    eboxsF<-t(eboxs)
    colnames(eboxsF)<-contexts
    colnames(eboxsRatio)<-contexts
    colnames(eboxsRank)<-contexts
    eboxTable<-do.call(cbind,lapply(contexts,function(cont)cbind(eboxsF[,cont],eboxsRatio[,cont],eboxsRank[,cont])))
    colnames(eboxTable)<-paste(sapply(contexts,rep,3),c("No","Ratio","Order") ,sep=".")
    eboxTable
}


addEboxSheet<-function(book,cellName,env,reg,motifs){
    whichBed<-function(motif,env,reg)
        env$bed[grepMotifs(motif,env$fasta[env$reg[,reg],]),]
    s<-createSheet(book,cellName)

    for(i in seq(length(motifs))){
        da<-whichBed(motifs[i],env,reg)
        column=(i-1)*3+1
        addDataFrame(data.frame(da),row.names=FALSE,sheet=s,startColumn = column ,startRow = 2)
        addDataFrame(motifs[i],col.names=FALSE,row.names=FALSE,sheet=s,startColumn = column ,startRow = 1)
    }
    i=length(motifs)+1
    no<-! makeLogic(unique(unlist(sapply(motifs,grepMotifs,env$fasta[env$reg[,reg],]))),length(env$fasta[env$reg[,reg],]))
    column=(i-1)*3+1
    addDataFrame(env$bed[no,],row.names=FALSE,sheet=s,startColumn = (i-1)*3+1,startRow = 2)
    addDataFrame(data.frame("None"),col.names=FALSE,row.names=FALSE,sheet=s,startColumn = column ,startRow = 1)
}

addPage<-function(book,name,df){
    s<-createSheet(book,name)
    addDataFrame(data.frame(df),sheet=s)
}

motifs<-c(genEboxCombs(),"CANNTG","NNNNNN")

eboxBED<-eboxPref(BED)
eboxUNI<-eboxPref(UNI)
eboxPCA<-eboxPref(PCA)

ebook<-xlsx::createWorkbook()
addPage(ebook,"PCA",eboxPCA)
addPage(ebook,"ALL",eboxBED)
addPage(ebook,"UNI",eboxUNI)

gc()
lapply(contexts,function(cont)
    addEboxSheet(ebook,paste0("PCA-",cont),PCA,cont,motifs[1:10]))

gc()
lapply(contexts,function(cont)
    addEboxSheet(ebook,paste0("ALL-",cont),BED,cont,motifs[1:10]))

gc()
lapply(contexts,function(cont)
    addEboxSheet(ebook,paste0("UNI-",cont),UNI,cont,motifs[1:10]))

gc()
saveWorkbook(ebook,path.expand("~/Desktop/test.xlsx"))

    
