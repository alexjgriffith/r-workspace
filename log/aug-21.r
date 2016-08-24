library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(CCCA)

geneFile<-system.file("exdata","hg19.RefSeqGenes.csv",package="CCCA")

geneList<-read.table(geneFile,comment="",header=T)




extN<-lapply(seq(1000), function(i) {
    ext<-geneList[i,c("exonStarts","exonEnds")]
    cbind(as.numeric(strsplit(as.vector(ext[[1]]) ,",")[[1]]),
      as.numeric(strsplit(as.vector(ext[[2]]) ,",")[[1]] ))})

lapply(extN, function(x) length(x))

exonGRanges<-getSeq(BSgenome.Hsapiens.UCSC.hg19,"chr1",extN[,1],extN[,2])

exon<-as.vector(exonGRanges)

options<-as.factor(apply(cbind(unlist(lapply(c("A","C","G","T"),rep,16))
,rep(sapply(c("A","C","G","T"),rep,4),4)
,rep(c("A","C","G","T"),16)),1,paste0,collapse=""))

factor(sapply(seq(floor(length(exon)/3)),
       function(x)
           paste0(as.character(exon[c((x-1)*3+1,(x-1)*3+2,(x-1)*3+3)]),collapse="")),levels=levels(options))

list(list(options[1],options[2],2,1),
     list(options[3],options[4],3,),
     list(options[5],options[6],1))


opts<-list("GGGNNNAGT",
           "TATCCT",
           "TGCTCT")

sapply(opts,function(x)do.call(cbind,gregexpr(IUPACtoBase(x),exonGRanges))!=-1
       )
