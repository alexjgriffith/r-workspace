library(CCCA)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

jpwmsp<-loadPWM("~/masters/r-workspace/CCCA/inst/exdata/jaspar.pwm","jaspar")[,3]

jpwms<-lapply(jpwmsp,function(x)matrix(c(x),nrow=4,byrow=TRUE))
jnames<-sapply(jpwms,CCCA::PWMtoCons)
jnames<-gsub("^N*||N*$","",jnames)


jl<-sapply(jnames,function(x)length(grepMotifs(x,PCA$fasta)))

jsublM<-jnames[order(jl,decreasing=TRUE)]

jid<-loadPWM("~/masters/r-workspace/CCCA/inst/exdata/jaspar.pwm","jaspar")[,2]
jfid<-unlist(lapply(strsplit(unlist(jid),"\t"),function(x) x[2]))


jaspar<-list(jnames=jnames,jfid=jfid,jl=jl,jsublM=jsublM)

save(jaspar,file="~/r-workspace/data/jaspar.RData")
