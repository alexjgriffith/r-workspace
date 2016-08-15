getwd()
setwd()

library("Biostrings")
library("stringr")
library("motifRG")


member<-function(posibList,mem,test="default"){
    fun<-switch(test,default=function(x,y){all(y %in% x)})
    sapply(posibList,fun,mem)}

memberWrapper<-function(bedData,lis ,n=4,sep="-"){
    member(strsplit(as.character(bedData[,n]),sep),lis)}

geneAssoc<-function(point,geneList,bounds){
    # Currently Returns the enhancer locations
    peak<-(as.numeric(point[[2]])+as.numeric(point[[3]]))/2
    tss<-geneList$txStart
    ess<-geneList$txEnd
    chrom<-geneList$chrom
    a<-which(chrom==as.character(point[[1]]))
    b<-a[which(tss[a]-bounds[1]<peak)]
    r<-b[which(ess[b]+bounds[3]>peak)]
    e<-r[c(which(tss[r]-bounds[2]>peak),which(ess[r]+bounds[4]<peak))]
    e}
       
minGene<-function(x,y,start,miz=1){
    x[which(order(abs(start[x]-y))<=miz)]}

geneAssociation<-function(bedData,geneList,bounds,miz=1){
    locations<-apply(bedData,1,function(x){geneAssoc(x,geneList,bounds)})
    peak<-(as.numeric(bedData[,2])+as.numeric(bedData[,3]))/2
    tssGenes<-mapply(minGene,locations,peak,MoreArgs=list(start=geneList$txStart,miz=miz))
    lapply(tssGenes, function(x){as.character(geneList$name2[x])})}



fileLocation="/home/agriffith/Dropbox/UTX-Alex/jan/"

bedData<-read.delim(paste(fileLocation,"combined_sorted.bed",sep=""),header=0)
geneList<-read.delim("hg19.RefSeqGenes.csv")
#allFasta<-readDNAStringSet(paste(fileLocation,"combined.fasta",sep=""),use.names=TRUE)
#foreg<-allFasta[member(strsplit(as.character(types[,4]),"-"),c("k562","eryt"))]
#backg<-allFasta[member(strsplit(as.character(types[,4]),"-"),c("jurk","cem","rpmi"))]

geneAssociation(bedData[1:100,],geneList,  c(50000,2000,2000,50000))

strList<-"(intersection jurk cem rpmi (not (union k562 eryt)))"

subList<-list("union","jurk","cem","rpmi",list("not",list("union","k562","eryt")))

for(i in subList[2:length(subList)]){if(is.list(i)){print i}}
