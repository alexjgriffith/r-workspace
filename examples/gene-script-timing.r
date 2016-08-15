heightFile<-"peaks/jan/combined_heights.bed"
fileLocation<-"~/"
geneList<-read.delim(paste(fileLocation,"hg19.RefSeqGenes.csv",sep=""))
data<-readSplitTable(heightFile)$data

chrom<-as.character(geneList$chrom)
tss<-as.numeric(geneList$txStart)

strand<-geneList$strand
levels(strand)<-c(-1,1)
strand<-as.numeric(strand)


genomicRegions<-function(chrom,tss,strand,proxUp,proxDown,distal)
    geneRangeToRegion(geneRanges(chrom,tss,strand,5000,1000,50000),chrom)



regions<-genomicRegions(chrom,tss,strand,5000,1000,100000)


reg<-mapply(function(pc,loc)buildRegions(data,pc,loc)[,1],
            list(1,1,3,c(3,5),c(3,7),c(3,7),7),
            list("top","bottom","top",c("top","top"),c("top","bottom"),c("top","bottom"),"top"))


Sys.time()-start
#
jurkGenes<-peakGeneRegions(bedData[reg[,2],],a,geneList)
filenames<-lapply(cbind("erythroid","t-all","ecfc","other","hspc","meka","diff"),paste, "-genes.txt",sep="")
for (i in seq(7)){
   genes<-peakGeneRegions(bedData[reg[,i],],a,geneList)
   write.table(genes,filenames[[i]],quote=FALSE,col.names=FALSE,row.names=FALSE)}
