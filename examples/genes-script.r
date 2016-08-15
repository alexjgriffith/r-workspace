source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/jan-functions.r")
###source("~/Masters/CCCA/inst/scipts/eboxFrequency.r")
source("~/Masters/mulcal/newR/rGREAT.r")
#source("~/Dropbox/R/makeLatexTable.R")


# generate variables
source("~/r-workspace/jan-variables.r")


geneFile<-"/home/agriffith/Dropbox/UTX-Alex/jan/hg19.RefSeqGenes.csv"
geneList<-read.delim(geneFile)
chrom<-as.character(geneList$chrom)
tss<-as.numeric(geneList$txStart)
strand<-geneList$strand
# double check to make sure that levels(strand) > c("-","+")
levels(strand)<-c(-1,1)
strand<-as.numeric(strand)
#statsAndData<-readSplitTable(heightFile)

regions<-genomicRegions(chrom,tss,strand,1000,5000,1000000)

#env<-getPRC5(2)
env20<-getPRC20(2)
env<-env20;


genes<-geneMatrix(env$over,env$reg,regions,geneList)

lapply(getGenes(genes),length)

lapply(names(getGenes(genes)),function(name){
    saveGenes(getGenes(genes),name,2,5)
})

getGenes<-function(gn){
    temp<-lapply(colnames(gn),function(x,gn) rownames(gn[gn[,x]==1,]),gn)
    names(temp)<-colnames(gn)
    temp
}

saveGenes<-function(gn,sou,sd,p,name="name",dir="~/Dropbox/UTX-Alex/br-data/genes"){
    filename<-paste0("gene_feb8_",name,"_",sou,"_pvalue=",p,"_sd=",sd,".txt")
    filename<-paste0(dir,"/",filename)
    write.table(gn[[sou]],filename,col.names = FALSE,row.names=FALSE,quote = FALSE)    
}


## feb 23
