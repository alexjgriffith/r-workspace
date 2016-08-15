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



id="name2"
genes20pca2<-geneMatrix(over20,reg20FSel(2),regions,geneList,id)

genes5pca2<-geneMatrix(over5,reg5FSel(2),regions,geneList,id)

getGenes<-function(gn){
    temp<-lapply(colnames(gn),function(x,gn) rownames(gn[gn[,x]==1,]),gn)
    names(temp)<-colnames(gn)
    temp
}

saveGenes<-function(gn,sou,sd,p,name="name",dir="~/Dropbox/UTX-Alex/br-data/genes"){
    filename<-paste0("gene_",name,"_",sou,"_pvalue=",p,"_sd=",sd,".txt")
    filename<-paste0(dir,"/",filename)
    write.table(gn[[sou]],filename,col.names = FALSE,row.names=FALSE,quote = FALSE)    
}

id="name2"
genes20pca2<-geneMatrix(over20,reg20FSel(2),regions,geneList,id)
genes5pca2<-geneMatrix(over5,reg5FSel(2),regions,geneList,id)
saveGenes(getGenes(genes20pca2),"ECFC",2,20)    
saveGenes(getGenes(genes20pca2),"Erythroid",2,20)    
saveGenes(getGenes(genes20pca2),"HSC",2,20)
saveGenes(getGenes(genes20pca2),"Leukemia",2,20)
saveGenes(getGenes(genes5pca2),"ECFC",2,5)    
saveGenes(getGenes(genes5pca2),"Erythroid",2,5)    
saveGenes(getGenes(genes5pca2),"HSC",2,5)
saveGenes(getGenes(genes5pca2),"Leukemia",2,5)


genes20pca2<-geneMatrix(over20,reg20FSel(3),regions,geneList,id)
genes5pca2<-geneMatrix(over5,reg5FSel(3),regions,geneList,id)
saveGenes(getGenes(genes20pca2),"ECFC",3,20)    
saveGenes(getGenes(genes20pca2),"Erythroid",3,20)    
saveGenes(getGenes(genes20pca2),"HSC",3,20)
saveGenes(getGenes(genes20pca2),"Leukemia",3,20)
saveGenes(getGenes(genes5pca2),"ECFC",3,5)    
saveGenes(getGenes(genes5pca2),"Erythroid",3,5)    
saveGenes(getGenes(genes5pca2),"HSC",3,5)
saveGenes(getGenes(genes5pca2),"Leukemia",3,5)

genes20pca2<-geneMatrix(over20,reg20FSel(1),regions,geneList,id)
genes5pca2<-geneMatrix(over5,reg5FSel(1),regions,geneList,id)
saveGenes(getGenes(genes20pca2),"ECFC",1,20)    
saveGenes(getGenes(genes20pca2),"Erythroid",1,20)    
saveGenes(getGenes(genes20pca2),"HSC",1,20)
saveGenes(getGenes(genes20pca2),"Leukemia",1,20)
saveGenes(getGenes(genes5pca2),"ECFC",1,5)    
saveGenes(getGenes(genes5pca2),"Erythroid",1,5)    
saveGenes(getGenes(genes5pca2),"HSC",1,5)
saveGenes(getGenes(genes5pca2),"Leukemia",1,5)


rownames(genes5pca2[genes5pca2[,"HSC"]==1,])



lapply(colnames(genes20pca),writegenes,genes20pca,paste("~/Dropbox/genes20pca-",id,"-",sep=""))
