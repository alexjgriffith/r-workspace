source("~/r-workspace/CCCA.r")
source("~/r-workspace/project.r")

source("~/Masters/CCCA/inst/scripts/eboxFrequency.r")

#source("~/masters/CCCA/inst/scripts/eboxFrequency.r")

# cluster only
# source("~/R/x86_64-unknown-linux-gnu-library/3.2/CCCA/scripts/eboxFrequency.r")

library(BSgenome.Hsapiens.UCSC.hg19)

## Build new all set of peaks
categories<-as.character(unlist(read.table("~/Dropbox/Data/categories/22-categories.txt")))

contexts<-unique(swapFunD(categories))

group<-function(a,b){
    print("reducing")
    print(dim(b))
    rbind(a,b)
}

shiftIfZero<-function(summit,size=300){
    sapply(summit,function(x) {if (x-size<0)
                                   return(size)
                               else
                                   return (x)
                           })
}

makePeaks<-function(peaks){
    print("making peaks")
    print(dim(peaks))
    summit<-shiftIfZero(peaks[,2])
    df<-data.frame(chr=peaks[,1],summit=summit,name=peaks[,3])
    df[order(df$chr,df$summit),]
}

contexts<-unique(swapFunD(categories))

peaks<-lapply(contexts,function(x) makePeaks(do.call(rbind,lapply(categories[x==swapFunD(categories)],function(cat) readPeaksXLS(makePeakFiles(cat),cat)))))



unifiedPeaks<-lapply(peaks,function(sa) unifyBedFile(sa, 700))

mapply(function(peaks,contexts){
    write.table(data.frame(chr=peaks$chr,start=peaks$summit-350,end=peaks$summit-350),paste0("~/Dropbox/UTX-Alex/Paper/all_peaks_",contexts,".bed"),quote=FALSE,col.names=FALSE)
},
    unifiedPeaks,contexts)


unifiedPeaks<-lapply(peaks,function(sa) unifyBedFile(sa, 700))

mapply(function(peaks,contexts){
    write.table(peaks,paste0("~/Dropbox/Data/all_over/all_peaks_",contexts,".bed"),quote=FALSE,col.names=FALSE)
},
       unifiedPeaks,contexts)



### make ebox comb table

normalize<-CCCA:::normalize

env<-getPRC20(2);

env<-addFasta(env)


combs<-c(genEboxCombs(),"CANNTG")

varLocs<-addNames(lapply(contexts,function(cont) addNames(lapply(combs,function(comb)CCCA::grepMotifs(comb,env$fasta[env$reg[,cont],])),combs,list=TRUE)),contexts,list=TRUE)
total<-sapply(varLocs,function(loc) sapply(loc,length))
total<-rbind(total, noEbox=mapply(function(locs,name){ sum(env$reg[,name])-length(locs$CANNTG)}, varLocs,contexts))
ratio<-sweep(total,2,sapply(contexts,function(x) sum(env$reg[,x])),function(a,b) a/b)
rank<-rbind(apply(ratio[1:10,],2,order),0,0)



PCAout<-do.call(cbind,addNames(lapply(contexts,function(conts) data.frame(rank=rank[,conts], no=total[,conts],ratio=ratio[,conts])),contexts,list=TRUE))


write.table(PCAout,"~/Dropbox/UTX-Alex/Paper/Analysis/PCA_ebox_var_table.tab",quote=FALSE)


summitToSE<-function(peaks){
    data.frame(chr=peaks$chr,start=peaks$summit-300,end=peaks$summit+300)
}

env<-list()
env$bed<-summitToSE(do.call(rbind,lapply(addNames(unifiedPeaks,contexts,list=TRUE),function(x) x[,c(1,2)])))


peakLength<-lapply(unifiedPeaks,function(x) dim(x)[1])
peakLengthSummed<-c(0,Reduce(function(a,b) c(a,a[length(a)]+b),peakLength))

env$reg<-addColnames(apply(cbind(peakLengthSummed[1:length(peakLengthSummed)-1],peakLengthSummed[2:length(peakLengthSummed)]),1,function(x) makeLogic(seq(x[1]+1,x[2]),peakLengthSummed[length(peakLengthSummed)])),contexts)

which(apply(env$reg,1,sum)>1)

env<-addFasta(env)

# repeat above process for env
BEDout<-do.call(cbind,addNames(lapply(contexts,function(conts) data.frame(rank=rank[,conts], no=total[,conts],ratio=ratio[,conts])),contexts,list=TRUE))

write.table(BEDout,"~/Dropbox/UTX-Alex/Paper/Analysis/BED_ebox_var_table.tab",quote=FALSE)


## find ebox locations and put them into excell file
library(xlsx)

ework<-createWorkbook()

addEboxSheet<-function(cellName,context){
s<-createSheet(ework,cellName)

for(i in seq(10)){
    da<-env$bed[env$reg[,context],][varLocs[[context]][[i]],]
    column=(i-1)*3+1
    addDataFrame(da,row.names=FALSE,sheet=s,startColumn = column ,startRow = 2)
    addDataFrame(data.frame(combs[i]),col.names=FALSE,row.names=FALSE,sheet=s,startColumn = column ,startRow = 1)
}

i=12
da<-env$bed[env$reg[,context],]
nn<-varLocs[[context]]$CANNTG
no<-setdiff(seq(dim(da)[1]),nn)
column=(i-1)*3+1
addDataFrame(env$bed[no,],row.names=FALSE,sheet=s,startColumn = (i-1)*3+1,startRow = 2)
    addDataFrame(data.frame("noEbox"),col.names=FALSE,row.names=FALSE,sheet=s,startColumn = column ,startRow = 1)
}

addEboxSheet("ENDO-all","ECFC")
addEboxSheet("ERY-all","Erythroid")
addEboxSheet("TALL-all","Leukemia")
addEboxSheet("HSPC-all","HSC")


addEboxSheet("ENDO-pca","ECFC")
addEboxSheet("ERY-pca","Erythroid")
addEboxSheet("TALL-pca","Leukemia")
addEboxSheet("HSPC-pca","HSC")

s<-createSheet(ework,"ebox-pca")
addDataFrame(PCAout,sheet=s,startColumn = 1 ,startRow = 1)

s<-createSheet(ework,"ebox-bed")
addDataFrame(BEDout,sheet=s,startColumn = 1 ,startRow = 1)

saveWorkbook(ework,path.expand("~/Dropbox/UTX-Alex/Paper/TabS4.xlsx"))


## make list of genes
normalize<-CCCA:::normalize

env<-getPRC20(2);

#env<-addFasta(env)

#combs<-c(genEboxCombs(),"CANNTG")

geneFile<-"/home/agriffith/Dropbox/UTX-Alex/jan/hg19.RefSeqGenes.csv"
geneList<-read.delim(geneFile)


chrom<-as.character(geneList$chrom)
tss<-as.numeric(geneList$txStart)
strand<-geneList$strand
# double check to make sure that levels(strand) > c("+","-")
levels(strand)<-c(-1,1)
strand<-as.numeric(strand)
#statsAndData<-readSplitTable(heightFile)

regions<-genomicRegions(chrom,tss,strand,1000,5000,1000000)

PCAgenes<-geneMatrix(env$over,env$reg,regions,geneList,id="name2")

BEDgenes<-geneMatrix(env$bed,env$reg,regions,geneList,id="name2")

## find gg ebox locations for genes on cluster

GGlocs<-grepMotifs("CAGGTG",env$fasta)

NGGlocs<-setdiff(seq(1,length(env$fasta)),GGlocs)


BEDgenesGG<-geneMatrix(env$bed[GGlocs,],env$reg[GGlocs,],regions,geneList,id="name2")

BEDgenesNGG<-geneMatrix(env$bed[NGGlocs,],env$reg[NGGlocs,],regions,geneList,id="name2")

PCAgenesGG<-geneMatrix(env$bed[GGlocs,],env$reg[GGlocs,],regions,geneList,id="name2")

PCAgenesNGG<-geneMatrix(env$bed[NGGlocs,],env$reg[NGGlocs,],regions,geneList,id="name2")

save(PCAgenes,BEDgenes,PCAgenesGG,BEDgenesGG,PCAgenesNGG,BEDgenesNGG,regions,geneList,geneFile,file="~/Dropbox/UTX-Alex/Paper/genes.RData")



selectGenes<-function(matrix,context){
    df<-data.frame(rownames(matrix)[matrix[,context]==1])
    colnames(df)<-context
    df
}


excellSheet<-function(matrix){
    ework<-createWorkbook()
    s<-createSheet(ework,"Genes")
    for(i in seq(dim(matrix)[2])){
        addDataFrame(selectGenes(matrix,colnames(matrix)[i]),col.names=TRUE,row.names=FALSE,sheet=s,startColumn = i)
    }
    for(i in colnames(matrix)){
        createSheet(ework,paste0("BP-",i))
        createSheet(ework,paste0("MF-",i))
    }
    ework    
}

tableGene<-function(matrix,file){
for(i in seq(dim(matrix)[2])){
    write.table(selectGenes(matrix,colnames(matrix)[i]),file,col.names=TRUE,row.names=FALSE,quote=FALSE)
}
}


ework<-excellSheet(PCAgenes)
saveWorkbook(ework,path.expand("~/Dropbox/UTX-Alex/Paper/TabS8(PCA).xlsx"))


tableGene(PCAgenes,path.expand("~/Dropbox/UTX-Alex/Paper/genes_PCA.txt"))
tableGene(PCAgenesGG,path.expand("~/Dropbox/UTX-Alex/Paper/genes_PCAGG.txt"))
tableGene(PCAgenesNGG,path.expand("~/Dropbox/UTX-Alex/Paper/genes_PCANGG.txt"))
tableGene(BEDgenes,path.expand("~/Dropbox/UTX-Alex/Paper/genes_BED.txt"))
tableGene(BEDgenesGG,path.expand("~/Dropbox/UTX-Alex/Paper/genes_BEDGG.txt"))
tableGene(BEDgenesNGG,path.expand("~/Dropbox/UTX-Alex/Paper/genes_BEDNGG.txt"))



ework<-excellSheet(BEDgenes)
saveWorkbook(ework,path.expand("~/Dropbox/UTX-Alex/Paper/TabS7(BED).xlsx"))

ework<-excellSheet(BEDgenesGG)
saveWorkbook(ework,path.expand("~/Dropbox/UTX-Alex/Paper/TabS9(BEDGG).xlsx"))

ework<-excellSheet(BEDgenesNGG)
saveWorkbook(ework,path.expand("~/Dropbox/UTX-Alex/Paper/TabS11(BEDNGG).xlsx"))

ework<-excellSheet(PCAgenesGG)
saveWorkbook(ework,path.expand("~/Dropbox/UTX-Alex/Paper/TabS10(PCAGG).xlsx"))

ework<-excellSheet(PCAgenesNGG)
saveWorkbook(ework,path.expand("~/Dropbox/UTX-Alex/Paper/TabS12(PCANGG).xlsx"))




## make motif files (on cluster)

for(i in contexts){
    homerWrapper(env$fasta,env$reg[,i],env$reg[,"NONE"],"/data/binaries/homer/bin/homer2",motifsFile=paste0("~/Dropbox/Data/homer_",i,"_PCA_6.txt"),opts="-S 25 -len 6")
    homerWrapper(env$fasta,env$reg[,i],env$reg[,"NONE"],"/data/binaries/homer/bin/homer2",motifsFile=paste0("~/Dropbox/Data/homer_",i,"_PCA_7.txt"),opts="-S 25 -len 7")
    homerWrapper(env$fasta,env$reg[,i],env$reg[,"NONE"],"/data/binaries/homer/bin/homer2",motifsFile=paste0("~/Dropbox/Data/homer_",i,"_PCA_8.txt"),opts="-S 25 -len 8")
}


bed<-list()
bed$bed<-summitToSE(do.call(rbind,lapply(addNames(unifiedPeaks,contexts,list=TRUE),function(x) x[,c(1,2)])))


peakLength<-lapply(unifiedPeaks,function(x) dim(x)[1])
peakLengthSummed<-c(0,Reduce(function(a,b) c(a,a[length(a)]+b),peakLength))

bed$reg<-addColnames(apply(cbind(peakLengthSummed[1:length(peakLengthSummed)-1],peakLengthSummed[2:length(peakLengthSummed)]),1,function(x) makeLogic(seq(x[1]+1,x[2]),peakLengthSummed[length(peakLengthSummed)])),contexts)

which(apply(bed$reg,1,sum)>1)

bed<-addFasta(bed)

for(i in contexts){
    tfast<-append(bed$fasta[bed$reg[,i],],env$fasta[env$reg[,"NONE"],])
    fg<-seq(1:sum(bed$reg[,i]))
    bg<-seq((sum(bed$reg[,i])+1):sum(env$reg[,"NONE"]))
    homerWrapper(tfast,fg,bg,"/data/binaries/homer/bin/homer2",motifsFile=paste0("~/Dropbox/Data/homer_",i,"_BED_6.txt"),opts="-S 25 -len 6")
    homerWrapper(tfast,fg,bg,"/data/binaries/homer/bin/homer2",motifsFile=paste0("~/Dropbox/Data/homer_",i,"_BED_7.txt"),opts="-S 25 -len 7")
    homerWrapper(tfast,fg,bg,"/data/binaries/homer/bin/homer2",motifsFile=paste0("~/Dropbox/Data/homer_",i,"_BED_8.txt"),opts="-S 25 -len 8")
}


### call stamp on motifs
stampDF<-function(combs,type,pvalue=20,sd=2,dir="~/Dropbox/Data/homer-paper/homer_")
    do.call(rbind,apply(combs,1,function(x) data.frame(file=paste(dir,x[1],"_",x[2],".txt",sep=""),compatator=x[3],name=x[1],size=x[2],pvalue=pvalue,sd=sd)))


combs<-stampDF(cbind(strsplit("Erythroid_PCA Erythroid_BED Leukemia_PCA Leukemia_BED HSC_PCA HSC_BED ECFC_PCA ECFC_BED"," ")[[1]],c(rep(6,8),rep(7,8),rep(8,8)),c(rep("TRANSFAC_Fams",24),rep("JASPAR_Fams",24))))

library(httr)

db<-multiStamp(combs)

write.table(db,"~/Dropbox/Data/homer-paper/motifAnnotations.txt",quote=FALSE,row.names=FALSE)

db$ann<-db$genes

# from tall-erythroid-motif.r
sortBy<-function(df,c1){
    df[order(df[[c1]]),]
}

find<-function(value,lin)
    makeLogic(grep(value,lin,ignore.case = TRUE),length(lin))


buildList<-function(x,y){
    do.call(rbind,lapply(x,function(x) do.call(rbind,lapply(y,function(y) cbind(x,y)))))
}

buildFM<-function(treat,control,Msize,cutoff){
    temp<-apply(buildList(c("TRANSFAC_Fams","JASPAR_Fams"),seq(3)),1,function(x) 
    sortBy(db[with(db,{find(treat,dataset) &find(control,dataset)&comparator==x[1]  &rank==x[2] & size==Msize &MACS==cutoff} ),cols],"order"))
cbind(temp[[1]][,c("motif","order","pvalue")],do.call(cbind,lapply(temp,function(x) x[,c("ann","escore")])))
}


file<-path.expand("~/Dropbox/Data/homer-paper/stamp_annot.xlsx")

wb<-createWorkbook()

cols<-c("ann","motif","dataset","order","rank","escore","MACS","pvalue")

for(treat in contexts)
    for(control in c("PCA","BED")){
        sn<-paste0(treat,"_",control)
        ts<-createSheet(wb,sn)
        srow=1
        for(size in c(6,7,8)){            
            bf<-buildFM(treat,control,size,20)            
            addDataFrame(bf,sheet=ts,row.names = FALSE,startRow = srow)
            srow=srow+dim(bf)[1]+2
            print(sn)
            print(srow)
            print(size)
        }
    }

saveWorkbook(wb,file)


## need to go through these by hand and build a subset

## David test
source("http://bioconductor.org/biocLite.R")
biocLite("RDAVIDWebService")

library(rJava)
.jinit()
.jcall("java/lang/System", "S", "getProperty", "java.runtime.version")


library("RDAVIDWebService")
library("BACA")

user<-"agriffith@ohri.ca"

david<-DAVIDWebService(email="agriffith@ohri.ca", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

fileName<-system.file("files/termClusterReport1.tab.tar.gz",
  package="RDAVIDWebService")
untar(fileName)

termCluster<-DAVIDTermCluster(untar(fileName, list=TRUE))
termCluster

head(summary(termCluster))

david<-DAVIDWebService$new(url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

getHttpProtocolVersion(david)
                                        #setHttpProtocolVersion(david, "HTTP/1.0")

setEmail(david,"agriffith@ohri.ca")
show(david)

connect(david)

setTimeOut(david, 50000);

RDAVIDWebService::getGeneListNames(david)

####
david<-DAVIDWebService$new()


source("~/r-workspace/CCCA.r")
source("~/r-workspace/project.r")

source("~/Masters/CCCA/inst/scripts/eboxFrequency.r")

library("RDAVIDWebService")
library("BACA")

user<-"agriffith@ohri.ca"

library(org.Hs.eg.db)

selectGenes<-function(matrix,context){
    df<-data.frame(rownames(matrix)[matrix[,context]==1])
    colnames(df)<-context
    df
}

categories<-as.character(unlist(read.table("~/Dropbox/Data/categories/22-categories.txt")))

contexts<-unique(swapFunD(categories))



load("~/Dropbox/UTX-Alex/Paper/genes.RData")
genelist<-addNames(lapply(contexts,function(x) selectGenes(PCAgenes,x)),contexts,list=TRUE)

symbolToENTREZID<-function(list)
    Filter(function(x)!is.na(x) ,select(org.Hs.eg.db, as.character(unlist(list)), "ENTREZID", "SYMBOL")[,2])

ez<-lapply(genelist,symbolToENTREZID)

#getIdTypes(david)

a<-local({
                                        #david<-DAVIDWebService$new(email="agriffith@ohri.ca")
    
result <- addList(david, ez[[4]], 
                  idType = "ENTREZ_GENE_ID", listName = contexts[4], 
                  listType = "Gene")

setCurrentSpecies(david, 1)
RDAVIDWebService::getSpecieNames(david)
setAnnotationCategories(david, "GOTERM_BP_FAT")
res1 <- getFunctionalAnnotationChart(david, threshold = 1, count = 2)
setAnnotationCategories(david, "GOTERM_MF_FAT")
res2 <- getFunctionalAnnotationChart(david, threshold = 1, count = 2)
list(res1,res2)
})
### composite motif

Erythroid<-list(res1,res2)


PCAGO<-lapply(contexts,function(con)globalenv()[[con]])

names(PCAGO)<-contexts

#PCAGO<-lapply(contexts,function(con)globalenv()[[con]])



excellSheet<-function(matrix, GO){
    ework<-createWorkbook()
    s<-createSheet(ework,"Genes")
    for(i in seq(dim(matrix)[2])){
        addDataFrame(selectGenes(matrix,colnames(matrix)[i]),col.names=TRUE,row.names=FALSE,sheet=s,startColumn = i)
    }
    for(i in colnames(matrix)){
        if(i %in% contexts){
        s<-createSheet(ework,paste0("BP-",i))
        r<-addDataFrame(GO[[i]][[1]],col.names=TRUE,row.names=FALSE,sheet=s)
        print(paste0("MF-",i))
        d<-createSheet(ework,paste0("MF-",i))
        r<-addDataFrame(GO[[i]][[2]],col.names=TRUE,row.names=FALSE,sheet=d)
    }
    }
    ework    
}

tableGene<-function(matrix,file){
for(i in seq(dim(matrix)[2])){
    write.table(selectGenes(matrix,colnames(matrix)[i]),file,col.names=TRUE,row.names=FALSE,quote=FALSE)
}
}

library(xlsx)

ework<-excellSheet(PCAgenes,PCAGO)

saveWorkbook(ework,path.expand("~/Dropbox/UTX-Alex/Paper/test-TabS8(PCA).xlsx"))



genGO<-function(contexts,ez){
david<-DAVIDWebService$new(url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
getHttpProtocolVersion(david)
setEmail(david,"agriffith@ohri.ca")
show(david)
connect(david)
setTimeOut(david, 100000);
RDAVIDWebService::getGeneListNames(david)

getBPMF<-function(context,david){
    print(context)
    result <- addList(david, ez[[context]], 
                      idType = "ENTREZ_GENE_ID", listName = context, 
                      listType = "Gene")

setCurrentSpecies(david, 1)
RDAVIDWebService::getSpecieNames(david)
setAnnotationCategories(david, "GOTERM_BP_FAT")
res1 <- getFunctionalAnnotationChart(david, threshold = 1, count = 2)
setAnnotationCategories(david, "GOTERM_MF_FAT")
res2 <- getFunctionalAnnotationChart(david, threshold = 1, count = 2)
list(res1,res2)
}

addNames(lapply(contexts, getBPMF,david),contexts,list=TRUE)
}


ez<-lapply(addNames(lapply(contexts,function(x) selectGenes(PCAgenesNGG,x)),contexts,list=TRUE),symbolToENTREZID)

PCANGGGO<-genGO(contexts,ez)

ework<-excellSheet(PCAgenesNGG,PCANGGGO)
saveWorkbook(ework,path.expand("~/Dropbox/UTX-Alex/Paper/test-TabS12(PCANGG).xlsx"))


ez<-lapply(addNames(lapply(contexts,function(x) selectGenes(PCAgenesGG,x)),contexts,list=TRUE),symbolToENTREZID)

PCAGGGO<-genGO(contexts,ez)

ework<-excellSheet(PCAgenesGG,PCAGGGO)
saveWorkbook(ework,path.expand("~/Dropbox/UTX-Alex/Paper/test-TabS10(PCAGG).xlsx"))


ez<-lapply(addNames(lapply(contexts,function(x) selectGenes(BEDgenesGG,x)),contexts,list=TRUE),symbolToENTREZID)

BEDGGGO<-genGO(contexts,ez)

ework<-excellSheet(BEDgenesGG,BEDGGGO)
saveWorkbook(ework,path.expand("~/Dropbox/UTX-Alex/Paper/test-TabS9(BEDGG).xlsx"))


ez<-lapply(addNames(lapply(contexts,function(x) selectGenes(BEDgenesNGG,x)),contexts,list=TRUE),symbolToENTREZID)

BEDNGGGO<-genGO(contexts,ez)

save(BEDNGGGO,file="~/Dropbox/UTX-Alex/Paper/BEDNGGGO(too-big).RData")

#ework<-excellSheet(BEDgenesNGG)

#saveWorkbook(ework,path.expand("~/Dropbox/UTX-Alex/Paper/test-TabS11(BEDNGG).xlsx"))

### cluster motifs

# see aug 1st
