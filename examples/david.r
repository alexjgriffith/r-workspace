## Install david web service

source("http://bioconductor.org/biocLite.R")
biocLite("RDAVIDWebService")

install.packages('rJava',,'http://www.rforge.net/',type="source")
library(rJava)
.jinit()
.jcall("java/lang/System", "S", "getProperty", "java.runtime.version")
## must be jre 1.8 if not update your system



library("RDAVIDWebService")
library(org.Hs.eg.db) # used to transform GENE NAMES to EZ IDS
library(CCCA)
source("~/r-workspace/project/project.r")
source("~/r-workspace/project/project-variables.r")
source("~/r-workspace/project/ccca.r")


# check if you can connect to david
david<-DAVIDWebService$new(url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
getHttpProtocolVersion(david)
setEmail(david,"agriffith@ohri.ca")
show(david)
connect(david)
setTimeOut(david, 100000);
RDAVIDWebService::getGeneListNames(david)

env<-getPRC20(2)
geneFile<-"~/Dropbox/UTX-Alex/jan/hg19.RefSeqGenes.csv"
geneList<-read.delim(geneFile)





t<-list()
t$chr<-as.vector(geneList$chrom)
t$strand<-geneList$strand
t$tss<-geneList$txStart
t$tss[t$strand=="-"]<-geneList$txEnd[t$strand=="-"]
t<-as.data.frame(t)
    
regionsB<-genomicRegions(t$chr,t$tss,t$strand,1000,5000,1000000)


strand<-t$strand
levels(strand)<-c(-1,1)

regions<-CCCA::genomicRegions(t$chr,t$tss,t$strand,1000,5000,100000)


eryt<-greatGeneAssoc(env$over[CCCA::normalize(env$prc$eigenVectors[,1])>2&CCCA::normalize(env$prc$eigenVectors[,2])>-2,c(1,2,3)],regions,geneList)

eryt<-greatGeneAssoc(env$over[env$reg[,"Erythroid"],c(1,2,3)],regionsB,geneList)

CCCA::genomicRanges(t$chr,t$tss,t$strand,1000,5000,1000000)

genes<-geneMatrix(env$over,env$reg[,contexts],regions,geneList,id="name2")

## Functions
symbolToENTREZID<-function(list)
    Filter(function(x)!is.na(x) ,select(org.Hs.eg.db, as.character(unlist(list)), "ENTREZID", "SYMBOL")[,2])

genesEZ<-addNames(lapply(contexts,function(x)symbolToENTREZID(unique(selectGenes(genes,x,df=FALSE)))),contexts,list=TRUE)

davidResults<-genGO(contexts,genesEZ)



davidResults[["Erythroid"]][[1]][,5]<(1e-1)

ezidGene<-function(ezid){
    ezidSplit<-strsplit(ezid,", ")[[1]]
    geneName<-select(org.Hs.eg.db,ezidSplit, "SYMBOL" ,"ENTREZID")[,2]
    paste(geneName,collapse=", ")
}

replaceEZID<-function(table,pvalue=1e-1){
    sigReg<-table[,5]<pvalue
    geneNames<-sapply(table[sigReg,6],ezidGene)
    stable<-table[sigReg,]
    stable[,6]<-geneNames
    stable
}


genesWNames<-addNames(lapply(contexts,
                             function(x) list(BP=replaceEZID(davidResults[[x]][[1]]),MF=replaceEZID(davidResults[[x]][[2]]))),contexts ,list=TRUE)
       
BPER<-replaceEZID(davidResults[["Erythroid"]][[1]])


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

library(xlsx)
ework<-createWorkbook()
ework<-excellSheet(genes,genesWNames)
saveWorkbook(ework,path.expand("~/Dropbox/UTX-Alex/Paper/fixed-TabS12(PCA).xlsx"))



lapply(davidResults[["Erythroid"]][[1]][1:10,6],ezidGene)
    
genGO<-function(user,contexts,ez){
    david<-DAVIDWebService$new(url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
    getHttpProtocolVersion(david)
    setEmail(david,user)
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
    #cs<-makeForkCluster(4)
    addNames(lapply(contexts, getBPMF,david),contexts,list=TRUE)
    #stopCluster(cs)
}

