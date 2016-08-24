#source("http://bioconductor.org/biocLite.R")
library(CCCA)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library('org.Hs.eg.db')
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
          makeLogic(locations,length(fasta))}
),
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


load("~/Dropbox/UTX-Alex/Paper/Raw Data/peakLocations.RData")

### Add genes and go GO
library(RDAVIDWebService)
library(xlsx)
library(org.Hs.eg.db)

user<-"agriffith@ohri.ca"


geneFile<-system.file("exdata","hg19.RefSeqGenes.csv",package="CCCA")
geneList<-read.delim(geneFile)
geneDF<-local({
    t<-list()
    t$chr<-as.vector(geneList$chrom)
    t$strand<-geneList$strand
    t$tss<-geneList$txStart
    t$tss[t$strand=="-"]<-geneList$txEnd[t$strand=="-"]
    as.data.frame(t)
})

regions<-CCCA::genomicRegions(geneDF$chr,geneDF$tss,
                        geneDF$strand,1000,5000,100000)

save(regions,file="~/Dropbox/Data/genomicRegions.RData")

load("~/Dropbox/Data/genomicRegions.RData")

genes<-geneMatrix(PCA$over,PCA$reg[,contexts],   
                  regions,geneList,id="name2")

colnames(genes)<-contexts

ez<-symbolToENTREZID(unique(selectGenes(genes,"Leukemia",df=FALSE)))

    david<-DAVIDWebService$new(url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
    getHttpProtocolVersion(david)
    setEmail(david,user)
    show(david)
    connect(david)
    setTimeOut(david, 100000);
RDAVIDWebService::getGeneListNames(david)

result <- addList(david, ez, 
                          idType = "ENTREZ_GENE_ID", listName = "Leukemia", 
                          listType = "Gene")
        
        setCurrentSpecies(david, 1)
        RDAVIDWebService::getSpecieNames(david)
        setAnnotationCategories(david, "GOTERM_BP_FAT")
        res1 <- getFunctionalAnnotationChart(david, threshold = 1, count = 2)

genesEZ<-local({
    fun<-function(x){
        symbolToENTREZID(unique(selectGenes(genes,x,df=FALSE)))
    }
    addNames(lapply(contexts,fun),contexts,list=TRUE)
})




addGenes<-function(env,contexts,geneList,regions){
    genes<-geneMatrix(env$over,env$reg[,contexts],                      
                      regions,geneList,id="name2")
    colnames(genes)<-contexts
    genesEZ<-local({
        fun<-function(x){
            symbolToENTREZID(unique(selectGenes(genes,x,df=FALSE)))
        }
        addNames(lapply(contexts,fun),contexts,list=TRUE)
    })
    env$genes<-genes
    env$genesEZ<-genesEZ
    env
}

makeEx<-function( env,name){
    geneGOExcell<-createWorkbook()
    genesWNames<-lapply(contexts,
                        function(x){ 
                            list(BP=replaceEZID(env$GO[[x]][[1]]),
                                 MF=replaceEZID(env$GO[[x]][[2]]))})
    names(genesWNames)<-contexts
    geneGOExcell<-excellSheet(env$genes,genesWNames)
    saveWorkbook(geneGOExcell,path.expand(name))
}


## GO terms for PCA
PCA$over<-PCA$bed
PCA<-addGenes(PCA,contexts,geneList,apply(PCA$reg,2,"&",PCA$motifs[,"GG"]))

PCA$GO<-genGO(user,contexts,PCA$genesEZ)
makeEx(PCA,"~/Dropbox/PCAGO_GG.xlsx")


## GO terms for UNI
UNI<-local({
UNI$over<-UNI$bed
UNI<-addGenes(UNI,contexts,geneList,regions)
UNI$GO<-genGO(user,contexts,UNI$genesEZ)
makeEx(UNI,"~/Dropbox/UNIGO.xlsx")
UNI
})

## GO terms for BED
BED<-local({
BED$over<-BED$bed
BED<-addGenes(BED,contexts,geneList,regions)
BED$GO<-genGO(user,contexts,BED$genesEZ)
makeEx(BED,"~/Dropbox/BEDGO.xlsx")
BED
})

save(PCA,BED,UNI,file="~/Dropbox/UTX-Alex/Paper/Raw Data/peakLocations.RData")


## Prefered distance plots


load("~/Dropbox/UTX-Alex/Paper/Raw Data/peakLocations.RData")

PCA$fasta[PCA$reg[,"Erythroid"]]
