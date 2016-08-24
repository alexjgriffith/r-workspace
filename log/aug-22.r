library(CCCA)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(org.Hs.eg.db)
library(RDAVIDWebService)
library(xlsx)
library(org.Hs.eg.db)
source(system.file("scripts","eboxFrequency.r",package="CCCA"))
source("~/r-workspace/project/ccca.r")
source("~/r-workspace/project/project.r")
source("~/r-workspace/project/project-variables.r")

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

addMotifs<-function(env,ops,name="motifs"){
    env[[name]]=
        do.call(cbind,lapply(ops,function(o)
        o[[2]](env$fasta)))
    colnames(env[[name]])<-lapply(ops, function(x) x[[1]])
    env
}
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

load("~/Dropbox/UTX-Alex/Paper/Raw Data/peakLocations.RData")
load("~/Dropbox/Data/genomicRegions.RData")


local({
    PCA$over<-PCA$bed
    PCA<-addGenes(PCA,contexts,geneList,PCA$reg)
    PCA$GO<-genGO(user,contexts,PCA$genesEZ)
    makeEx(PCA,"~/Dropbox/PCAGO_GG.xlsx")

})

local({
    PCA$over<-PCA$bed
    PCA<-addGenes(PCA,contexts,geneList,apply(PCA$reg,2,"&",PCA$motifs[,"GG"]))
    PCA$GO<-genGO(user,contexts,PCA$genesEZ)
    makeEx(PCA,"~/Dropbox/PCAGO_GG.xlsx")
})


local({
    PCA$over<-PCA$bed
    PCA<-addGenes(PCA,contexts,geneList,apply(PCA$reg,2,"&",PCA$motifs[,"GNN"]))
    PCA$GO<-genGO(user,contexts,PCA$genesEZ)
    makeEx(PCA,"~/Dropbox/PCAGO_GNN.xlsx")
})

cat(unique(as.character(unlist(selectGenes(PCA$genes,"HSC")))))
cat(unique(as.character(unlist(selectGenes(PCA$genes,"Erythroid")))))
cat(unique(as.character(unlist(selectGenes(PCA$genes,"Leukemia")))))
cat(unique(as.character(unlist(selectGenes(PCA$genes,"ECFC")))))

PCA$GO[,5]


t<-PCA$GO[["Leukemia"]][[1]][PCA$GO[["Leukemia"]][[1]][,5]<0.1,]

Ids<-lapply(t[,6],function(x) as.numeric(strsplit(x,", ")[[1]]))

mat<-outer(ids,ids,Vectorize(function(a,b) {l<-length(intersect(a,b)); ifelse(l>5,l,0)}))

h<-hclust(as.dist(mat))
opts<-cutree(h,h=1)


df<-data.frame(t(combn(t[,2],2)),mat[upper.tri(mat)],t(combn(opts,2)))

r<-df[,3]>0

write.table(data.frame(df[r,1],df[r,4],df[r,2]),file="~/Desktop/test.sif",col.names=FALSE,row.names=FALSE,sep="\t")



length()

na.omit(match(c(1,2,3),c(2,3)))

### Reperesent the number of peaks in each context
rbind(sapply(contexts,function(c) sum(PCA$reg[,c])),
sapply(contexts,function(c) sum(BED$reg[,c])),
sapply(contexts,function(c) sum(UNI$reg[,c])))

