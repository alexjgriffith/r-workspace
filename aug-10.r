#options(max.print=200)
library(CCCA)
source("~/r-workspace/project.r")
source("~/r-workspace/project-variables.r")
source("~/r-workspace/ccca.r")


## Test if geneRegions are roughly equivelent
## There are a few more regions found with the new
## Method, need to check that out
env<-getPRC20(2)
list1<-genGeneTSS()
geneFile<-"~/Dropbox/UTX-Alex/jan/hg19.RefSeqGenes.csv"
list2<-read.delim(geneFile)
## Make sure that there the gene sets are roughly equivelent
length(intersect(as.character(list1$name),as.character(list2$name2)))/length(union(as.character(list1$name),as.character(list2$name2)))
## Fix strandedness naming convention for old method
strand<-list1$strand
levels(strand)<-c(-1,1)
strand<-as.numeric(strand)
## compare the defunct genomicRegions in CCCA with the new faster version
regions1<-CCCA::genomicRegions(list1$chr,list1$tss,strand,1000,5000,100000)
regions2<-genomicRegions(list1$chr,list1$tss,list1$strand,1000,5000,100000)
## turn the list of regions into data frames
treg1<-addColnames(as.data.frame(do.call(rbind,lapply(regions1,function(x) c(x[c(1,2,3)],length(x)-3)))),strsplit("chr start end num"," ")[[1]])
treg2<-addColnames(as.data.frame(do.call(rbind,lapply(regions2,function(x) c(x[c(1,2,3)],length(x)-3)))),strsplit("chr start end num"," ")[[1]])
## sort and visualize
treg1[with(treg1,order(chr,as.numeric(as.vector(start))))[1:10],]
treg2[with(treg2,order(chr,as.numeric(as.vector(start))))[1:10],]


## Rerun David using new geneRegions function
env<-getPRC20(2)

glist<-genGeneTSS()

glistr<-glist[with(glist,order(chr,tss)),]

regions<-with(glist,genomicRegions(chr,tss,strands,1000,5000,100000))

genes<-geneMatrix(env$over,env$reg[,contexts],regions,glist,"id")


library("RDAVIDWebService")
david<-DAVIDWebService$new(url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
getHttpProtocolVersion(david)
setEmail(david,"agriffith@ohri.ca")
show(david)
connect(david)
setTimeOut(david, 50000);
RDAVIDWebService::getGeneListNames(david)

symbolToENTREZID<-function(list)
    Filter(function(x)!is.na(x) ,select(org.Hs.eg.db, as.character(unlist(list)), "ENTREZID", "SYMBOL")[,2])


myDavid<-function(david,genelist,name){
    result <- addList(david, genelist, 
                  idType = "ENTREZ_GENE_ID", listName = name, 
                  listType = "Gene")
    RDAVIDWebService::getSpecieNames(david)
    setAnnotationCategories(david, "GOTERM_BP_FAT")
    res1 <- getFunctionalAnnotationChart(david, threshold = 1, count = 2)
    setAnnotationCategories(david, "GOTERM_MF_FAT")
    res2 <- getFunctionalAnnotationChart(david, threshold = 1, count = 2)
    list(res1,res2)
}

PCAgenesNew<-list()

cat<-"Leukemia"
PCAgenesNew[[cat]]<-myDavid(david,symbolToENTREZID(list(selectGenes(genesNAMEC,cat,FALSE))),cat)

cat<-"Erythroid"
PCAgenesNew[[cat]]<-myDavid(david,symbolToENTREZID(list(selectGenes(genesNAMEC,cat,FALSE))),cat)


cat<-"HSC"
PCAgenesNew[[cat]]<-myDavid(david,symbolToENTREZID(list(selectGenes(genesNAMEC,cat,FALSE))),cat)


cat<-"ECFC"
PCAgenesNew[[cat]]<-myDavid(david,symbolToENTREZID(list(selectGenes(genesNAMEC,cat,FALSE))),cat)

PCAGONEW<-PCAgenesNew

genesPCANEW<-genesNAMEC

save(genesPCANEW,PCAGONEW,file="~/Desktop/PCAGONEW.RData")




## My personal GO approach
library(compiler)
innerO<-cmpfun(Vectorize(function(a,b){
    y<-.Internal(pmatch(as.character(a), as.character(b), 0, TRUE))
    n<-.Internal(unique(y, FALSE, FALSE, min(length(y), length(attr(y,"levels")) + 1L)))
    la<-length(.Internal(unique(a, FALSE, FALSE, min(length(a), length(attr(a,"levels")) + 1L))))
    lb<-length(.Internal(unique(b, FALSE, FALSE, min(length(b), length(attr(a,"levels")) + 1L))))
    mm<-min(la,lb)
    rem<-length(.Internal(which(n==0)))
    1==(length(n)-rem)/mm#/min(length(unique(a)),length(unique(b)))
    }),options=list(optimize=3,supressALL=TRUE))


genesNAME<-geneMatrix(env$over,env$reg[,contexts],regions,glist,"name")

genesNAMEB<-geneMatrix(env$over,env$reg[,contexts],regions2,list2,"name2")

x<-as.list(org.Hs.egGO2EG)
rm<-as.list(org.Hs.egGO[mappedkeys(org.Hs.egGO)])
numGenes<-length(rm)
names(rm)<-select(org.Hs.eg.db,names(rm), "SYMBOL","ENTREZID")$SYMBOL
gosize<-sapply(as.list(x),length)

agenes<-selectGenes(genesNAME,"Erythroid",TRUE)

lgenes<-length(agenes)
gt<-as.data.frame(do.call(rbind,lapply(agenes,function(i)
    do.call(rbind,lapply(rm[[i]],function(x)cbind(GOID=x[["GOID"]],Ontology=x[["Ontology"]],Evidence=x[["Evidence"]],gene=i))))))
##gt$Evidence=="TAS"

temp<-sapply(as.character(unique(gt[gt$Ontology=="BP","GOID"])),function(sec)gt[gt$GOID==sec,"gene"])
temp<-Filter(function(x) length(x)>2,temp)
subt<-names(temp)
temp2<-as.data.frame(t(sapply(subt,function(sec) c(length(gt[gt$GOID==sec,"gene"]),gosize[sec],numGenes,lgenes))))
p<-dbinom(temp2[,1],temp2[,4],temp2[,2]/temp2[,3])
l<-sum(sort(p<0.05))

cn<-names(temp)[order(p)[1:l]]
temp<-lapply(order(p)[1:l],function(x) temp[[x]])
names(temp )<-cn
    


res<-outer(temp,temp,innerO)
subt<-unique(apply(res,1,function(x) names(sort(sapply(which(x),function(i)length(temp[[i]])),decreasing=FALSE))[1]))



temp2<-as.data.frame(t(sapply(subt,function(sec) c(length(gt[gt$GOID==sec,"gene"]),gosize[sec],numGenes,lgenes))))
p<-dbinom(temp2[,1],temp2[,4],temp2[,2]/temp2[,3])
l<-sum(sort(p<0.01))
nms<-do.call(rbind,lapply(rownames(temp2[order(p)[1:l],]),function(x) paste(temp[[x]],collapse=",")))
terse<-data.frame(GOID=rownames(temp2[order(p)[1:l],]),
                  term=select(GO.db,rownames(temp2[order(p)[1:l],]),"TERM","GOID")[,2],
           p=p[order(p)[1:l]],genes=nms,
           addColnames(temp2[rownames(temp2[order(p)[1:l],]),],c("a","b","c","d")))
terse[,2]

## check pcagenes with geneNAMES
outer(contexts,contexts,Vectorize(function(a,b) do.call(function(...)length(intersect(...))/length(union(...)),list(selectGenes(genesNAMEB,a,FALSE),selectGenes(genesNAMEC,b,FALSE)))))

at<-lapply(list(PCAgenes,genesNAME),selectGenes,"Erythroid",FALSE)


do.call(rbind,lapply(as.character(glistr$name[1:20]),function(x){t<-list2[which(x==list2$name2),c("chrom","txEnd","strand","name2")];t[order(t$txEnd),];t[1,]}))$txEnd - glistr$tss[1:20]

neg<-list2$strand=="-"
pos<-list2$strand=="+"
out<-rep(FALSE,length(neg))
out[neg]<-list2$txEnd[neg]
out[pos]<-list2$txStart[pos]

strand<-list2$strand
levels(strand)<-c(-1,1)
strand<-as.numeric(strand)

regions2<-genomicRegions(list2$chrom,as.numeric(as.vector(out)),strand,1000,5000,100000)

genesNAMEC<-geneMatrix(env$over,env$reg[,contexts],regions2,list2,"name2")
