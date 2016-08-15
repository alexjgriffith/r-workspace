library(CCCA)

source("~/r-workspace/project.r")

source("~/Masters/CCCA/inst/scripts/eboxFrequency.r")


library("BSgenome.Hsapiens.UCSC.hg19")

env<-getPRC20(2)

env<-addFasta(env)

genEboxCombs()
sum(env$reg[,"Leukemia"])


length(a)

a<-grepMotifs("CAGGTG",env$fasta[env$reg[,"Leukemia"]])

cclocs<-cbind(env$bed[env$reg[,"Leukemia"],][a,],PC1=env$prc$eigenVectors[env$reg[,"Leukemia"],][a,1],env$prc$normData[env$reg[,"Leukemia"],][a,])

write.table(cclocs, "~/Dropbox/Data/cc_locations.tab",quote=FALSE,row.names=FALSE,sep="\t")

env$prc$eigenVectors

ll<-sum(env$reg[,"Leukemia"]==1)

(ll-length(a))/ll


### RNA
rna<-addColnames(read.table("~/Dropbox/Data/rna/TAL1-KD_DESeq2_pvalFC.rnk"),c("gene","fc"))
rnaSle<-data.frame(dr=log(abs(rna$fc),2)>1 & rna$fc>0,
                   ur=log(abs(rna$fc),2)>1 & rna$fc<0)

##CCCA

regions<-makeGeneRegions(
    "/home/agriffith/Dropbox/UTX-Alex/jan/hg19.RefSeqGenes.csv",rna=rna)
geneList<-makeGeneList("/home/agriffith/Dropbox/UTX-Alex/jan/hg19.RefSeqGenes.csv",rna=rna)

library(functional)

ggebox<-makeLogic(which(env$reg[,"Leukemia"])[a],dim(env$reg)[1])
#genesGNN<-geneMatrix(env$over,,dim(env$reg)[1]),regions,geneList)

genes<-greatGeneAssoc(env$over[ggebox,c(1,2,3)],regions,geneList)

gss<-sapply(as.character(unlist(genes$name2)), function(a) rna[rna$gene==a,2])

regions

aChr <- as.character(lapply(regions, "[[", 1))
aStart <- as.numeric(lapply(regions, "[[", 2))
aEnd <- as.numeric(lapply(regions, "[[", 3))
peaks <- apply(env$bed[ggebox, 2:3], 1, mean)

chr<-env$bed[ggebox, 1]
library(within)

peaks[chr=="chr1"] %within% cbind(aStart[aChr=="chr1"],aEnd[aChr=="chr1"])

geneV<-do.call(rbind,lapply(levels(chr),function(ch){
    s<-aStart[aChr==ch]
    e<-aEnd[aChr==ch]
    do.call(rbind,lapply(peaks[chr==ch],function(peak){
        reg<-which(peak >= s & peak <=e)
        if(length(reg)>0){
            g1<-regions[[reg]]
            t1<-as.character(geneList[as.numeric(g1[4:length(g1)]),]$name2)
            gss<-sapply(t1, function(a) rna[rna$gene==a,2])
            data.frame(name=t1[which.max(abs(gss))], value=gss[which.max(abs(gss))])
        }
        else
            data.frame(name="None",value =0)
}))}))

#data.frame(name=lapply(gss,"[[",1),value=lapply(gss,"[[",2))


out<-cbind(env$bed[ggebox,],geneV,PC1=env$prc$eigenVectors[ggebox,"PC1"])

data<-out[order(abs(out$PC1),decreasing=TRUE),]

down<-data[data$value>2,]
up<-data[data$value<(-2),]

write.table(down, "~/Dropbox/Data/cc_locations_down.tab",quote=FALSE,row.names=FALSE,sep="\t")

write.table(up, "~/Dropbox/Data/cc_locations_up.tab",quote=FALSE,row.names=FALSE,sep="\t")

write.table(data, "~/Dropbox/Data/cc_locations.tab",quote=FALSE,row.names=FALSE,sep="\t")

hist(log(out$value))

plot(seq(8,20),sapply(seq(8,20),function(i)length(findMotifDist(env$fasta,"GATAA","CANNTG",env$reg[,"Erythroid"],i+1,i+2,i+3))))

findMotifDist(env$fasta,"GATAA","CAGGTG",env$reg[,"Leukemia"],13,14,15)

gnnebox<-findMotifDist(env$fasta,"GATAA","CACCTG",env$reg[,"Leukemia"],13,14,15)


peaks <- apply(env$bed[gnnebox, 2:3], 1, mean)
chr<-env$bed[gnnebox, 1]
library(within)

peaks[chr=="chr1"] %within% cbind(aStart[aChr=="chr1"],aEnd[aChr=="chr1"])

geneV<-do.call(rbind,lapply(levels(chr),function(ch){
    s<-aStart[aChr==ch]
    e<-aEnd[aChr==ch]
    do.call(rbind,lapply(peaks[chr==ch],function(peak){
        reg<-which(peak >= s & peak <=e)
        if(length(reg)>0){
            g1<-regions[[reg]]
            t1<-as.character(geneList[as.numeric(g1[4:length(g1)]),]$name2)
            gss<-sapply(t1, function(a) rna[rna$gene==a,2])
            data.frame(name=t1[which.max(abs(gss))], value=gss[which.max(abs(gss))])
        }
        else
            data.frame(name="None",value =0)
}))}))

#data.frame(name=lapply(gss,"[[",1),value=lapply(gss,"[[",2))


gnnout<-cbind(env$bed[gnnebox,],geneV,PC1=env$prc$eigenVectors[gnnebox,"PC1"])


write.table(gnnout, "~/Dropbox/Data/gcc_locations.tab",quote=FALSE,row.names=FALSE,sep="\t")
