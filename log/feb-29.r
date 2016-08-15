source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/jan-functions.r")
source("~/r-workspace/feb-functions.r")
source("~/Masters/CCCA/inst/scipts/eboxFrequency.r")
source("~/r-workspace/preferedDistance.r")

rna<-addColnames(read.table("~/Dropbox/Data/rna/TAL1-KD_DESeq2_pvalFC.rnk"),c("gene","fc"))

rnaSle<-data.frame(dr=log(abs(rna$fc),2)>1 & rna$fc>0,
                   ur=log(abs(rna$fc),2)>1 & rna$fc<0)



env<-getMotifInfo(addFasta(getPRC20(2)),c(genEboxCombs()))


LeukemiaBed<-unionRegs(env,swapFunB,"Leukemia")

    

gnnbed<-findMotifDist(LeukemiaBed,"CAGATG","GATAA",LeukemiaBed,-14,-15,-16)
gnnpca<-findMotifDist(env$Leukemia.alt,"CAGATG","GATAA",env$Leukemia.alt,-14,-15,-16)
ccbed<-intersect(union(grep("CACCTG",env$fasta,ignore.case=TRUE),grep("CAGGTG",env$fasta,ignore.case=TRUE)),which(LeukemiaBed))
ccpca<-intersect(union(grep("CACCTG",env$fasta,ignore.case=TRUE),grep("CAGGTG",env$fasta,ignore.case=TRUE)),which(env$Leukemia.alt))


locGNN<-cbind(makeLogic(gnnbed,dim(env$reg)[1]),makeLogic(gnnpca,dim(env$reg)[1]),makeLogic(ccbed,dim(env$reg)[1]),makeLogic(ccpca,dim(env$reg)[1]))


regions<-makeGeneRegions(rna=rna)

geneList<-makeGeneList(rna=rna)

genesGNN<-geneMatrix(env$over,locGNN,regions,geneList)

ga<-unique(rownames(genesGNN))
fc<-addNames(rna$fc,rna$gene,list=TRUE)
gc<-fc[ga]

selsUp<-lapply(seq(4),function(i) intersect(as.character(rna$gene[rnaSle[,2]]),rownames(genesGNN[genesGNN[,i]==1,])))
selsDown<-lapply(seq(4),function(i) intersect(as.character(rna$gene[rnaSle[,1]]),rownames(genesGNN[genesGNN[,i]==1,])))
allUp<-list(log(abs(rna$fc[rnaSle[,2]]),2))
allDown<-list(log(abs(rna$fc[rnaSle[,2]]),2))

abslog2<-function(x) log(abs(x),2)

aal2<-function(x,gc) lapply(x,function(sel) abslog2(gc[sel]))

upR<-list(aal2(selsUp,gc)[[2]],aal2(selsUp,gc)[[4]],unlist(allUp))

downR<-list(aal2(selsDown,gc)[[2]],aal2(selsDown,gc)[[4]],unlist(allDown))


mapply(function(x,name){
    print(name)
    write.table(x, name,col.names = FALSE,row.names=FALSE,quote=FALSE)
},selsUp,paste0("~/Dropbox/UTX-Alex/br-data/genes/",c("bed_gnn","pca_gnn","bed_cc","pca_cc"),"_up.txt"))
mapply(function(x,name){
    print(name)
    write.table(x, name,col.names = FALSE,row.names=FALSE,quote=FALSE)
},selsDown,paste0("~/Dropbox/UTX-Alex/br-data/genes/",c("bed_gnn","pca_gnn","bed_cc","pca_cc"),"_down.txt"))



#png("~/Dropbox/UTX-Alex/br-data/downR.png")
boxplot(downR)
#dev.off()
#png("~/Dropbox/UTX-Alex/br-data/upR.png")
boxplot(upR)
#dev.off()

outer(upR,upR,Vectorize(function(x,y) t.test(x,y)$p.value))

outer(downR,downR,Vectorize(function(x,y) t.test(x,y)$p.value))

