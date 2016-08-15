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
env<-addFasta(env20)

gd<-addColnames(read.table("~/Dropbox/Data/rna/TAL1-KD_DESeq2_pval.rnk"),c("genes","diff"))
gdf<-gd[is.na(gd$genes)==FALSE,]



genes<-geneMatrix(env$over,env$reg,regions,geneList)

a<-sapply(rownames(genes)[genes[,"Erythroid"]==1],function(x) gdf$diff[x==gdf$genes])

a<-Filter(length,sapply(rownames(genes)[genes[,"Leukemia"]==1],function(x) gdf$diff[x==gdf$genes]))



gata<-slowMotLoc("GATAA",compliment("GATAA"),env$fasta,env$reg[,"Leukemia"])

canntg<-slowMotLoc(IUPACtoBase("CANNTG"),compliment(IUPACtoBase("CANNTG")),env$fasta,env$reg[,"Leukemia"])

cacctg<-slowMotLoc("CACCTG",compliment("CACCTG"),env$fasta,env$reg[,"Leukemia"])

length(intersect(cacctg,gata))

length(intersect(canntg,gata))

mList<-c(fox,ebox,"GATAA",nr)
motifs<-c(genEboxCombs(),mList)
cList<-sapply(motifs,compliment)
locationsM<-addNames(lapply(motifs,grep,env$fasta),motifs,list=TRUE)
locationsC<-addNames(lapply(cList,grep,env$fasta),motifs,list=TRUE)


fun<-function(env)
    function(x,y) {
    length(findSharedRegions(locationsM,locationsC,x,"GATAA",env$reg[,y]))/length(intersect(locationsM[[x]],which(env$reg[,y])))
}

tb<-addNames(outer(seq(11),colnames(env$reg),Vectorize(fun(env))),colnames(env$reg),motifs[seq(11)])

reg<-addColnames(do.call(cbind,lapply(1:12,function(x) makeLogic(findSharedRegions(locationsM,locationsC,x,"GATAA") ,length(env$fasta)))),motifs[1:12])




mdToMatrix<-function(md,env){
    nf<-function(data){
        ldata<-which(env$reg[,data])
        mapply(function(x,y){
            length(intersect(x,ldata))/
                length(ldata)
        },md20,reg)
    }
    list<-lapply(colnames(env$reg),nf)
    matrix<-do.call(cbind,list)
    addColnames(matrix ,colnames(env$reg))
}

distanceBetweenMotifs<-function(env,motifs,locationsM,locationsC,m1,m2){
    reg<-addNames(lapply(m1,function(x) findSharedRegions(locationsM,locationsC,x,m2)),sapply(motifs[m1],consenusIUPAC),list=TRUE)         
    md20<-addNames(motifDistances(reg,env$fasta),paste0(names(reg)),list=TRUE)
    md20
}

md20<-distanceBetweenMotifs(env,motifs,locationsM ,locationsC,1:10,"GATAA")
    
mdToMatrix(md20,env)



## note: reason for different sizes between methods
## motifHist includes both + and - strands with all options (its good to see if there are some distant enrichments)
## md20 requires the closest possible option out of + and - such
## that only one occures on each string
## findSharedRegions uses the union region of +/-
plotFastMot(env,motifs,cList,locationsM,locationsC,12,13)

plotFastMot(env,motifs,cList,locationsM,locationsC,6,13)

plotFastMot(env,motifs,cList,locationsM,locationsC,5,13)

length(motifHist(env$fasta,motifs,cList,locationsM,locationsC,6,13,env$reg[,"Leukemia"]))

length(motifHist(env$fasta,motifs,cList,locationsM,locationsC,8,13,env$reg[,"Leukemia"]))

length(intersect(md20[["CACCTG"]],which(env$Leukemia.alt)))

length(intersect(md20[["CAGATG"]],which(env$Leukemia.alt)))

a<-addNames(lapply(1:11,function(x) findSharedRegions(locationsM,locationsC,x,"GATAA")),sapply(motifs[1:11],consenusIUPAC),list=TRUE)

length(intersect(a[["CACCTG"]],which(env$Leukemia.alt)))

length(intersect(a[["CAGATG"]],which(env$Leukemia.alt)))



a<-selectRegionsByDistance(IUPACtoBase("CANNTG"),"GATAA", intersect(reg$CACCTG,which(env$Erythroid.alt)),env$fasta)

stem(min(a):max(a),getHeights(a),xlim(-32,32))



selectRegionsByDistance()
motif1="CACCTG"
motif2="GATAA"
locations=reg$CACCTG
fasta=env$fasta


