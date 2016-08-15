library(CCCA)
library(parallel)
library(grid)
library(BSgenome.Hsapiens.UCSC.hg19)

categories<-c("jurk","eryt","cd34","cem_1")

peakFiles<-paste0("peakFiles/",categories,"_combined_mock.bed")
readFiles<-paste0("readFiles/",categories,"_unique_nodupes.bed")

env<-generate(peakFiles,rawFiles,macsCutOff=20,
              cl=20,categories=categories)

env$prc<-pca(env$heights)


swapFun<-createSwapFun("jurk Jurkat eryt Erythroid cd34 CD34 cem_1 CEM")

swapFunB<-createSwapFun("jurk Leukemia eryt Erythrod cd32 HSC cem_1 Leuekemia")

swapFunC<-createSwapFun("Leukemia blue Erythroid red HSC orange")

sepFun<-sepAxis(list(name="Erythroid",pc=1,sd=2,fun=">"),
                list(name="Leukemia",pc=1,sd=-2,fun="<"),
                list(name="HSC",pc=2,sd=2,fun=">"))

center<-(normalize(env$prc$eigenVectors[,1])<inner) &
    (normalize(env$prc$eigenVectors[,1])>-inner) &
        (normalize(env$prc$eigenVectors[,2])>-inner) &
            (normalize(env$prc$eigenVectors[,4])>-inner) &
                (normalize(env$prc$eigenVectors[,4])<inner)

env$reg<-cbind(sepFun(env),center)



dend(x,norm=pass,n=3,linkage="average",colours=c("green","red","blue"))


plotPCs(env$prc$eigenVectors,c(1,2),env$prc$normData,categories)

matr<-pca2Matr(env$prc)
plotPCMat2D(matr,c("PC1","PC2"),categories,swapFun,swapFunB,swapFunC)


p<-stackedContribWrapper(env$heights,env$over,categories,list(1),swapFunB,swapFunC)


env<-addFasta(env)
motifFile<-"test.motif"
motifs<-with(env,{
    homerWrapper(fasta,reg[,"Leukemia"],r[,"Center"],"inst/lib/homer-4.7/bin/homer2",motifFile,opts=paste("-S 25 -len ",len," > /tmp/homerTrash 2>&1 ",sep=""))
})
annot<-stampWrapper(motifFile,"TRANSFAC_Fams","Leukemia_Center")


geneFile<-"hg19.RefSeqGenes.csv"
geneList<-read.delim(geneFile)
chrom<-as.character(geneList$chrom)
tss<-as.numeric(geneList$txStart)
strand<-geneList$strand
levels(strand)<-c(-1,1)
strand<-as.numeric(strand)
regions<-genomicRegions(chrom,tss,strand,1000,5000,1000000)
genes<-with(env,{geneMatrix(over,reg,regions,geneList)})
