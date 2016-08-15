# these are for project

library(CCCA)
library(parallel)
library(grid)
library(BSgenome.Hsapiens.UCSC.hg19)

categories<-CCCA::readCategories("~/Dropbox/UTX-Alex/jan/catagories")
controls<-CCCA::readCategories("~/Dropbox/UTX-Alex/jan/contcats")
categoriesb<-c("jurk","eryt","cd34","cem_1")

makePeakFiles<-function(categories,mock="combined")
     paste0("~/Dropbox/Data/august_peaks/",categories,
           "~",mock,"_mock_peaks.xls")
makeRDFiles<-function(categories)
     paste0("/mnt/brand01-00/mbrand_analysis/data_sets/",categories,"/",categories,
           "_unique_nodupes.bed")

genReg<-function(env,sd,inner){
    ecfcfun<-function(pcs,dis){        
        andM(mapply(function(a,b) {
            pcs[,a]<b },
            seq(dim(pcs)[2]),dis))
    }
    reg<-(sepAxis(list(name="Erythroid",pc=1,sd=sd,fun=">"),
             list(name="Leukemia",pc=1,sd=-sd,fun="<"),
             list(name="HSC",pc=4,sd=sd,fun=">"),
             list(name="ECFC",pc=c(2,4),sd=c(-sd,-sd),fun=ecfcfun),
             inclusive=FALSE))(env)
    regn<-(CCCA:::normalize(env$prc$eigenVectors[,1])<inner) &
    (CCCA:::normalize(env$prc$eigenVectors[,1])>-inner) &
        (CCCA:::normalize(env$prc$eigenVectors[,2])>-inner) &
            (CCCA:::normalize(env$prc$eigenVectors[,4])>-inner) &
                (CCCA:::normalize(env$prc$eigenVectors[,4])<inner)
    cbind(reg,NONE=regn,ALL=rep(TRUE,dim(reg)[1]))    
}

## need to add all the swapfuns from project.r

pf<-makePeakFiles(categories)
rawData<-makeRDFiles(controls)
env<-generate(pf,rawData,macsCutOff=20,cl=20,categories=categories)

env<-addFasta(env)

env<-getMotifInfo(env,c(genEboxCombs(),"CANNTG","GATAA"))

#saveEnv(env,"~/ccca-test")



env<-list(categories=readCategories("~/Dropbox/Data/categories/22-categories.txt"),
     heights=read.table("~/Dropbox/Data/UDM/22_treatment_pvalue_20_control_combined.txt",header=T),
     over=read.table("~/Dropbox/Data/AFS/22_pvalue_20_control_combined.txt",header=T),
     bed=read.table("~/Dropbox/Data/AFS/22_pvalue_20_control_combined.txt",header=T)[,1:3]
     )

env$prc<-pca(env$heights)

env$reg<-genReg(env,2,0.25)

env<-addFasta(env)

env<-getMotifInfo(env,c(genEboxCombs(),"CANNTG","GATAA"))

ds<-list(composite=findMotifDist(env$fasta,"CANNTG","GATAA",
             env$reg[,"Leukemia"],-14,-15,-16),
         gnn=intersectN(grepMotifs("CANNTG",env$fasta),
             grepMotifs("GATAA",env$fasta),
             which(env$reg[,"Leukemia"])),
         cc=intersect(grepMotifs("CACCTG",env$fasta),
             which(env$reg[,"Leukemia"])),
         ga=intersect(grepMotifs("CAGATG",env$fasta),
             which(env$reg[,"Leukemia"])),
         nn=intersect(grepMotifs("CANNTG",env$fasta),
             which(env$reg[,"Leukemia"])),
         runx=intersect(grepMotifs("ACCACA",env$fasta),
             which(env$reg[,"Leukemia"])),
         leukemic=which(env$reg[,"Leukemia"]))

locM<-do.call(cbind,lapply(ds,makeLogic,length(env$fasta)))

rna<-addColnames(read.table("~/Dropbox/Data/rna/TAL1-KD_DESeq2_pvalFC.rnk"),c("gene","fc"))
rnaSle<-data.frame(dr=log(abs(rna$fc),2)>1 & rna$fc>0,
                   ur=log(abs(rna$fc),2)>1 & rna$fc<0)



makeGeneRegions<-function(
    geneFile,
    rna=NULL){
    geneList<-makeGeneList(geneFile,rna)
    chrom<-as.character(geneList$chrom)
    tss<-as.numeric(geneList$txStart)
    strand<-geneList$strand
    levels(strand)<-c(-1,1)
    strand<-as.numeric(strand)
    genomicRegions(chrom,tss,strand,1000,5000,1000000)
}

makeGeneList<-function(geneFile,rna=NULL){
    onlyOne<-function(geneNames){
        geneOrder<-order(geneNames)
        rorder<-order(geneOrder)
        as<-geneNames[geneOrder[1:length(geneOrder)-1]]
        bs<-geneNames[geneOrder[2:length(geneOrder)]]
        mult<-c(TRUE,as!=bs)
        mult[rorder]
    }
    geneListP<-read.delim(geneFile)
    if(is.null(rna))
        rna<-data.frame(gene=geneListP$name2,fc=0)
    genesInList<-sapply(geneListP$name2,
                        function(x)
                            x %in% rna$gene)&onlyOne(geneListP$name2)
    geneListP[genesInList,]
}

geneMatrix<-function(over,reg,regions,geneList,id="name2"){
    geneCount<-function(y,allGenes){
        d<-setdiff(allGenes,y)
        t<-rbind(data.frame(names=y,count=rep(1,length(y))),data.frame(names=d,count=rep(0,length(d))))
        t[order(as.character(t[,1])),2]
    }
    
    genes<-lapply(seq(dim(reg)[2]),function(x) greatGeneAssoc(over[reg[,x],c(1,2,3)],regions,geneList)[,id])
    allGenes<-sort(unique(as.character(unlist(genes))))
    geneCountWrapper<-Compose(as.character,unique,function(x)geneCount(x,allGenes))
    geneMatrix<-addNames(do.call(cbind,lapply(genes,geneCountWrapper)),colnames(reg),allGenes)
    return(geneMatrix)
}

regions<-makeGeneRegions(
    "/home/agriffith/Dropbox/UTX-Alex/jan/hg19.RefSeqGenes.csv",rna=rna)
geneList<-makeGeneList("/home/agriffith/Dropbox/UTX-Alex/jan/hg19.RefSeqGenes.csv",rna=rna)

genesGNN<-geneMatrix(env$over,locM,regions,geneList)

ga<-unique(rownames(genesGNN))
fc<-addNames(rna$fc,rna$gene,list=TRUE)

gcUp<-addNames(lapply(colnames(genesGNN),function(name)Filter(function(x){!is.na(x)},log(abs(fc[rnaSle$ur][rownames(genesGNN[genesGNN[,name]==1,])]),2))),colnames(genesGNN),list=TRUE)

gcDown<-addNames(lapply(colnames(genesGNN),function(name)Filter(function(x){!is.na(x)},log(abs(fc[rnaSle$dr][rownames(genesGNN[genesGNN[,name]==1,])]),2))),colnames(genesGNN),list=TRUE)

png("~/Dropbox/dr-genes.png")
boxplot(gcDown)
dev.off()

png("~/Dropbox/ur-genes.png")
boxplot(gcUp)
dev.off()


p1<-plotPCMat2D(pca2Matr(env$prc),c("PC1","PC2"),controls,Compose(mockToNorm,swapFun),Compose(mockToNorm,swapFunD),swapFunC)
p2<-plotPCMat2D(pca2Matr(env$prc),c("PC1","PC4"),controls,Compose(mockToNorm,swapFun),Compose(mockToNorm,swapFunD),swapFunC)

png("~/Dropbox/controls.png",width=480,height=240)
multiPlot(list(p1,p2),1,2)
dev.off()


gataPeaks<-addColnames(read.table("~/Dropbox/UTX-Alex/br-data/jurk-gata/combined~jurk_gata3_peaks.bed"),c("chr","start","end"))

egata<-addFasta(list(bed=gataPeaks))

length(env)-length(grepMotifs("CANNTG",egata$fasta))

bedOverlaps<-function(a,b,l=700){
    pu<-sortDataFrame(rbind(cbind(a,name="a",summit=(a$end+a$start)/2),
                        cbind(b,name="b",summit=(b$end+b$start)/2)),"chr","summit")
    afs<-unifyBedFile(pu,l)
#    afs[afs$a==afs$b,]
    afs
}

do.call(rbind, lapply(colnames(env$reg),function(name){
    com<-bedOverlaps(egata$bed,env$bed[env$reg[,name],])
    a<-com[(com[,4]==1)&(com[,3]==1),]
    tenv<-addFasta(list(bed=addNames(list(a[,1],a[,2]-350,a[,2]+350),c("chr","start","end"),list=TRUE)))
    sapply(CCCA::genEboxCombs(),function(i)
    findMotifDist(tenv$fasta,i,"GATAA",
             rep(TRUE,length(tenv$fasta)),-14,-15,-16))
    #dim(a)[1]-length(grepMotifs("CANNTG",tenv$fasta))
}))


lapply( list(function(x){length(which(x[,3]==1))},function(x){length(which(x[,4]==1))},function(x){length(which((x[,3]==1)&(x[,4]==1)))}),function(x) x(com))

c(7322-687,1449-687)

dim(env$bed)

a<-sepAxis(list(list(name="a",pc=1,fun=">",sd=2),list(name="b",pc=1,fun="<",sd=-2)))(env)

apply(a,1,sum)==1

uniqueSepAxis(a)



local({
t<-cbind(c(TRUE,TRUE,TRUE),c(FALSE,FALSE,TRUE))
apply(t,1,function(x) sum(x))
dim(t)
})


#

local({single<-generate(makePeakFiles(categories,"single"),
                 makeRDFiles(categories),
                 macsCutOff=5,cl=20,categories=categories)
no<-generate(makePeakFiles(categories,"no"),
                 makeRDFiles(categories),
                 macsCutOff=5,cl=20,categories=categories)
save(single,no,file="~/Dropbox/changebg.RData")})

library("functional")
library("grid")
library("ggplot2")

p1<-plotPCMat2D(pca2Matr(single$prc),c("PC1","PC2"),categories,swapFun,swapFunD,swapFunC)
p2<-plotPCMat2D(pca2Matr(single$prc),c("PC1","PC3"),categories,swapFun,swapFunD,swapFunC)
p3<-plotPCMat2D(pca2Matr(single$prc),c("PC3","PC5"),categories,swapFun,swapFunD,swapFunC)
p4<-plotPCMat2D(pca2Matr(single$prc),c("PC3","PC7"),categories,swapFun,swapFunD,swapFunC)
png("~/Dropbox/single-s.png",width=480,height=480)
multiPlot(list(p1,p2,p3,p4),2,2)
dev.off()
p1<-plotPCMat2D(pca2Matr(no$prc),c("PC1","PC2"),categories,swapFun,swapFunD,swapFunC)
p2<-plotPCMat2D(pca2Matr(no$prc),c("PC1","PC3"),categories,swapFun,swapFunD,swapFunC)
p3<-plotPCMat2D(pca2Matr(no$prc),c("PC3","PC6"),categories,swapFun,swapFunD,swapFunC)
p4<-plotPCMat2D(pca2Matr(no$prc),c("PC3","PC7"),categories,swapFun,swapFunD,swapFunC)
png("~/Dropbox/no-s.png",width=480,height=480)
multiPlot(list(p1,p2,p3,p4),2,2)
dev.off()


#' PCA tp Matrix Transformation
#' 
#' @param x the normalized input data and eigenvectors 
#' @param n normalizing function applied to the eigevectors
#' @return the dot product of the normalized data and eigenvectors
#' @export
pca2Matr<-function(x,n=pass){
    if(is.null(x$normData) | is.null(x$eigenVectors))
        stop("x needs variables normData and eigenVectors")
    # need to add tests for dimensions
    t(x$normData)%*%apply(x$eigenVectors,2,n)
}


plotPCMat2D<-function(matr,pcs,categories,swapFun,swapFunB,swapFunC,...){
        df<-data.frame(x=matr[,pcs[1]],y=matr[,pcs[2]],categories=swapFun(categories),Conditions=swapFunB(categories))
        plotPCMatAux(df,pcs,categories,swapFunC(unique(sort(swapFunB(categories)))),...)
    }


plotPCMatAux<-function(df,pcs,categories,colours,sdf=NULL,text=NULL,legend=NULL){
    if(is.null(sdf))
       postext<-df
   else
       postext<-shiftCols(df$x,df$y,categories,sdf)
    p<-ggplot(df,aes(x=x,y=y,col=Conditions,label=categories))+ylab(pcs[2])+xlab(pcs[1])+geom_point(size=6,shape=20)+theme_bw()
    if(! is.null(text))
        p<-p+geom_text(x=postext$x,y=postext$y,show_guide=F,size=5)
    if(! is.null(colours))
        p<-p+scale_color_manual(values=colours)
    if(is.null(legend))
        p<-p+theme(legend.position="none")+scale_x_continuous(breaks=NULL)+scale_y_continuous(breaks=NULL)
    p
}



plotPCMat2D(pca2Matr(env$prc),c("PC1","PC2"),categories,swapFun,swapFunB,swapFunC)



lapply(colnames(env$reg),function(x) write.table(env$fasta[env$reg[,x]],paste0("~/Dropbox/",x,".fasta"),quote=FALSE,col.names=FALSE,row.names=FALSE))
