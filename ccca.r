# Functions to be added to CCCA
library(CCCA)
library(functional)
library(ggplot2)
library(parallel)
library(grid)
#library(BSgenome.Hsapiens.UCSC.hg19)

readPeaksXLS <- function(file, name = file,pValue=20) {
    bedData <- read.table(file, header = TRUE, skip = "#")
    bedData <- bedData[bedData$X.log10.pvalue > pValue, ]
    if (dim(bedData)[1] > 0) {
        ret <- data.frame(bedData$chr, bedData$abs_summit, 
                          name)
        colnames(ret) <- c("chr", "summit", "name")
    }
    else ret <- NULL
    ret
}


mergeFun<-function(ma,swapFun){
    newCols<-swapFun(colnames(ma))
    #print(newCols)
    unc<-unique(newCols)
    outL<-c()
    for(i in unc){
        pos<-which(i==newCols)
        #print(pos)
        #print(data.frame(c=colnames(ma)[pos],v=do.call(rbind,lapply(pos,function(x) sum(ma[,x])))))
        if(length(pos)>1)
            temp<-orM(ma[,pos])
        else{
            temp<-ma[,pos]
            print(sum(temp))
        }
        outL<-cbind(outL,temp)}
    
    colnames(outL)<-unc
    outL
}

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

bedOverlaps<-function(a,b,l=700){
    pu<-sortDataFrame(rbind(cbind(a,name="a",summit=(a$end+a$start)/2),
                        cbind(b,name="b",summit=(b$end+b$start)/2)),"chr","summit")
    afs<-unifyBedFile(pu,l)
#    afs[afs$a==afs$b,]
    afs
}


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
    print(dim(x$normData))
    print(dim(x$eigenVectors))
    t(x$normData)%*%apply(x$eigenVectors,2,n)
}


plotPCMat2D<-function(matr,pcs,categories,swapFun,swapFunB,swapFunC,...){
    print(swapFun(categories))
        df<-data.frame(x=matr[,pcs[1]],y=matr[,pcs[2]],categories=swapFun(categories),Conditions=swapFunB(categories))
        plotPCMatAux(df,pcs,categories,swapFunC(unique(sort(swapFunB(categories)))),...)
    }


plotPCMatAux<-function(df,pcs,categories,colours,sdf=NULL,text=NULL,legend=NULL,label=NULL,blank=NULL){
    if(is.null(sdf))
       postext<-df
   else
       postext<-shiftCols(df$x,df$y,categories,sdf)
    p<-ggplot(df,aes(x=x,y=y,col=Conditions,label=categories))+geom_point(size=10,shape=20)+theme_bw()+ylab(pcs[2])+xlab(pcs[1])
    if(! is.null(blank))
    p<-p+ theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
    if(! is.null(label))
        p<-p+ylab("")+xlab("")
    if(! is.null(text))
        p<-p+geom_text(x=postext$x,y=postext$y,show_guide=F,size=5)
    if(! is.null(colours))
        p<-p+scale_color_manual(values=colours)
    if(is.null(legend))
        p<-p+theme(legend.position="none")+scale_x_continuous(breaks=NULL)+scale_y_continuous(breaks=NULL)
    p
}

# need to give proper variable names
removeDataSets<-function(x,env){
    a<-apply(cbind(env$over[,c(x)],0),1,sum)
    b<-apply(cbind(env$over[,4:dim(env$over)[2]],0),1,sum)
    toKeep<-!Reduce("|",lapply(seq(length(x)),function(x) (b==1 & a==1)))
    print(x)
    print(sum(toKeep))
    envMeka<-list()
    envMeka$categories<-Filter(function(y) { ! y %in% x},env$categories)
    colnames(env$heights)<-env$categories
    envMeka$prc$normData<-qn(env$heights[toKeep,envMeka$categorie])
    envMeka$prc$eigenVectors<-prcomp(t(envMeka$prc$normData))$rotation
    envMeka
}
