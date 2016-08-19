# Functions to be added to CCCA
library(CCCA)
library(functional)
library(ggplot2)
library(parallel)
library(grid)
#library(BSgenome.Hsapiens.UCSC.hg19)



#' Replace EZ ID
#' @description takes a table returned by david and swaps the EZ IDs
#' with official gene symbols
#' @param table David GO table
#' @param pvalue minimum cut off for keeping GO terms
#' @returns a table of GO terms with IDs replaced with names
#' @examples
#' genesEZ<-list(list("GATA","TAL1","RUNX1"))
#' contexts<-c("test")
#' user<-testuser@example.com
#' davidResults<-genGO(user,contexts,genesEZ)
#' ## Replace the ids with names and only keep the 
#' ## GO terms that have a p-value less than 1e-3
#' replaceEZID(davidResults[[1]][[1]],pvalue=1e-3)
#' @export
replaceEZID<-function(table,pvalue=1e-1){
	## Takes a string of EZ ids and splits it. Then the ids are 
	## swapped with gene names and returns a string of gene names
	ezidGene<-function(ezid){
		ezidSplit<-strsplit(ezid,", ")[[1]]
		geneName<-select(org.Hs.eg.db,ezidSplit, 
			"SYMBOL" ,"ENTREZID")[,2]
		paste(geneName,collapse=", ")
	}

    sigReg<-table[,5]<pvalue
    geneNames<-sapply(table[sigReg,6],ezidGene)
    stable<-table[sigReg,]
    stable[,6]<-geneNames
    stable
}

genGeneTSS<-function(exons=NULL){
    if(!(require(org.Hs.eg.db) & require(TxDb.Hsapiens.UCSC.hg19.knownGene)))
        NULL
    if(is.null(exons))
    exons<-exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by='gene')
    t<-c(0,exons@partitioning@end)
    regs<-cbind(t[1:(length(t)-1)]+1,t[2:length(t)])
    df<-data.frame(chr=as.character(exons@unlistData@seqnames[regs[,1]]),
                   tss=0,
                   strands=c("+","-")[as.numeric(exons@unlistData@strand[regs[,1]])],
                   id=exons@partitioning@NAMES,
                       name=select(org.Hs.eg.db,exons@partitioning@NAMES, "SYMBOL","ENTREZID")$SYMBOL)
 
    hold<-rep(0,dim(regs)[[1]])
    hold[df$strands=="-"]<-regs[df$strands=="-",2]
    df$tss[df$strands=="-"]=exons@unlistData@ranges@start[hold[df$strands=="-"]]+exons@unlistData@ranges@width[hold[df$strands=="-"]]
    hold[df$strands=="+"]<-regs[df$strands=="+",1]
    df$tss[df$strands=="+"]=exons@unlistData@ranges@start[hold[df$strands=="+"]]
    df[!is.na(as.character(df$name)),]
}


selectGenes<-function(matrix,context,df=NULL){
    if(is.null(df)){
        df<-data.frame(rownames(matrix)[matrix[,context]==1])
        colnames(df)<-context
        return(df)
    }
    else
        return(rownames(matrix)[matrix[,context]==1])
}

## significantly faster than genomicRegions in CCCA
## also i think it is a bit clearer than the previous version
genomicRegions<-function(chr,tss,strand,proxUp,proxDown,distal){
    swapif<-function(x){
        if((x[1]<=x[2]))
            x
        else
            cbind(x[2],x[1])
    }
    levels(strand)<-c(-1,1)
    strand<-as.numeric(strand)
    basalDomains<-t(apply(
    cbind(tss-strand*proxUp,tss+strand*proxDown),1,swapif))
    bound<-na.omit(do.call(rbind,lapply(levels(chr),function(lev){
        y<-basalDomains[chr==lev,]
        if(sum(chr==lev)==1){
            return (cbind(max(0,y[1]-distal),(y[2]+distal)))
        }
        else if(sum(chr==lev)<1){
            return (cbind(NA,NA))
        }
        else {
            y<-y[order(y[,1]),]
            len<-dim(y)[1]
            lower<-c(0,y[1:(len-1),2])
            upper<-c(y[2:len,1],y[len,1]+distal)
            extBD<-cbind(y,lower,upper,y[,1]-distal,y[,2]+distal)
            lbl<-rowSums(cbind(extBD[,1]>extBD[,3],
                               extBD[,1] > extBD[,3]& extBD[,3] < extBD[,5],
                               extBD[,1] > extBD[,3]& extBD[,3] < extBD[,5]&extBD[,5]<0))
            lb<-rep(0,len)
            lb[lbl==0]=extBD[lbl==0,1]
            lb[lbl==1]=extBD[lbl==1,3]
            lb[lbl==2]=extBD[lbl==2,5]
            lb[lbl==3]=0
            ubl<-rowSums(cbind(extBD[,2]<extBD[,4],
                               extBD[,2]<extBD[,4]&extBD[,4]>extBD[,6]))
            ub<-rep(0,len)
            ub[ubl==0]=extBD[ubl==0,2]
            ub[ubl==1]=extBD[ubl==1,4]
            ub[ubl==2]=extBD[ubl==2,6]
            return(cbind(lb,ub))
            ##return(cbind(NA,NA))
        }
})))
    bc<-na.omit(do.call(c,sapply(levels(chr),function(lev){
        if(sum(chr==lev)<1)
            NA
        else
            as.character(chr[chr==lev])
    })))
    geneRangeToRegion(bound,bc)
}

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

geneMatrix<-function(over,reg,geneRegions,geneList,id="name"){
    contexts<-colnames(reg)
    a<-lapply(contexts,function(x) greatGeneAssoc(env$over[reg[,x],c(1,2,3)],regions,geneList))
    names(a)<-contexts
    b<-sort(unique(unlist(lapply(a,function(x) x[,id]))))

    om<-matrix(FALSE,length(b),length(a))
    colnames(om)<-contexts
    rownames(om)<-b

    for(i in contexts){
        om[as.character(a[[i]][,id]),i]=TRUE
    }
    om
}


geneMatrix_old<-function(over,reg,regions,geneList,id="name2"){
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
