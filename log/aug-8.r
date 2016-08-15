library(CCCA)
library(Biostrings)
source("~/r-workspace/project.r")
source("~/r-workspace/project-variables.r")
source("~/r-workspace/ccca.r")

source("http://bioconductor.org/biocLite.R")
biocLite(org.Hs.eg.db)


source("http://bioconductor.org/biocLite.R")
biocLite("GO")

library(org.Hs.eg.db)

x <- org.Hs.egCHRLOC
# Get the entrez gene identifiers that are mapped to chromosome locations
mapped_genes <- mappedkeys(x)

biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")

require(TxDb.Hsapiens.UCSC.hg19.knownGene)

exons.db <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by='gene')

egs <- unlist(  mget(symbols[ symbols %in% keys(org.Hs.egSYMBOL2EG) ],org.Hs.egSYMBOL2EG) )

str(exons.db)


exons <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by='gene')

t.db<-transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene,"gene")

txdb<-TxDb.Hsapiens.UCSC.hg19.knownGene
cds_by_tx1 <- cdsBy(txdb, "tx", use.names=TRUE)


select(org.Hs.eg.db,exons.db@partitioning@NAMES, "SYMBOL","ENTREZID")

# ensure all genes have appropriate strand values


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

t<-genGeneTSS(exons)

regions<-genomicRegions(t$chr,t$tss,t$strand,1000,5000,1000000)

a<-lapply(contexts,function(x) greatGeneAssoc(env$over[env$reg[,x],c(1,2,3)],regions,t))

names(a)<-contexts
b<-sort(unique(unlist(lapply(a,function(x) x[,"name"]))))

env<-getPRC20(2)



setdiff(b,a[[1]]$name)

om<-matrix(FALSE,length(b),length(a))
colnames(om)<-contexts
rownames(om)<-b

for(i in contexts){
    om[as.character(a[[i]]$name),i]=TRUE
}
    

regs[,as.numeric(exons.db@unlistData@strand[regs[,1]])]

## speed up genomicRegions function
    swapif<-function(x){
        if((x[1]<=x[2]))
            x
        else
            cbind(x[2],x[1])
    }

chr<-t$chr

tss<-t$tss
strand<-t$strand
levels(strand)<-c(-1,1)
strand<-as.numeric(strand)
proxUp<-1000
proxDown<-5000
distal<-1000000
basalDomains<-t(apply(
    cbind(tss-strand*proxUp,tss+strand*proxDown),1,swapif))




bound<-na.omit(do.call(rbind,lapply(levels(chr),function(lev){
    y<-basalDomains[chr==lev,]
    if(sum(chr==lev)==1){
        print(y)
        return (cbind(y[1],y[2],max(0,y[1]-distal),(y[2]+distal)))
        }
    else if(sum(chr==lev)<1){
        return (cbind(NA,NA,NA,NA))
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
        return(cbind(y[,1],y[,2],lb,ub))
        ##return(cbind(NA,NA))
    }
})))
bc<-na.omit(do.call(c,sapply(levels(chr),function(lev){
    if(sum(chr==lev)<1)
         NA
    else
        as.character(chr[chr==lev])
})))
data.frame(bc,bound)


## if true current minimum
## if false continue

## if true minimum = lower
## if false continue

## if true  then 0
## if false then distal



genomeDomains<-do.call(
        rbind,
        lapply(seq(n),extendsMax,chrom,basalDomains,tss,distal))
    return(genomeDomains)



exons.db@unlistData@strand

exons.db@unlistData@ranges[exons.db@partitioning[[1]]]

length(exons.db@unlistData@strand@values)*2

str(exons.db@unlistData@ranges)

t<-c(0,exons.db@partitioning@end)



(regs[,2]-regs[,1])[1:10]


exons.db@unlistData@strand@lengths[1:10]

inverse.rle(exons.db@unlistData@strand[regs[1:10,1]])

sapply(seq(10),function(i) exons.db[[i]]@strand@values)

x<-Rle(factor(rep(c(1,2,1,1,3,4,5),each=2),levels=0:5))

inverse.factor.rle<-function (x) 
{
    if (is.null(le <- x$lengths) || is.null(v <- x$values) || 
        length(le) != length(v)) 
        stop("invalid 'rle' structure")
    rep.int(v, le)
}

as.numeric(exons.db@unlistData@strand@values[1:10])




exons.db@unlistData@strand@values[as.numeric(exons.db@unlistData@strand)]

sapply(seq(10),function(i) as.numeric(exons.db@strand@values))
