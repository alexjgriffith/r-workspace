### Do the prefered distance stuff for Eboxs
library(CCCA)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
#library(org.Hs.eg.db)
#library(RDAVIDWebService)
#library(xlsx)
#library(org.Hs.eg.db)
source(system.file("scripts","eboxFrequency.r",package="CCCA"))
source("~/r-workspace/project/ccca.r")
source("~/r-workspace/project/project.r")
source("~/r-workspace/project/project-variables.r")

load("~/Dropbox/UTX-Alex/Paper/Raw Data/peakLocations.RData")

motifs<-c(genEboxCombs(),"CANNTG","GATAA")

locs<-lapply(motifs,function(motif) CCCA::grepMotifs(motif,PCA$fasta))

mloc<-lapply(motifs,function(x) grep(IUPACtoBase(x), PCA$fasta))

cloc<-lapply(motifs,function(x) grep(compliment(IUPACtoBase(x)), PCA$fasta))

names(mloc)<-motifs
names(cloc)<-motifs

easyDiffMotif<-function(fasta,motif1,motif2){
    motifs<-c(motif1,motif2)
    mloc<-lapply(motifs,function(x) grep(IUPACtoBase(x), fasta))
    cloc<-lapply(motifs,function(x) grep(compliment(IUPACtoBase(x)),fasta))
    names(mloc)<-motifs
    names(cloc)<-motifs
    diffMotif(fasta,motifs,mloc,cloc)
}

diffMotif<-function(fasta,motifs,mloc,cloc,i=1,j=2,min=200){
    n1<-motifs[[i]]
    n2<-motifs[[j]]

    sharedm<-intersect(mloc[[n1]],mloc[[n2]])
    sharedc<-intersect(cloc[[n1]],cloc[[n2]])
    if(length(union(sharedm,sharedc))<min)
        return(NA);
    l1m<-gregexpr(IUPACtoBase(n1), fasta[sharedm,],ignore.case=TRUE)
    l2m<-gregexpr(IUPACtoBase(n2), fasta[sharedm,],ignore.case=TRUE)
    l1c<-gregexpr(compliment(IUPACtoBase(n1)), fasta[sharedc,],ignore.case=TRUE)
    l2c<-gregexpr(compliment(IUPACtoBase(n2)), fasta[sharedc,],ignore.case=TRUE)
    mh<-cbind(loc=sharedm,pos=mapply(function(a,b) {x<-unique(c(outer(a,b,Vectorize(function(i,j) i-j))));x[which.min(abs(x))]},l1m,l2m))
    ch<-cbind(loc=sharedc,pos=mapply(function(a,b) {x<-unique(c(outer(a,b,Vectorize(function(i,j) j-i))));x[which.min(abs(x))]},l1c,l2c))
    both<-intersect(sharedc,sharedm)
    p2<-ch[!ch[,1] %in% both,]
    p3<-mh[!mh[,1] %in% both,]
    ret<-rbind(p2,p3)
    if(length(both)>0){
        p1<-cbind(both,apply(cbind(ch[ch[,1] %in% both,2],mh[mh[,1] %in% both,2]),1,function(x) x[which.min(abs(x))]))
        ret<-rbind(ret,p1)
    }

    colnames(ret)<-c("loc","pos")
    clean<-function(df){
        x<-df[,2]
        change<- x>= (- nchar(n1)) & x<= nchar(n2)        
        #bounds<- x > -100 & x < 100
        df[(!change), ]#& bounds,]
    }
    clean(ret[order(ret[,1]),])
}


h<-diffMotif(PCA$fasta,motifs,mloc,cloc,11,12)

plotH<-function(h,motif1,motif2,x=NULL,...){
    if(is.null(x))
        x<-do.call(seq,as.list(range(h)))
    y<-sapply(x,function(x) sum(h==x))
    y[x> (- nchar(motif1)) & x< nchar(motif2)]<-0
    stem(x,y,main=paste(motif1,motif2),...)
}


moreMotifs<-gsub(">","",unlist(CCCA::loadPWM("~/Dropbox/UTX-Alex/Paper/Analysis/homer_Leukemia_PCA_6.txt")[,"name"]))

motifs<-c("CANNTG","GATTA",moreMotifs)


jpwms<-loadPWM("~/masters/r-workspace/CCCA/inst/exdata/jaspar.pwm","jaspar")[,3]

jid<-loadPWM("~/masters/r-workspace/CCCA/inst/exdata/jaspar.pwm","jaspar")[,2]

jfid<-unlist(lapply(strsplit(unlist(jid),"\t"),function(x) x[2]))



jnames<-sapply(jpwms,CCCA::PWMtoCons)
jnames<-gsub("^N*||N*$","",jnames)

TBP<-jnames[which(jfid=="TBP")]

## filter motifs to make sure thay map to at least 50 locs


jl<-sapply(jnames,function(x)length(grepMotifs(x,PCA$fasta[PCA$reg[,"Leukemia"]])))

jsublM<-jnames[jl>100]

motifs<-c("GATAA","CANNTG",jsublM[1:20])

mloc<-lapply(motifs,function(x) grep(IUPACtoBase(x), PCA$fasta[PCA$reg[,"Leukemia"]]))
cloc<-lapply(motifs,function(x) grep(compliment(IUPACtoBase(x)), PCA$fasta[PCA$reg[,"Leukemia"]]))
names(mloc)<-motifs
names(cloc)<-motifs

tofind<-t(combn(seq(length(motifs)),2))


Hall<-mT1(PCA$fasta,motifs)

Hall<-apply(tofind,1,function(x) diffMotif(PCA$fasta[PCA$reg[,"Leukemia"]],motifs,mloc,cloc,x[1],x[2]))


whichHall<-sapply(Hall,function(x)length(x)>200)


df<-data.frame(Hall[[which(whichHall)[4]]])

plotH(df,motifs[tofind[which(whichHall)[4],1]],motifs[tofind[which(whichHall)[4],2]])

plot(density(df[,2],1))

ggplot(df,aes(x=pos))+geom_density()


findLocs<-function(fasta,mloc,cloc,n1){
    l1m<-gregexpr(IUPACtoBase(n1), fasta[mloc],ignore.case=TRUE)
    l1c<-gregexpr(compliment(IUPACtoBase(n1)),fasta[cloc],ignore.case=TRUE)
    
    mh<-cbind(loc=mloc,sapply(l1m,function(x) x[which.min(abs(x))]))
    ch<-cbind(loc=cloc,sapply(l1c,function(x) x[which.min(abs(x))]))    
    both<-intersect(mloc,cloc)
    p2<-ch[!ch[,1] %in% both,]
    p3<-mh[!mh[,1] %in% both,]
    ret<-rbind(p2,p3)
    if(length(both)>0){        
        p1<-cbind(both,apply(cbind(ch[ch[,1] %in% both,2],mh[mh[,1] %in% both,2]),1,function(x) x[which.min(abs(x))]))
        ret<-rbind(ret,p1)
    }
    colnames(ret)<-c("loc","pos")
    ret[order(ret[,1]),]        
}




quickPlotPD<-function(fasta,motif1,motif2){                             
    motifs<-c(motif1,motif2)
    mloc<-lapply(motifs,function(x) grep(IUPACtoBase(x), fasta))
    cloc<-lapply(motifs,function(x) grep(compliment(IUPACtoBase(x)),fasta))
    names(mloc)<-motifs
    names(cloc)<-motifs
    t1<-diffMotif(fasta,motifs,mloc,cloc,)[,2]
    hs<-getHeights(t1[t1>-100&t1<99],c(-99,99))    
    a<-findLocs(fasta,mloc[[1]],cloc[[1]],motifs[1])
    b<-findLocs(fasta,mloc[[2]],cloc[[2]],motifs[2])
    ah<-cbind(do.call(seq,as.list(range(a[,2]))),getHeights(a[,2]))
    bh<-cbind(do.call(seq,as.list(range(b[,2]))),getHeights(b[,2]))
    mod<-as.integer(convolve(ah[ah[,1]<201,2],bh[bh[,1]<201,2],
                             type="open"))[]
    mp<-mod[101:299]/sum(mod)
    par(mfrow=c(3,2))
    plot(density(b[,2]))
    plot(density(a[,2]))
    plot(mp,type="l")
    plot(hs,type="l")
    plot(getHeights(hs),type="l")
    plot(mapply(function(a,b)dbinom(b,sum(hs),a,log=TRUE),mp,hs),type="l")
}


print.mT1<-function(obj){
    #cat("mT1\n\n")
    cat("motifs:  ")
    cat(paste0(obj$motifs[1],"\n"))
    cat(paste("        ",obj$motifs[2:length(obj$motifs)],collapse="\n"))
    cat("\n\n")
    cat(paste0("combs: ",length(obj$combs), "\n"))
    cat(paste0("sufficent: ",sum(sapply(obj$diff,function(x) all(is.na(x)))), "\n"))
    main<-data.frame(t(apply(obj$combs[!obj$sig,],1,function(x) obj$motifs[x])))    
    colnames(main)<-c("motif1","motif2")
    mpv<-sapply(which(!obj$sig),function(x) min(obj$pvalue[[x]]))
    mpvl<-sapply(which(!obj$sig),function(x) which.min(obj$pvalue[[x]])-100)
    print(data.frame(main,mpv,mpvl))
    
}

plot.mT1<-function(obj,motif1=NULL,motif2=NULL,i=NULL){
    if(!is.null(i)){
        main<-data.frame(t(apply(obj$combs[!obj$sig,],1,function(x) obj$motifs[x])))
        motif1<-main[i,1]
        motif2<-main[i,2]
    }
    with(obj,{
        n1<-which(motifs==motif1)
        n2<-which(motifs==motif2)
        i<-union(which(combs[,1]==n1 & combs[,2] ==n2),
        which(combs[,1]==n2 & combs[,2] ==n1))[1]
    par(mfrow=c(3,2))
    plot(density(dens[[n1]][,2]))
    plot(density(dens[[n2]][,2]))
    plot(mp[[i]],type="l")
    plot(hs[[i]],type="l")
        plot(getHeights(hs[[i]]),type="l")
        print(log(sum(hs[[i]]),10))
    plot(pvalue[[i]],type="l")
    })
}


mT1<-function(fasta,motifs){
    mloc<-lapply(motifs,function(x) grep(IUPACtoBase(x), fasta))
    cloc<-lapply(motifs,function(x) grep(compliment(IUPACtoBase(x)),fasta))
    names(mloc)<-motifs
    names(cloc)<-motifs
    tofind<-t(combn(seq(length(motifs)),2))
    t1<-apply(tofind,1,function(x) diffMotif(fasta,motifs,mloc,cloc,x[1],x[2]))
    large<-sapply(t1,function(x) all(unlist(is.na(x))))
    a<-lapply(motifs,function(x) findLocs(fasta,mloc[[x]],cloc[[x]],x))
    mT1<-list(diff=t1,dens=a,sig=large,motifs=motifs,combs=tofind)
    prob<-apply(cbind(tofind,seq(dim(tofind)[1])),1,function(x) ePD(t1[[x[3]]],a[[x[1]]],a[[x[2]]],nchar(motifs[[x[1]]]),nchar(motifs[[x[2]]])))
    mT1<-append(mT1,list(hs=lapply(prob,function(x) x[["hs"]]),
                mp=lapply(prob,function(x) x[["mp"]]),
                pvalue=lapply(prob,function(x) x[["pvalue"]])))
    attr(mT1,"class")<-"mT1"
    mT1
}

ePD<-function(t1,a,b,n1=0,n2=0){
    block<-function(x,min,max){
        c(x[1:(min-1)],rep(0,max-min+1),x[(max+1):length(x)])
    }
    test<-function(k,n,p){
        if(k>4)
            return(dbinom(k,n,p,log=TRUE))
        else
            return(1)
        }
    if(any(is.na(t1))){        
        return(list(hs=NA,mp=NA,pvalue=NA))
    }
    hs<-getHeights(t1[t1>-100&t1<99],c(-99,99))
    n<-length(t1)
    ah<-cbind(do.call(seq,as.list(range(a[,2]))),getHeights(a[,2]))
    bh<-cbind(do.call(seq,as.list(range(b[,2]))),getHeights(b[,2]))
    mod<-as.integer(convolve(ah[ah[,1]<201,2],bh[bh[,1]<201,2],
                             type="open"))[]
    mp<-mod[101:299]/sum(mod)
    pvalue<-block(mapply(function(a,b)test(b,n,a),mp,hs),100-n2,100+n1)
    list(hs=hs,mp=mp,pvalue=pvalue)
}

motifs<-c("GATAA","CANNTG",jsublM[1:30])
system.time(test<-mT1(PCA$fasta,motifs))

plot(test,"CANNTG","GATAA")

quickPlotPD(PCA$fasta,"GATAA","CANNTG")

quickPlotPD(BED$fasta,"AGGCGG","SCACTG")

quickPlotPD(BED$fasta,"CANNTG","AGGCCG")

quickPlotPD(UNI$fasta,TBP,"CACATG")

dbinom()

mod[104]

plot(mod/sum(mod))

1,

n<-(c(seq(1:100)^2, rev(seq(1:99)^2)))


as.integer(convolve(getHeights(a[,2])[1:100],getHeights(b[,2])[1:100],type="open"))


getHeights<-function (h, range = c(min(h), max(h))) 
{
    rep <- rep(0, (range[2] - range[1] + 1))
    for (i in h){
        rep[i - range[1] + 1] <- rep[i - range[1] + 
                                     1] + 1
        }
    rep
}



plot(getHeights(a[,2]))


### Ebox prefferences
#### moved to examples/eboxPref.r

