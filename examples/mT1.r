######## Preferred Distances
######## Libraries
library(CCCA)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

######## Sources
source(system.file("scripts","eboxFrequency.r",package="CCCA"))
source("~/r-workspace/project/ccca.r")
source("~/r-workspace/project/project.r")
source("~/r-workspace/project/project-variables.r")

######## Variables
## Load PCA UDM UNI
load("~/Dropbox/UTX-Alex/Paper/Raw Data/peakLocations.RData")
## Load JASPAR motifs

jpwmsp<-loadPWM("~/masters/r-workspace/CCCA/inst/exdata/jaspar.pwm","jaspar")[,3]

jpwms<-lapply(jpwmsp,function(x)matrix(c(x),nrow=4,byrow=TRUE))
jnames<-sapply(jpwms,CCCA::PWMtoCons)
jnames<-gsub("^N*||N*$","",jnames)


jl<-sapply(jnames,function(x)length(grepMotifs(x,PCA$fasta)))

jsublM<-jnames[order(jl,decreasing=TRUE)]

jid<-loadPWM("~/masters/r-workspace/CCCA/inst/exdata/jaspar.pwm","jaspar")[,2]
jfid<-unlist(lapply(strsplit(unlist(jid),"\t"),function(x) x[2]))

getJASPAR<-function(name){
    jfid[jnames==name]
}

motifs<-c("GATAA","CANNTG",jsublM[1:10])

######## Functions




diffMotif<-function(fasta,motifs,mloc,cloc,i=1,j=2,min=200,combine="Merged"){
    n1<-motifs[[i]]
    n2<-motifs[[j]]
    sharedm<-intersect(mloc[[n1]],mloc[[n2]])
    sharedc<-intersect(cloc[[n1]],cloc[[n2]])
    both<-intersect(sharedc,sharedm)
    ## if there are not enough overlaping locations return NA
    if(length(union(sharedm,sharedc))<min)
        return(NA);
    l1m<-gregexpr(IUPACtoBase(n1), fasta[sharedm,],ignore.case=TRUE)
    l2m<-gregexpr(IUPACtoBase(n2), fasta[sharedm,],ignore.case=TRUE)
    l1c<-gregexpr(compliment(IUPACtoBase(n1)), fasta[sharedc,],ignore.case=TRUE)
    l2c<-gregexpr(compliment(IUPACtoBase(n2)), fasta[sharedc,],ignore.case=TRUE)
    ## Locations under the same peak are combined in two ways
    ## First: Keep the smallest distance
    ## Merged: Keep all distances   
    if(combine=="First"){
        mh<-cbind(loc=sharedm,pos=mapply(function(a,b) {x<-unique(c(outer(a,b,Vectorize(function(i,j) i-j))));x[which.min(abs(x))]},l1m,l2m))
        ch<-cbind(loc=sharedc,pos=mapply(function(a,b) {x<-unique(c(outer(a,b,Vectorize(function(i,j) j-i))));x[which.min(abs(x))]},l1c,l2c))
        p2<-ch[!ch[,1] %in% both,]
        p3<-mh[!mh[,1] %in% both,]        
        ret<-rbind(p2,p3)
        if(length(both)>0){
            p1<-cbind(both,apply(cbind(ch[ch[,1] %in% both,2],mh[mh[,1] %in% both,2]),1,function(x) x[which.min(abs(x))]))
            ret<-rbind(ret,p1)
        }
    }
    else if(combine=="Merged"){
        keep<-function(l1,l2,fun,shared){
            pos<-mapply(function(a,b) {unique(c(outer(a,b,Vectorize(fun))))},l1,l2,SIMPLIFY=FALSE)
            lens<-lapply(pos,length)
            loc<-unlist(mapply(function(x,n){rep(x,n)} ,shared,lens,SIMPLIFY=FALSE))
            cbind(loc=loc,pos=unlist(pos))
        }
        mh<-keep(l1m,l2m,function(i,j){i-j},sharedm)
        ch<-keep(l1c,l2c,function(i,j){j-i},sharedc)
        p2<-ch[!ch[,1] %in% both,]
        p3<-mh[!mh[,1] %in% both,]        
        ret<-rbind(p2,p3)

        if(length(both)>0){
            p1<-rbind(ch[!ch[,1] %in% both,],mh[mh[,1] %in% both,])
            ret<-rbind(ret,p1)
        }
    }        
    colnames(ret)<-c("loc","pos")
    clean<-function(df){
        x<-df[,2]
        change<- x>= (- nchar(n2)) & x<= nchar(n1)        
        #bounds<- x > -100 & x < 100
        df[(!change), ]#& bounds,]
    }
    clean(ret[order(ret[,1]),])
}



findLocs<-function(fasta,mloc,cloc,n1,combine="Merged"){
    l1m<-gregexpr(IUPACtoBase(n1), fasta[mloc],ignore.case=TRUE)
    l1c<-gregexpr(compliment(IUPACtoBase(n1)),fasta[cloc],ignore.case=TRUE)
    both<-intersect(mloc,cloc)
    if(combine=="First"){
        mh<-cbind(loc=mloc,sapply(l1m,function(x) x[which.min(abs(x))]))
        ch<-cbind(loc=cloc,sapply(l1c,function(x) x[which.min(abs(x))]))    
        p2<-ch[!ch[,1] %in% both,]
        p3<-mh[!mh[,1] %in% both,]
        ret<-rbind(p2,p3)
        if(length(both)>0){        
            p1<-cbind(both,apply(cbind(ch[ch[,1] %in% both,2],mh[mh[,1] %in% both,2]),1,function(x) x[which.min(abs(x))]))
            ret<-p1
        }
    }
    else if(combine=="Merged"){
        keep<-function(l,loc){
            pos=unlist(lapply(l,unique))
            len<-sapply(l,length)
            eloc<-unlist(mapply(function(x,n)rep(x,n),loc,len,SIMPLIFY=FALSE))
            cbind(loc=eloc,pos=pos)
        }
        
        mh<-keep(l1m,mloc)
        ch<-keep(l1c,cloc)
        p2<-ch[!ch[,1] %in% both,]
        p3<-mh[!mh[,1] %in% both,]
        ## ret<-rbind(p2,p3)
        ## ch<<-ch
        ## mh<<-mh
        ## if(length(both)>0){
        ##     p1<-merge.data.frame(as.data.frame(ch[!ch[,1] %in% both,]),as.data.frame(mh[mh[,1] %in% both,]))
        ##     print(p1)
        ##     ##p1<-do.call(rbind,lapply(do.call(seq,as.list(range(p1[,1]))),function(x){ y<-x[p1[,2]==x,]; y[!duplicated(y,FALSE),]} ))
            
        ##     ret<-rbind(ret,p1)
        ## }
        x<-merge(as.data.frame(ch),as.data.frame(mh))
        ret<-rbind(do.call(cbind,c(x[order(x[,1]),])),p2,p3)
        print(head(ret))
    }
    colnames(ret)<-c("loc","pos")
    ret[order(ret[,1]),]        
}

main.mT1<-function(obj){
    main<-data.frame(t(apply(obj$combs[!obj$sig,],1,function(x) obj$motifs[x])))
    if(length(main)>0){
        colnames(main)<-c("motif1","motif2")
        mpv<-sapply(which(!obj$sig),function(x) min(obj$pvalue[[x]]))
        mpvl<-sapply(which(!obj$sig),function(x) which.min(obj$pvalue[[x]])-300)
        return(data.frame(main,mpv,mpvl))
    }
    else
        return (NULL)
    }

print.mT1<-function(obj){
    #cat("mT1\n\n")
    cat("motifs:  ")
    cat(paste0(obj$motifs[1],"\n"))
    cat(paste("        ",obj$motifs[2:length(obj$motifs)],collapse="\n"))
    cat("\n\n")
    cat(paste0("combs: ",length(obj$combs), "\n"))
    cat(paste0("sufficent: ",sum(!test$sig)),"\n")
    main<-main(obj)
    print(main)
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
    plot(density(dens[[n1]][,2]),main=motif1)
    plot(density(dens[[n2]][,2]),main=motif2)
    plot(mp[[i]],type="l",ylab="p",main="Convolution")
    plot(hs[[i]],type="l",ylab="Frequency",main= paste0(motifs[combs[i,1]],"-",motifs[combs[i,2]]))
        plot(getHeights(hs[[i]]),type="l",xlab="Height",ylab="Frequency",main="Freq")
    plot(pvalue[[i]],type="l",ylab="p-value",xlab="Index",main="Tests")
    })
}



mT1<-function(fasta,motifs,cl=NULL){
    mloc<-lapply(motifs,function(x) grep(IUPACtoBase(x), fasta))
    cloc<-lapply(motifs,function(x) grep(compliment(IUPACtoBase(x)),fasta))
    names(mloc)<-motifs
    names(cloc)<-motifs
    tofind<-t(combn(seq(length(motifs)),2))
    if(is.null(cl)){
        t1<-apply(tofind,1,function(x){
            print(unlist(lapply(x,function(x) motifs[x])) )
            diffMotif(fasta,motifs,mloc,cloc,x[1],x[2])
        })
    }
    else{
        t1<-parApply(cl,tofind,1,function(x){
            print(unlist(lapply(x,function(x) motifs[x])) )
            diffMotif(fasta,motifs,mloc,cloc,x[1],x[2])
        })
    }
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
        #print(c(k,n,p))
        if(n==0)
            return(0)
        else if(k>10)
            return(log(binom.test(k,n,p)$p.value,10))
           #return(log(1-pbinom(k,n,p),10))
        else
            return(0)
    }    
    combHeights<-function(x,...){
        min<-min(x)
        max<-max(x)
        lapply(list(...),function(h){        
            y<- rep(0,length(x))
            for(i in h){
                if(i>=min & i<=max)
                    y[i-min +1]<-y[i-min +1]+1
            }
            y
        })
    }
    if(any(is.na(t1))){        
        return(list(hs=NA,mp=NA,pvalue=NA))
    }
    ##hs<-getHeights(t1[t1>-100&t1<99],c(-99,99))
    hs<-combHeights(seq(-299,299),t1[,2])[[1]]
    n<-length(t1[,2])
    if(max(hs)>n){
        stop("n >hs")
    }
    #ah<-cbind(do.call(seq,as.list(range(a[,2]))),getHeights(a[,2]))
    #bh<-cbind(do.call(seq,as.list(range(b[,2]))),getHeights(b[,2]))
    mod<-as.integer(do.call(convolve,append(combHeights(seq(1,300),a,b),list(type="open"))))
    mp<-mod/sum(mod)
    ##mod<-as.integer(convolve(ah[ah[,1]<201,2],bh[bh[,1]<201,2],
    ##                         type="open"))[]
    ##mp<-mod[101:299]/sum(mod)
    pvalue<-block(mapply(function(a,b)test(b,n,a),mp,hs),100-n2,100+n1)
    list(hs=hs,mp=mp,pvalue=pvalue)
}



########### Application
plot(test,"CANNTG","GATAA")

plot(test,i=12)

plot(combHeights(seq(-299,299),test$diff[[12]])[[1]],type="l")

plot(test$hs[[12]],type="l")

cl<-makeForkCluster(4)

stopCluster(cl)

motifs<-unique(c("CANNTG","GATAA",jsublM[1:300]))

system.time(test<-mT1(PCA$fasta,motifs[1:10]))

## adjust the flip
x<-main.mT1(testR)

x[order(x[,3])[1:20],]

CANNTG           ASMGGAA
ASMGGAA            CACCTG

"YCTTATCTS"
"RTCTGGHW"

getJASPAR("TAATCC")

plot(testR,i=67)
