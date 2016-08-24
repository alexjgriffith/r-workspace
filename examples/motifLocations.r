## to be used with mT1 currently broken

motifLocations<-function(fasta,motifs){
    forwards<-function(motif,sub)
        gregexpr(IUPACtoBase(motifs),
                 fasta,ignore.case=TRUE)
    backwards<-function(motif,sub)
        gregexpr(compliment(IUPACtoBase(motifs)),
                 fasta,ignore.case=TRUE)
    mloc<-lapply(motifs,function(x) grep(IUPACtoBase(x), fasta))
    cloc<-lapply(motifs,function(x) grep(compliment(IUPACtoBase(x)),fasta))  
    forw<-mapply(forwards,motifs,mloc,SIMPLIFY=FALSE)
    back<-mapply(backwards,motifs,cloc,SIMPLIFY=FALSE)
    names(mloc)<-motifs
    names(cloc)<-motifs
    ret<-list(mloc=mloc,cloc=cloc,forw=forw,back=back,motifs=motifs)
    class(ret)<-"motifLocations"
    ret
}


mT1.motifLocations<-function(locs){
    with(locs,{
    tofind<-t(combn(seq(length(motifs)),2))
    t1<-apply(tofind,1,function(x){
        print(unlist(lapply(x,function(x) motifs[x])) )
        diffMotif.motifLocations(locs,x[1],x[2])
        })
    large<-sapply(t1,function(x) all(unlist(is.na(x))))
    a<-lapply(motifs,function(x) findLocs.motifLocations(locs,x))
    mT1<-list(diff=t1,dens=a,sig=large,motifs=motifs,combs=tofind)

    prob<-apply(cbind(tofind,seq(dim(tofind)[1])),1,function(x) ePD(t1[[x[3]]],a[[x[1]]],a[[x[2]]],nchar(motifs[[x[1]]]),nchar(motifs[[x[2]]])))
    mT1<-append(mT1,list(hs=lapply(prob,function(x) x[["hs"]]),
                mp=lapply(prob,function(x) x[["mp"]]),
                pvalue=lapply(prob,function(x) x[["pvalue"]])))
    attr(mT1,"class")<-"mT1"
    mT1
    })
}

findLocs.motifLocations<-function(obj,n1,combine="Merged"){
    with(obj,{
        print(n1)
        l1m<-forw[[n1]][mloc[[n1]]]    
        l1c<-back[[n1]][cloc[[n1]]]
    both<-intersect(mloc[[n1]],cloc[[n1]])
    if(combine=="First"){
        mh<-cbind(loc=mloc[[n1]],sapply(l1m,function(x) x[which.min(abs(x))]))
        ch<-cbind(loc=cloc[[n1]],sapply(l1c,function(x) x[which.min(abs(x))]))    
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
            pos=unlist(l)
            len<-sapply(l,length)
            eloc<-unlist(mapply(function(x,n)rep(x,n),loc,len,SIMPLIFY=FALSE))
            cbind(loc=eloc,pos=pos)
        }
        
        mh<-keep(l1m,mloc[[n1]])
        ch<-keep(l1c,cloc[[n1]])
        p2<-ch[!ch[,1] %in% both,]
        p3<-mh[!mh[,1] %in% both,]
        ret<-rbind(p2,p3)      
        if(length(both)>0){
            p1<-rbind(ch[!ch[,1] %in% both,],mh[mh[,1] %in% both,])
            ret<-rbind(ret,p1)
        }
    }
    colnames(ret)<-c("loc","pos")
    ret[order(ret[,1]),]        
})}

diffMotif.motifLocations<-function(obj,i=1,j=2,min=200,combine="Merged"){
    with(obj,{
    n1<-motifs[[i]]
    n2<-motifs[[j]]
    sharedm<-intersect(mloc[[n1]],mloc[[n2]])
    sharedc<-intersect(cloc[[n1]],cloc[[n2]])
    both<-intersect(sharedc,sharedm)
    ## if there are not enough overlaping locations return NA
    if(length(union(sharedm,sharedc))<min)
        return(NA);
    l1m<-forw[[n1]][sharedm]    
    l2m<-forw[[n2]][sharedm]
    l1c<-back[[n1]][sharedc]
    l2c<-back[[n1]][sharedc]
    ## Locations under the same peak are combined in two ways
    ## First: Keep the smallest distance
    ## Merged: Keep all distances
    #combine="First"
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
        change<- x>= (- nchar(n1)) & x<= nchar(n2)        
        df[(!change), ]
    }
    clean(ret[order(ret[,1]),])
})}
