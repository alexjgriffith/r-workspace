######## Preferred Distances
######## Libraries
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

######## Variables
## Load PCA UDM UNI
load("~/Dropbox/UTX-Alex/Paper/Raw Data/peakLocations.RData")
## Load JASPAR motifs

load("~/r-workspace/data/jaspar.RData")

######## Functions

#' IUPAC to base
#'
#' Transforms an IUPAC dna sequence into a string of neucleotides.
#' Combo characters are replaced by their compoenents in bracket.
#' For example CNC -> C[ACGT]C. With the rl flag set to TRUE the
#' expanded set of sequences will be returnd. For example CNC
#' would become c("CAC","CCC","CGC","CTC").
#' @param char the input string of IUPAC characters
#' @param rl flag to return all variants
#' @return Character containing neucleotides + brackets
#' @examples
#' IUPACtoBase("HGATAA")
#' IUPACtoBase("CANNTG")
#' IUPACtoBase("CANNTG",TRUE)
IUPACtoBase<-function (char, rl = FALSE){
    IUPAC <- strsplit(char, "")[[1]]
    IUPACCharacters <- list("A", "C", "G", "T", c("A", "G"), 
        c("C", "T"), c("C", "G"), c("A", "T"), c("G", "T"), c("A", 
            "C"), c("C", "G", "T"), c("A", "G", "T"), c("A", 
            "C", "T"), c("A", "C", "G"), c("A", "C", "G", "T"))
    names(IUPACCharacters) <- c("A", "C", "G", "T", "R", "Y", 
        "S", "W", "K", "M", "B", "D", "H", "V", "N")
    if (any(!IUPAC %in% names(IUPACCharacters))) {
        stop("Input string contains non IUPAC characters.")
    }
    vals <- sapply(IUPAC, function(x) {
        (IUPACCharacters[x])
    })
    if (rl) {
        Base <- c("")
        for (i in vals) {
            kit <- c()
            for (j in Base) kit <- c(kit, paste(j, i, sep = ""))
            Base <- kit
        }
    }
    else {
        Base <- c("")
        for (i in vals) if (length(i) == 1) 
            Base <- paste(Base, i, sep = "")
        else Base <- paste(Base, "[", do.call(paste, as.list(c(i, 
            sep = ""))), "]", sep = "")
    }
    return(Base)
}

#' Complement
#'
#' Takes a string of neucleotides (ACTG) and returns the complement
#' (TGAC). This works on strings only, for DNAStringSet please refer
#' to the `Biostrings::complement` function. To transform IUPAC characters
#' into neucleotide form refer to IUPACtoBase.
#' @param string the input string of nucleotides
#' @return The composite string
#' @examples
#' complement(IUPACtoBase("CANNTG"))
#' complement("CA[AC]TT[ACGT]GG")
complement<-function (string){
    chars <- c("A", "G", "C", "T", "[", "]")
    names(chars) <- c("T", "C", "G", "A", "]", "[")
    paste(rev(sapply(strsplit(string, "")[[1]], function(x) {
        (chars[x])
    })), collapse = "")
}

#' Get Jaspar
#'
#' To be used in conjuction with the jaspar.RData found in `extdata`.
#' Pass `getJASPAR` a composite motif from the jaspar database
#' and it will return the name of the motif.
#' @param name composite motifs IUPAC characters only
#' @return the name of the motif if it came from the jaspar database
#' @examples
#' # load the required Jaspar motif data
#' load(system.file("exdata","jaspar.RData",package="mT1"))
#' getJASPAR(jaspar$names[1])
#' @export
getJASPAR<-function(name){
    with(jaspar,{
        if(is.null(jfid)||is.null(jnames))
            stop(paste0("jaspar.RData must be loaded. ",
                        "load(system.file(\"exdata\",\"jaspar.RData\"",
                        ",package=\"mT1\"))"))
        jfid[jnames==name]
    })
}

#' Motif Distance
#'
#' Find the distances between two motifs within a subset of the genome
#' represented by fasta. 
#' 
#' @param fasta  a DNAStringSet
#' @param motif1 Character of IUPAC characters
#' @param motif2 Character of IUPAC characters
#' @return two columns , the first being the index the second being the
#' distance
#' @examples
#' ## load fasta file
#' load(system.file("exdata","fasta.RData",package="mT1"))
#' ## Find distances
#' distances<-motifDistance(fasta,"CANNTG","HGATAA")
#'
#' ## Determine the range of the motif locations and plot the
#' ## frequency
#' width<-nchar(fasta[1]) # all fasta strings should be the same width
#' r<-seq(-width+1,width-1)
#' y<-combHeights(r,distances[,2])[[1]]
#' plot(r,y,main="CANNTG-HGATAA")
#' @export
motifDistance<-function(fasta,motif1,motif2){
    motifs<-c(motif1,motif2)
    mloc<-lapply(motifs,function(x) grep(IUPACtoBase(x), fasta))
    cloc<-lapply(motifs,function(x) grep(complement(IUPACtoBase(x)),fasta))
    names(mloc)<-motifs
    names(cloc)<-motifs
    diffMotif(fasta,motifs,mloc,cloc)
}

#' Motif Difference
#' 
#' Find the distances between two motifs within a subset of the genome
#' represented by fasta. 
#' @param fasta DNAStringSet
#' @param motifs list of motif Character Atomics
#' @param mloc list of locations of each motif in fasta
#' @param cloc list of locations of each compositemotif in fasta
#' @param i index 1
#' @param 2 index 2
#' @param min numerical atomic, minimum length of mloc[[i or j]]
#' @param combine Merged|First
#' @return two columns , the first being the index the second being the
#' distance
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
    l1c<-gregexpr(complement(IUPACtoBase(n1)), fasta[sharedc,],ignore.case=TRUE)
    l2c<-gregexpr(complement(IUPACtoBase(n2)), fasta[sharedc,],ignore.case=TRUE)
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
        dchar<-nchar(n2)-nchar(n1)
        mh<-keep(l1m,l2m,function(i,j){i-j},sharedm)
        ch<-keep(l1c,l2c,function(i,j){j-i+dchar},sharedc)
        p2<-ch[!ch[,1] %in% both,]
        p3<-mh[!mh[,1] %in% both,]        
        x<-merge(as.data.frame(ch),as.data.frame(mh))
        ret<-rbind(do.call(cbind,c(x[order(x[,1]),])),p2,p3)
    }        
    colnames(ret)<-c("loc","pos")
    clean<-function(df){
        x<-df[,2]
        change<- x>= (- nchar(n2)) & x<= nchar(n1)        
        df[(!change), ]
    }
    clean(ret[order(ret[,1]),])
}

#' Motif PDF
#'
#' Generate a list of locations for a single motif.
#' @param fasta DNAStringSet from Biostrings
#' @param motif single string, containing only DNA IUPAC characters
#' @return two columns , the first being the instance the second being the
#' distance
#' @examples
#' ## load fasta file
#' load(system.file("exdata","fasta.RData",package="mT1"))
#' pdf<-motifPDF(fasta,"CANNTG")
#' x<-seq(1:100)
#' y<-combHeights(x,pdf[,2])[[1]]
#' plot(x,y,main="CANNTG",ylab="Frequency",xlab="Index")
#' @export
motifPDF<-function(fasta,motif){    
    mloc<-grep(IUPACtoBase(motif), fasta)
    cloc<-grep(complement(IUPACtoBase(motif)),fasta)
    findLocs(fasta,mloc,cloc,motif)
}

#' Find Motif Locations
#'
#' Finds the location of motifs on each indici of fasta. The results can
#' be used to build the emperical PDF of motif locations.
#' @param fasta DNAStringSet
#' @param mloc fasta indicies that have  motifs
#' @param cloc fasta indicies that have compliment motifs
#' @param n1 motif name, Atomic Character of only IUPAC characters
#' @param combine If two motifs are found on the same index Merged|First
#' @return two columns , the first being the instance the second being the
#' distance
findLocs<-function(fasta,mloc,cloc,n1,combine="Merged"){
    l1m<-gregexpr(IUPACtoBase(n1), fasta[mloc],ignore.case=TRUE)
    l1c<-gregexpr(complement(IUPACtoBase(n1)),fasta[cloc],ignore.case=TRUE)
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
        x<-merge(as.data.frame(ch),as.data.frame(mh))
        ret<-rbind(do.call(cbind,c(x[order(x[,1]),])),p2,p3)
    }
    colnames(ret)<-c("loc","pos")
    ret[order(ret[,1]),]        
}

#' Make Title
#' 
#' generic makes title for a plot
makeTitle<-function(x,...){
    UseMethod("makeTitle",x)
}

#' Make Title mT1
#'
#' Generates a title based on the composite motifs of a combination
#' @param obj mT1 object
#' @param i indici of obj$combs
#' @return atomic Character with motif info
makeTitle.mT1<-function(obj,i){
    with(obj,{
        do.call(paste,append(lapply(combs[i,],function(x) motifs[x]),list(sep="-")))
    })
}

#' main
#'
#' returns the pertanent information in an object
#' @examples
#' load(system.file("exdata","objMT1.RData",package="mT1"))
#' main(1)
#' main(objMT1)
#' @export
main<-function(x,...){
    UseMethod("main",x)
}

#' main default
#' 
#' without a class simply return x
#' @examples
#' main(1)
#' @export
main.default<-function(x){
    x
}

#' main mT1
#' 
#' returns the pertanent infromation concerning obj
#' @param obj mT1 object
#' @return data.frame(<motif1>,<motif2>,<max pvalue><loc max pvalue>)
#' @examples
#' load(system.file("exdata","objMT1.RData",package="mT1"))
#' main(objMT1)
#' @export
main.mT1<-function(obj){
    main<-data.frame(t(apply(obj$combs[!obj$sig,],1,
                             function(x) obj$motifs[x])))
    if(length(main)>0){
        colnames(main)<-c("motif1","motif2")
        mpv<-sapply(which(!obj$sig),function(x)
            min(obj$pvalue[[x]]))
        mpvl<-sapply(which(!obj$sig),function(x)
            which.min(obj$pvalue[[x]])-300)
        return(data.frame(main,mpv,mpvl))
    }
    else
        return (NULL)
}


#' print mT1
#'
#' Prints a summary of the mT1 object.
#' @param obj mT1 object
#' @export
print.mT1<-function(obj){
    #cat("mT1\n\n")
    cat("motifs:  ")
    cat(paste0(obj$motifs[1],"\n"))
    cat(paste("        ",obj$motifs[2:length(obj$motifs)],collapse="\n"))
    cat("\n\n")
    cat(paste0("combs: ",length(obj$combs), "\n"))
    cat(paste0("sufficent: ",sum(!test$sig)),"\n")
    main<-main.mT1(obj)
    print(main)
}

#' plot mT1
#'
#' Used to plot the results of mT1. Note this does not accept variable
#' arguments to `plot.default`.
#' @param obj mT1 object
#' @param motif1 motif in set analyzed
#' @param motif2 motif in set analyzed
#' @param i index value
#' @export
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
        plot(hs[[i]],type="l",ylab="Frequency",
             Main(motifs[combs[i,1]],"-",motifs[combs[i,2]]))
        plot(getHeights(hs[[i]]),type="l",xlab="Height",
             ylab="Frequency",main="Freq")
        plot(pvalue[[i]],type="l",ylab="p-value",xlab="Index",
             main="Tests")
    })
}


#' mT1
#'
#' Generates a mT1 object from DNAStringSet and a vector of motifs
#'
#' @param fasta DNAStringSet
#' @param motifs list motifs being compared, IUPAC chars only
#' @param verbose flag to print motifs as completed
#' @param cl cluster from makeForkCluster
#' @return a mT1 object
#' @examples
#' library(Biostrings) # needed to get fasta data from genomic co-ords
#' library(BSgenome.Hsapiens.UCSC.hg19) # genome
#' ## load a set of Jaspar motifs as strings
#' load(system.file("exdata","jaspar.RData",package="mT1"))
#' 
#' ## Example set of peaks
#' load(system.file("exdata","peaks.RData",package="mT1"))
#' 
#' ## Transform genomic co-ords into neucleotides
#' genome<-BSgenome.Hsapiens.UCSC.hg19
#' fasta<-getSeq(genome,peaks$chr,start=peaks$start,end=peaks$end)
#' 
#' ## Motifs to compare
#' motifs<-c("CANNTG","GATAA",jaspar$jsublM[1:8])
#' 
#' ## Find the preferred distances between `motifs` under the peaks
#' objMT1<-mT1(fasta,motifs)
#' 
#' ## plot based on string
#' plot(objMT1,"CANNTG","GATAA")
#' 
#' ## Plot based on printed order
#' print(objMT)
#' plot(objMT1,i=1)
#' @export
mT1<-function(fasta,motifs,verbose=FALSE,cl=NULL){
    ## find which fasta indicies have the motifs of interest
    mloc<-lapply(motifs,function(x) grep(IUPACtoBase(x), fasta))
    cloc<-lapply(motifs,function(x) grep(complement(IUPACtoBase(x)),fasta))
    names(mloc)<-motifs
    names(cloc)<-motifs
    ## combination of all motifs
    tofind<-t(combn(seq(length(motifs)),2))
    ## check to see if diffMotif should be preformed in parallell
    if(is.null(cl)){
        t1<-apply(tofind,1,function(x){
            if(verbose)
                print(unlist(lapply(x,function(x) motifs[x])) )
            diffMotif(fasta,motifs,mloc,cloc,x[1],x[2])
        })
    }
    else{
        t1<-parApply(cl,tofind,1,function(x){
            if(verbose)
                print(unlist(lapply(x,function(x) motifs[x])) )
            diffMotif(fasta,motifs,mloc,cloc,x[1],x[2])
        })
    }
    ## Which motif combs had co-occuances
    large<-sapply(t1,function(x) all(unlist(is.na(x))))
    ## Determine the individual PDFs for each motif
    a<-lapply(motifs,function(x) findLocs(fasta,mloc[[x]],cloc[[x]],x))
    mT1<-list(diff=t1,dens=a,sig=large,motifs=motifs,combs=tofind)
    ## determine the p-values for each motif comb
    prob<-apply(cbind(tofind,seq(dim(tofind)[1])),1,
                function(x) ePD(t1[[x[3]]],a[[x[1]]],a[[x[2]]],
                                nchar(motifs[[x[1]]]),
                                Nchar[[x[2]]])))
    ## Build and return the mT1 object
    mT1<-append(mT1,list(hs=lapply(prob,function(x) x[["hs"]]),
                mp=lapply(prob,function(x) x[["mp"]]),
                pvalue=lapply(prob,function(x) x[["pvalue"]])))
    attr(mT1,"class")<-"mT1"
    mT1

}

#' binomial test
#'
#' A wrapper around binom.test that ensures that n!=0 and that k is at
#' least a minimum number.
#' @param k count
#' @param n total
#' @param p prob
#' @param min minimum k cut off
#' @return p-value
#' @examples
#' @export
btest<-function(k,n,p,min=10){
    if(n>0) # n must be greater than 0
        return(0)
    else if(k>min) # ensure that there are more k than min
        return(log(binom.test(k,n,p)$p.value,10))
    else
        return(0)
}    

#' Expected Motif Probobility
#'
#' Determine the expected probobility of finding two motifs at a specific
#' distance based on the emparical pdfs of the individual motifs
#' @param a numerical vector PWM of motif a
#' @param b numerical vector PWM of motif b
#' @param width width of fasta indicies
#' @return expectation for each motif distance
#' @examples
#' 
#' @export
eMP<-function(a,b,width){
    y<-combHeights(seq(1,width),a,b)
    mod<-as.integer(do.call(convolve,append(y,list(type="open"))))
    mod/sum(mod)
}

#' ePD
#'
#' Determine the pvalue from a vetor of distances and PWMs from the
#' two motifs being compared
#' @param t1 width 2 numerical vector of distances
#' @param a numerical vector PWM of motif a
#' @param b numerical vector PWM of motif b
#' @param width width of fasta indicies
#' @return
ePD<-function(t1,a,b,width){
    block<-function(x,min,max){
        c(x[1:(min-1)],rep(0,max-min+1),x[(max+1):length(x)])
    }
    if(any(is.na(t1))){        
        return(list(hs=NA,mp=NA,pvalue=NA))
    }
    hs<-combHeights(seq(-with+1,width-1),t1[,2])[[1]]
    n<-length(t1[,2])
    if(max(hs)>n){
        stop("n >hs")
    }
    mp<-eMP(a,b,width)
    pvalue<-mapply(function(a,b)btest(a,n,b),hs,mp)
    list(hs=hs,mp=mp,pvalue=pvalue)
}

#' Combination Heights
#'
#' Determines the number of occurances in each member ... along x.
#' @param x seqence of possible values represented as indicies
#' @param ... vectors of values to be mapped
#' @examples
#' x<-seq(20)
#' y<-runif(200,1,20)
#' par(mfrow=c(1,2))
#' plot(y)
#' plot(x,combHeights(x,y)[[1]])
#' @export
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

########### Application
## plot(test,"CANNTG","GATAA")

## plot(test,i=10)

## plot(combHeights(seq(-299,299),test$diff[[12]])[[1]],type="l")

## plot(test$hs[[12]],type="l")

## cl<-makeForkCluster(4)

## stopCluster(cl)

#motifs<-unique(c("CANNTG","GATAA",jaspar$jsublM[1:3]))

#system.time(test<-mT1(PCA$fasta,motifs))

