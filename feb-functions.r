source("~/r-workspace/nov-functions.r")
#source("~/Masters/CCCA/inst/scipts/eboxFrequency.r")


getPRC<-function(pvalue=5,control="single",treatment="treatment",forward="22"){
    categories<-readCategories(paste0("~/Dropbox/Data/categories/",forward,"-categories.txt"))
    over<-orderBed(readAFS(pvalue,control,forward=forward))
    bed<-over[,1:3]
    heights<-readUDM(pvalue,treatment,control,forward=forward)
    prc<-pca(heights)    
    list(over=over,bed=bed,heights=heights,prc=prc,categories=categories)
}

getPRC20<-function(n=1,inner=0.25){
    out<-getPRC(20,"combined")
    add<-with(out,{
    ECFC<-normalize(prc$eigenVectors[,4])<(-(n))&normalize(prc$eigenVectors[,2])<(-(n))
    ECFCA<-normalize(prc$eigenVectors[,2])<(-n)
    ECFCB<-normalize(prc$eigenVectors[,4])<(-n)
    HSC<-normalize(prc$eigenVectors[,4])>(n)
    Leukemia<-normalize(prc$eigenVectors[,1])<(-n)
    Erythroid<-normalize(prc$eigenVectors[,1])>(n)
    ALL<-rep(TRUE,length(Erythroid))
    NONE<-normalize(prc$eigenVectors[,4])>(-inner)&
        normalize(prc$eigenVectors[,4])<(inner) &
            normalize(prc$eigenVectors[,1])>(-inner) &
                normalize(prc$eigenVectors[,1])<(inner) &
                    normalize(prc$eigenVectors[,2])>(-inner)     
    ECFC.alt<-ECFC & !(Erythroid | Leukemia | HSC)
    HSC.alt<-HSC & !(Erythroid | Leukemia | ECFC)
    Leukemia.alt<-Leukemia & !(HSC | Erythroid | ECFC)
    Erythroid.alt<-Erythroid & !(HSC | Leukemia | ECFC)
    reg<-addColnames(cbind(Erythroid.alt,Leukemia.alt,HSC.alt,ECFC.alt,NONE,ALL),c("Erythroid","Leukemia","HSC","ECFC","NONE","ALL"))
    list(reg=reg,ALL=ALL,NONE=NONE,ECFC=ECFC,HSC=HSC,Leukemia=Leukemia,Erythroid=Erythroid,
         ECFC.alt=ECFC.alt,HSC.alt=HSC.alt,Leukemia.alt=Leukemia.alt,Erythroid.alt=Erythroid.alt,
         ECFCB=ECFCB,ECFCA=ECFCA)
})
    append(out,add)
}

getPRC20F<-function(out,n=1,inner=0.25){
    add<-with(out,{
    #ECFC<-function(prc=prc,n=n,inner=inner)
    #    {normalize(prc$eigenVectors[,4])<(-(n^2/2))&
    #         normalize(prc$eigenVectors[,2])<(-(n/2))}
    ECFC<-function(prc=prc,n=n,inner=inner)
        {normalize(prc$eigenVectors[,4])<(-n)&
             normalize(prc$eigenVectors[,2])<(-n)}
    
    HSC<-function(prc=prc, n=n,inner=inner){
        normalize(prc$eigenVectors[,4])>(n)
    }
    Leukemia<-function(prc=prc2, n=n,inner=inner){
        normalize(prc$eigenVectors[,1])<(-n)
    }
    Erythroid<-function(prc=prc, n=n,inner=inner){
        normalize(prc$eigenVectors[,1])>(n)
    }
    ALL<-function(prc=prc, n=n,inner=inner){
        rep(TRUE,length(Erythroid))
    }
    NONE<-function(prc, n=n,inner=inner){
        normalize(prc$eigenVectors[,4])>(-inner)&
        normalize(prc$eigenVectors[,4])<(inner) &
            normalize(prc$eigenVectors[,1])>(-inner) &
                normalize(prc$eigenVectors[,1])<(inner) &
                    normalize(prc$eigenVectors[,2])>(-inner)
    }
    ECFC.alt<-function(prc=prc, n=n,inner=inner){
        ECFC(prc,n,inner) &
            !(Erythroid(prc,n,inner) | Leukemia(prc,n,inner) | HSC(prc,n,inner))
    }
    HSC.alt<-function(prc=prc, n=n,inner=inner){
        HSC(prc,n,inner) &
            !(Erythroid(prc,n,inner) | Leukemia(prc,n,inner) | ECFC(prc,n,inner))
    }
    Leukemia.alt<-function(prc=prc, n=n,inner=inner){
        Leukemia(prc,n,inner) &
            !(HSC(prc,n,inner) | Erythroid(prc,n,inner) | ECFC(prc,n,inner))
    }
    Erythroid.alt<-function(prc=prc, n=n,inner=inner){
        Erythroid(prc,n,inner) &
            !(HSC(prc,n,inner) | Leukemia(prc,n,inner) | ECFC(prc,n,inner))
    }
    reg<-addColnames(cbind(Erythroid.alt(prc,n,inner),
                           Leukemia.alt(prc,n,inner),
                           HSC.alt(prc,n,inner),
                           ECFC.alt(prc,n,inner),
                           NONE(prc,n,inner),
                           ALL(prc,n,inner)),c("Erythroid","Leukemia","HSC","ECFC","NONE","ALL"))
    list(reg=reg)
})
    append(out,add)
}


getPRC5<-function(n=1,inner=0.25){
    out<-getPRC(5,"combined")
    add<-with(out,{
    ECFC<-normalize(prc$eigenVectors[,3])>n & normalize(prc$eigenVectors[,5])<(-n)
    #ECFC<-normalize(prc$eigenVectors[,4])<(-n)
    HSC<-normalize(prc$eigenVectors[,3])>(n) & normalize(prc$eigenVectors[,5])>(n)
    Leukemia<-normalize(prc$eigenVectors[,1])<(-n)
    Erythroid<-normalize(prc$eigenVectors[,1])>(n)
    ALL<-rep(TRUE,length(Erythroid))
    NONE<-normalize(prc$eigenVectors[,3])>(-inner)&
        normalize(prc$eigenVectors[,3])<(inner) &
            normalize(prc$eigenVectors[,1])>(-inner) &
                normalize(prc$eigenVectors[,1])<(inner) &
                    normalize(prc$eigenVectors[,5])>(-inner)&
                        normalize(prc$eigenVectors[,5])<(inner)     
    ECFC.alt<-ECFC & !(Erythroid | Leukemia | HSC)
    HSC.alt<-HSC & !(Erythroid | Leukemia | ECFC)
    Leukemia.alt<-Leukemia & !(HSC | Erythroid | ECFC)
    Erythroid.alt<-Erythroid & !(HSC | Leukemia | ECFC)
    reg<-addColnames(cbind(Erythroid.alt,Leukemia.alt,HSC.alt,ECFC.alt,NONE,ALL),c("Erythroid","Leukemia","HSC","ECFC","NONE","ALL"))
    list(reg=reg,ALL=ALL,NONE=NONE,ECFC=ECFC,HSC=HSC,Leukemia=Leukemia,Erythroid=Erythroid,
         ECFC.alt=ECFC.alt,HSC.alt=HSC.alt,Leukemia.alt=Leukemia.alt,Erythroid.alt=Erythroid.alt)
})
    append(out,add)
}

addFasta<-function(env){
    env$fasta<-getSeq(BSgenome.Hsapiens.UCSC.hg19,env$bed$chr,start=env$bed$start+150,width=300)
    env
}

getEboxInfo<-function(env){
    add<-with(env,{        
        eboxLoc<-sapply(genEboxCombs(),grepMotifs,fasta)
        peakLoc<-lapply(seq(4,dim(over)[2]),Compose(sel2(over),as.logical,which))
        names(peakLoc)<-categories
        uniqueLocs<-apply(over[,4:dim(over)[2]],1,sum)==1
        regLoc<-lapply(seq(dim(reg)[2]),Compose(sel2(reg),which))
        return(list(eboxLoc=eboxLoc,peakLoc=peakLoc,uniqueLocs=uniqueLocs,
                    regLoc=regLoc))
    })
    append(env,add)
}

getMotifInfo<-function(env,motifs){
    add<-with(env,{        
        eboxLoc<-sapply(motifs,grepMotifs,fasta)
        peakLoc<-lapply(seq(4,dim(over)[2]),Compose(sel2(over),as.logical,which))
        names(peakLoc)<-categories
        uniqueLocs<-apply(over[,4:dim(over)[2]],1,sum)==1
        regLoc<-lapply(seq(dim(reg)[2]),Compose(sel2(reg),which))
        return(list(eboxLoc=eboxLoc,peakLoc=peakLoc,uniqueLocs=uniqueLocs,
                    regLoc=regLoc))
    })
    append(env,add)
}


getConEboxInfo<-function(env){
    append(env,with(env,{
        conOver<-do.call(cbind,lapply(lapply(unique(swapFunD(categories)),function(x) do.call(cbind,lapply(which(x==swapFunD(categories)), function(x) over[,4:dim(over)[2]][,x]))),function(x) apply(cbind(x,0),1,function(x) if(sum(x)>0)1 else 0 )))
        colnames(conOver)<-unique(swapFunD(categories))
        conUniqueLocs<-apply(over[,4:dim(over)[2]],1,sum)==1
        conPeakLoc<-lapply(unique(swapFunD(categories)),function(x) do.call(unionN,lapply(which(x==swapFunD(categories)), function(x) peakLoc[[x]])))
        names(conPeakLoc)<-colnames(conOver)
        list(conOver=conOver,conUniqueLocs=conUniqueLocs,conPeakLoc=conPeakLoc)
    }))}


motifFrequency<-function(regI,mlocI,sel,subs=NULL,motifs=NULL,cols=NULL){
    getNC<-function(x)
        if(is.null(colnames(x))){names(x)}else{colnames(x) }
    applySel<-function(x,swapFunD){
        sapply(unique(swapFunD(colnames(x))),function(y) {apply(cbind(x[,swapFunD(colnames(x))==y],0),1, sum)})
    }

    if(is.null(motifs))
        motifs=names(mlocI)
    if(is.null(cols))
        cols=getNC(regI)
    motifAllFrequencyAlt<-motifFreqGen(function(){function(x,y){
        length(intersect(x,y))
    }})
    if (! is.list(regI))
        temp<-lapply(seq(dim(regI)[2]),function(x) which(regI[,x]))
    else
        temp<-regI
    if(is.null(subs))
        subs=union(unlist(temp),unlist(mlocI))
    reg<-lapply(temp,intersect,subs)
    mloc<-lapply(mlocI,intersect,subs)
    all<-sapply(mloc,length)

    count<-motifAllFrequencyAlt(motifs,cols,mloc,reg)
    countM<-rbind(sapply(reg,length),0)
    colnames(countM)<-cols;
    list(matrix=applySel(count,sel),
         colSize=applySel(countM,sel),
         rowSbize=all)   
}

findEboxByD<-function(env,sel=swapFunD){
    with(env,{
        eboxD<-motifDistances(eboxLoc,fasta)
        list(a=motifFrequency(peakLoc,eboxD,sel,which(uniqueLocs)),
             b=motifFrequency(reg,eboxD,pass))
    })}

findEboxByA<-function(env,sel=swapFunD){
    with(env,{
        list(a=motifFrequency(peakLoc,eboxLoc,sel,which(uniqueLocs)),
             b=motifFrequency(reg,eboxLoc,pass))
})}

findEboxBoth<-function(env){
    BA<-findEboxByA(env)
    BD<-findEboxByD(env)
    list(bedD=BD$a,bedA=BA$a,pcaD=BD$b,pcaA=BA$b)
    
}

swapFunM<-simpleSwapFun("tall_p1 P2B tall_p2_1 P1 tall_p2_2 P1_R tall_p3_1 P2A tall_p3_2 P2A_R jurk_sandar_1 Jurk1 jurk_sandar Jurk1_R jurk_1 Jurk2 jurk Jurk2 jurk_2 Jurk2_R rpmi_1 RPMI rpmi_2 RPMI_R cem_1 CEM cem_2 CEM_R ecfc-tsa ECFC ecfc.tsa ECFC ecfc ECFC meka MEKA cd133 CD133 cd34 CD34_1 cd34_new CD34_2 eryt_1 ProEB_A_1 eryt_2 ProEB_A_1_R eryt ProEB_A_1 eryt_f ProEB_F eryt_a ProEB_A_2 k562_1 K562 k562_2 K562_R")


readAFS<-function(pvalue=5,control="single",directory="~/Dropbox/Data/AFS/",forward=22){
    fname<-paste0(directory,forward,"_pvalue_",pvalue,"_control_",control,".txt")
    read.table(fname,header=T)
}

readUDM<-function(pvalue=5,tr="treatment",control="single",directory="~/Dropbox/Data/UDM/",forward=22){
    fname<-paste0(directory,forward,"_",tr,"_pvalue_",pvalue,"_control_",control,".txt")
    read.table(fname,header=T)
}


rawDataFiles<-function(categories)paste("/mnt/brand01-00/mbrand_analysis/data_sets/",categories,"/",categories,"_sorted.bed",sep="")


subsetPRC<-function(prc,pc,n=1,fun=">")
    do.call(fun,list(normalize(prc$eigenVectors[,pc]),n))



plotFastMot<-function(env20,mList,cList,locationsM,locationsC,n1,n2,reg=c("ECFC","Erythroid","HSC","Leukemia")){
par(mfrow=c(2,2),mar=c(2,2,3,2), oma=c(0, 0, 2, 0))
makeStem(motifHist(env20$fasta,mList,cList,locationsM,locationsC,n1,n2,env20$reg[,reg[1]]),main=reg[1])
makeStem(motifHist(env20$fasta,mList,cList,locationsM,locationsC,n1,n2,env20$reg[,reg[2]]),main=reg[2])
makeStem(motifHist(env20$fasta,mList,cList,locationsM,locationsC,n1,n2,env20$reg[,reg[3]]),main=reg[3])
makeStem(motifHist(env20$fasta,mList,cList,locationsM,locationsC,n1,n2,env20$reg[,reg[4]]),main=reg[4])
title(paste(consenusIUPAC(mList[n1]),"-",consenusIUPAC(mList[n2])),outer=TRUE)
paste0(consenusIUPAC(mList[n1]),"-",consenusIUPAC(mList[n2]))
}



makeStem<-function(y,xlim=c(-32,32),...){
    yn<-y[y>xlim[1]&y<xlim[2]]
    yt<-getHeights(yn);
    
    stem(do.call(seq,as.list(range(yn))),yt,xlim=xlim,...)
    #stem(xlim,yt,xlim=xlim,...)
}

#' select regions by distance
#' altered version,
selectRegionsByDistance<-function(motif1,motif2,locations,fasta){
    m1<-gregexpr(motif1,fasta[locations],ignore.case=TRUE)
    m2<-gregexpr(motif2,fasta[locations],ignore.case=TRUE)
    c1<-gregexpr(compliment(motif1),fasta[locations],ignore.case=TRUE)
    c2<-gregexpr(compliment(motif2),fasta[locations],ignore.case=TRUE)
    findDmin<-function(a,b){
        c<-NA
        if(a[1]==-1 & b[1]!=-1){
            c<-NA#b[which.min(abs(b))]
        }
        else if (a[1]!=-1 & b[1]==-1){
            c<-NA#a[which.min(abs(a))]
        }
        else if (a[1]!=-1 & b[1]!=-1){
            c<-outer(a,b,"-")
            c<-c[which.min(abs(c))]
        }
        else{
            c<-NA
        }
        c}
    #cd<-(mapply(findDmin,c1,c2)+1)-1
    cd<-(mapply(findDmin,c1,c2)+1)*-1
    md<-mapply(findDmin,m1,m2)
    #mapply(function(a,b) if(abs(a)<abs(b) & b!=0){ a } else {b},md,cd)
    a<-apply(cbind(md,cd),1,function(x) if(length(which(is.na(x)))>1){NA}else{min(x,na.rm=TRUE)})
    #stem(min(a):max(a),getHeights(a),xlim=c(-32,32))
    a
}


stem <- function(x,y,pch=16,linecol=1,clinecol=1,...){
if (missing(y)){
    y = x
    x = 1:length(x) }
    plot(x,y,pch=pch,...)
    for (i in 1:length(x)){
       lines(c(x[i],x[i]), c(0,y[i]),col=linecol)
    }
    lines(c(x[1]-2,x[length(x)]+2), c(0,0),col=clinecol)
}

tts<-function(x){
    paste0(x[,1],":",x[,2]-2000,"-",x[,3]+2000)
}

onlyOne<-function(geneNames){
    geneOrder<-order(geneNames)
    rorder<-order(geneOrder)
    as<-geneNames[geneOrder[1:length(geneOrder)-1]]
    bs<-geneNames[geneOrder[2:length(geneOrder)]]
    mult<-c(TRUE,as!=bs)
    mult[rorder]
}

makeGeneList<-function(geneFile="/home/agriffith/Dropbox/UTX-Alex/jan/hg19.RefSeqGenes.csv",rna=NULL){
    geneListP<-read.delim(geneFile)
    if(is.null(rna))
        rna<-data.frame(gene=geneListP$name2,fc=0)
    genesInList<-sapply(geneListP$name2,
                        function(x)
                            x %in% rna$gene)&onlyOne(geneListP$name2)
    geneListP[genesInList,]
}

makeGeneRegions<-function(
    geneFile="/home/agriffith/Dropbox/UTX-Alex/jan/hg19.RefSeqGenes.csv",
    rna=NULL){
    geneList<-makeGeneList(geneFile,rna)
    chrom<-as.character(geneList$chrom)
    tss<-as.numeric(geneList$txStart)
    strand<-geneList$strand
    levels(strand)<-c(-1,1)
    strand<-as.numeric(strand)
    genomicRegions(chrom,tss,strand,1000,5000,1000000)
}

unionRegs<-function(env,swapFun,reg){
    as.logical(mergeFun(env$over[,4:dim(env$over)[2]],swapFun)[,reg])
}

myrmna<-function(x){Filter(function(x)!is.na(x), x)}
