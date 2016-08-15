fastMot2<-function(env,m1,m2,xlim=c(-64,64)){
    motifs2<-c(IUPACtoBase(m1),IUPACtoBase(m2))
    cList2<-sapply(motifs2,compliment)
    lM2<-lapply(motifs2,grep,env$fasta)
    lC2<-lapply(cList2,grep,env$fasta)
    plotFastMot(env,motifs2,cList2,lM2,lC2,1,2,xlim=xlim)
}


heightStat<-function(inList){
    with(inList,{
        if(!is.na(distribution[[1]])){
        distribution<-distribution[ !(distribution> (-1 *nchar(consenusIUPAC(mota)))
                     & distribution< (  nchar(consenusIUPAC(motb))))]

        height<-getHeights(distribution)
        data.frame(num=num,
                   max=max(height),
            max10=mean(sort(height,decreasing = TRUE)[2:10]),
                   mean=mean(height),
                   sd=sd(height),
                   loc=which.max(height),
                   adjloc=which.max(height)+min(distribution)
                   )
    }
        else{
            data.frame(num=num,
                   max=NA,
                   max10=NA,
                   mean=NA,
                   sd=NA,
                   loc=NA,
                   adjloc=NA
                   )
        }
    })
}

qstem<-function(b,title=NULL,xlim=c(-32,32),q=FALSE){

    b<-b[b<xlim[2]&b>xlim[1]]
    if(! q){
            par(mar=c(2,2,2,0))
        stem(min(b):max(b),getHeights(b),xlim=xlim,main=title)
    }
    b
}



qstemNorm<-function(b,norm,title=NULL,xlim=c(-32,32),q=FALSE){

    b<-b[b<xlim[2]&b>xlim[1]]
    if(! q){
            par(mar=c(2,2,2,0))
        stem(min(b):max(b),norm(getHeights(b)),xlim=xlim,main=title)
    }
    b
}

hist2Motifs<-function(env,m1,m2,reg){
    motifs2<-c(IUPACtoBase(m1),IUPACtoBase(m2))
    cList2<-sapply(motifs2,compliment)
    lM2<-lapply(motifs2,grep,env$fasta)
    lC2<-lapply(cList2,grep,env$fasta)
    motifHist(env$fasta,motifs2,cList2,lM2,lC2,1,2,env$reg[,reg])
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

makeStem<-function(y,xlim=c(-32,32),...){
    yn<-y[y>xlim[1]&y<xlim[2]]
    yt<-getHeights(yn);
    
    stem(do.call(seq,as.list(range(yn))),yt,xlim=xlim,...)
    #stem(xlim,yt,xlim=xlim,...)
}

#' motifs<-c(genEboxCombs(),mList)
#' cList<-sapply(motifs,compliment)
#' locationsM<-lapply(motifs,grep,env$fasta)
#' locationsC<-lapply(cList,grep,env$fasta)
#' 
#' i=11
#' name<-plotFastMot(env,motifs,cList,locationsM,locationsC,7,13)
#'
#' @export 
plotFastMot<-function(env20,mList,cList,locationsM,locationsC,n1,n2,reg=c("ECFC","Erythroid","HSC","Leukemia"),xlim=c(-32,32)){
par(mfrow=c(2,2),mar=c(2,2,3,2), oma=c(0, 0, 2, 0))
makeStem(motifHist(env20$fasta,mList,cList,locationsM,locationsC,n1,n2,env20$reg[,reg[1]]),main=reg[1],xlim)
makeStem(motifHist(env20$fasta,mList,cList,locationsM,locationsC,n1,n2,env20$reg[,reg[2]]),main=reg[2],xlim)
makeStem(motifHist(env20$fasta,mList,cList,locationsM,locationsC,n1,n2,env20$reg[,reg[3]]),main=reg[3],xlim)
makeStem(motifHist(env20$fasta,mList,cList,locationsM,locationsC,n1,n2,env20$reg[,reg[4]]),main=reg[4],xlim)
title(paste(consenusIUPAC(mList[n1]),"-",consenusIUPAC(mList[n2])),outer=TRUE)
paste0(consenusIUPAC(mList[n1]),"-",consenusIUPAC(mList[n2]))
}


mplotFastMot<-function(env20,mList,cList,locationsM,locationsC,n1,n2,reg){
par(mfrow=c(ceiling(length(n2)/2),2),mar=c(2,2,3,2), oma=c(0, 0, 2, 0))
for (n in n2){
    makeStem(motifHist(env20$fasta,mList,cList,locationsM,locationsC,n1,n,env20$reg[,reg]),main=paste(consenusIUPAC(mList[n1]),"-",consenusIUPAC(mList[n])))
}
title(reg,outer=TRUE)
#paste0(consenusIUPAC(mList[n1]),"-",consenusIUPAC(mList[n2]))
}


getMinDistance<-function(fasta,motifa,motifb,loc){
    a<-cbind(which(loc),selectRegionsByDistance(IUPACtoBase(motifa),IUPACtoBase(motifb),loc,env$fasta))
    b<-a[!is.na(a[,2]),]
    b    
}

getLoc2Motif<-function(minDistance){
    t<-getHeights(minDistance[,2])
    b<-which.max(t)
    tm<-t
    tm[b]<-0
    c<-which.max(tm)
    c(b+min(minDistance)-1,c+min(minDistance)-1)    
}

selMotifReg<-function(minDistance,...){
    getMin<-function(loc,minDistance)
        minDistance[which(minDistance[,2]==loc),1]
    sort(unique(unlist(lapply(list(...),getMin,minDistance))))
}

findMotifDist<-function(fasta,motifa,motifb,loc,...){
    a<-getMinDistance(fasta,motifa,motifb,loc)
    #getLoc2Motif(a))
    selMotifReg(a,...)
}
