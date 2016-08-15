library("functional")

readCategories<-Compose(read.table,unlist,as.character)

getSwapCats<-function(cats,names){    
    swCat<-cats
    names(swCat)<-names
    function(x){as.character(unlist(lapply(x,function(x) swCat[x])))}
}

stringToSwap<-function(x)do.call(getSwapCats,splitZip(createZip(strsplit(x," ")[[1]])))

simpleSwapFun<-function(char){
    stringToSwap(paste (rev(strsplit(char, " ")[[1]]),collapse=" "))
}

orderBed<-function(ret)
    ret[order(as.character(ret[,1]),ret[,3]),]

qstem<-function(y,xlim=c(-32,32),...){
    yn<-y[y>xlim[1]&y<xlim[2]]
    yt<-getHeights(yn);    
    stem(do.call(seq,as.list(range(yn))),yt,xlim=xlim,...)
}

# utils
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



lzip<-function(...){
    apply(mapply(function(...)list(...),...),2,as.list)
}

pass<-function(x) x

createZip<-function(x)
    Map(function(i,j)cbind(x[i],x[j]),seq(1,length(x),2),seq(2,length(x),2))

splitZip<-function(inList)
    list(sapply(inList,"[",1),sapply(inList,"[",2))


addColnames<-function(matrix,colnames){
    colnames(matrix)<-colnames
    matrix
}

addRownames<-function(matrix,rownames){
    rownames(matrix)<-rownames
    matrix
}

addNames<-function(matrix,colnames,rownames=colnames,list=NULL){
    if(is.null(list)){
        colnames(matrix)<-colnames
        rownames(matrix)<-rownames
    }
    else{
        names(matrix)<-colnames
    }
    matrix        
}

normalize<-function(x)(x-mean(x))/(sqrt(var(x)))

makeLogic<-function(loc,size){
    x=rep(FALSE,size)
    x[loc]<-TRUE
    x
}


qn <-function(data){
    shape<-dim(data)
    sequence<-apply(data,2,order)
    reverseSequence<-unlist(apply(sequence,2,order))
    ranks<-apply(matrix(unlist(lapply(seq(shape[2]),function(i,x,y) x[y[,i],i],data,sequence)),ncol=shape[2]),1,sum)/shape[2]
    apply(reverseSequence,2,function(x) ranks[x])}


pca<-function(data,norm="qn"){
    if( is.function(norm))
       normData<-norm(data)
    else
        normData<-switch(norm,
               rowSumOne=t(apply(data,1, function(x) {x/sum(x)})),
               colSumOne=apply(data,2, function(x) x/sum(x)),
               normRow=t(apply(data,1, function(x) (x-mean(x))/var(x))),
               normCol=apply(data,2, function(x) (x-mean(x))/var(x)),
               rowsVarOne=t(apply(data,1, function(x) x/var(x))),
               colsVarOne=apply(data,2, function(x) x/var(x)),
               qn=qn(data),
               none=data,
               data)
    prc<-prcomp(t(normData))$rotation
    list(normData=normData,eigenVectors=prc)
}


