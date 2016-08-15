main<-function(){
    spec = matrix(c(
        'fileLocation','f', 1,"character",
        'saveFile','o', 1,"character",
        'pc1','a', 1,"integer",
        'pc2','b', 1,"integer",
        'nlc','n', 1,"integer"
    ),byrow=TRUE,ncol=4)
    args=getopt::getopt(spec=spec)
    makePCA(args$fileLocation,c(args$pc1,args$pc2),args$nlc)

    (function(x){length(which((x-mean(x))/sqrt(var(x))<(-1)))})(pcs$rotation[,2])
}


qn <-function(data){
    shape<-dim(data)
    sequence<-apply(data,2,order)
    reverseSequence<-unlist(apply(sequence,2,order))
    ranks<-apply(matrix(unlist(lapply(seq(shape[2]),function(i,x,y) x[y[,i],i],data,sequence)),ncol=shape[2]),1,sum)/shape[2]
    apply(reverseSequence,2,function(x) ranks[x])
}

standPCAPrep <-function(data,v="non"){
  switch(v,
         rowSumOne=t(apply(data,1, function(x) {x/sum(x)})),
         colSumOne=apply(data,2, function(x) x/sum(x)),
         row=t(apply(data,1, function(x) (x-mean(x))/var(x))),
         col=apply(data,2, function(x) (x-mean(x))/var(x)),
         rows1=t(apply(data,1, function(x) x/var(x))),
         cols1=apply(data,2, function(x) x/var(x)),
         colQn=qn(data),
         log10=log10(data),
         log2=log2(data),
         non=t(data))
}

makePCA<-function(dataFile,pc,n=3,title="",...){
## Load Data Into a data frame
cdata<-read.table(dataFile,header=TRUE)
## Get the length of the data
l<-length(cdata)
## n is the number of proceding columns you want for non data

stats<-cdata[seq(n)]
data<-as.matrix(apply(cdata[seq((n+1),l)],2, function(x) as.vector((unlist(x)))))
labels<-colnames(data)
normData<-qn(data)
pcs=prcomp(t(normData))

x<-t(as.matrix(pcs$rotation)) %*% as.matrix(normData)

## Ploting the Principle Components Crossed with the init data
#d1<-data.frame(x[pc[1],])
#d2<-data.frame(x[pc[2],])
#plot(t(d1),t(d2),...)
#title(main=title)
#text(x[pc[1],],x[pc[2],],labels=labels,cex=0.7,pos=3)
x
}
