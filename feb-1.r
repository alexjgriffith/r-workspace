source("~/r-workspace/nov-functions.r")
source("~/r-workspace/jan-functions.r")
source("~/r-workspace/jan-variables.r")
source("~/Dropbox/thesis/R/contrib_new.R")

bed<-orderBed(bed20)

heights<-addColnames(read.table("~/r-workspace/test.txt"),categories)
prc=pca(heights)

stackedContrib(prc$eigenVectors[,4],"contrib2",mergeFun(over20,swapFunD),swapFun=swapFunB,colourOveride =swapFunC,n=6)



plotPCMat2D(pca2Matr(prc),c(1,4),categories ,swapFun, swapFunD ,swapFunC)

stemRegions=normalize(prc$eigenVectors[,4])>3

erytRegions=normalize(prc$eigenVectors[,1])>2


prc$normData[stemRegions,]

test<-addColnames(prc$normData[stemRegions,],swapFunD(categories))

corder<-addColnames(do.call(cbind,lapply(seq(22),function(i) order(heights[,i]))),categories)

test<-addColnames(apply(heights,2,rord)[stemRegions,],swapFunD(categories))

test<-addColnames(apply(heights,2,rord)[erytRegions,],swapFunD(categories))

mcol<-addColnames(do.call(cbind,lapply(unique(swapFunD(categories)), function(x) apply(cbind(test[,x],0), 1, max))),unique(swapFunD(categories)))

length(which(mcol[,"Erythroid"]>0.5))

heights[,"cd34"]

rord<-function(x) {
    y<-sort(unique(unlist(x)))
    z<-sapply(y,function(w) length(which(x>=w)))/length(x)
    names(z)<-y
    z[as.character(x)]
}




x<-matrix(floor(runif(20,1,10)),ncol=2)

qn <-function(data){
    shape<-dim(data)
    sequence<-apply(data,2,order)
    reverseSequence<-unlist(apply(sequence,2,order))
    ranks<-apply(matrix(unlist(lapply(seq(shape[2]),function(i,x,y) x[y[,i],i],data,sequence)),ncol=shape[2]),1,sum)/shape[2]
    apply(reverseSequence,2,function(x) ranks[x])}

z<-x
shape<-dim(z)
sequence<-apply(z,2,order)
reverseSequence<-unlist(apply(sequence,2,order))
ranks<-apply(matrix(unlist(lapply(seq(shape[2]),
                                  function(i,x,y)
                                      x[y[,i],i],z,sequence)),ncol=shape[2]),1,sum)/shape[2]
apply(reverseSequence,2,function(x) ranks[x])

