### nov-30
rename3D<-function(matrix,colnames,rownames){
    for(i in seq(dim(matrix)[3])){
        matrix[,,i]<-addNames(matrix,colnames,rownames)        
    }
    matrix
}

quot<-function(x)x[1]/x[2]
norm<-function(x)apply(x,1,function(x-mean(x)/sqrt(var(x))))
norm<-function(x)  apply(x,1,log)

o<-addNames(((1-correl[,,30])+(1-correl[,,1]))/70,categories)
on<-o#norm(o)
no<-hclust(as.dist(on),"single")$order
ono<-on[no,no]
par(mar=c(1,1,7,7))
image(x=seq(1:22),y=seq(1:22),t(ono),col=colorRampPalette(c("red","white","blue"))(1024) ,ylab="",xlab="",yaxis=NULL,axes=FALSE,zlim = range(sort(on)[((length(o)*0)+1):(length(o)*1)]))
axis(3,at=seq(22),labels=FALSE)
axis(4,at=seq(22),labels=FALSE)
text(seq(22)-0.33,par("usr")[4]+.5,pos=4,labels=swapFun(categories[no]),srt=75,xpd=TRUE,cex=1.25)
text(par("usr")[2]+.1,seq(22),pos=4,labels=swapFun(categories[no]),xpd=TRUE,cex=1.25)


i=2
j=10
l<-lm(y~x,data.frame(y=correl[i,j,],x=pvalues))$coef[2]
plot(pvalues,correl[i,j,])
lines(pvalues,pvalues* l[2] + l[1])

y<-seq(getRange(ol,.99)[1],getRange(ol,.99)[2],length=12)
image(y,0,matrix(seq(12)),ylab="",xlab="")

hist(ol)

## dec 1
fastaLocations<-paste(">",bed20$chro,":",bed20$start,"-",bed20$end,sep="")
