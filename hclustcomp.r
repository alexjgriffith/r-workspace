o<-lapply( seq(10000),function(i) dist(matrix(runif(10000),ncol=100)))

sl<-do.call(rbind,lapply(seq(10)*10,function(k){
    j<-lapply(o,function(i) as.dist(as.matrix(i)[1:k,1:k]))
    start=proc.time()
    for(i in j){
        fastcluster::hclust(i)
    }
    step1=proc.time()
    for(i in j){
        stats::hclust(i)
    }
    step2=proc.time()
    print(rbind(step1-start,step2-step1))
    cbind((step1-start)[[1]],(step2-step1)[[1]])
}))

#png("~/Dropbox/clustcomp.png")
plot(seq(10)*10,sl[,2],type="l",ylim=range(unlist(sl)),xlab="Matrix Size", ylab="Time (s) for 10,000 Matrix",col="red")
lines(seq(10)*10,sl[,1],type="l",col="blue")
#dev.off()
