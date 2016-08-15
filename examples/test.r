modulous<-function(x,m)
    {t1<-floor(x/m)
     (x-t1*m)}

not<-function(x){
    sapply(x,function(x){if (x) FALSE else TRUE })}

buildTests<-function(tests){
    function(args){with(as.list(args),all(eval(tests)))}}

combinTests<-buildTests(tests=quote(c(
                             n>=c,
                             n>0,
                             modulous(c,1)==0,
                             modulous(n,1)==0)))

permutation<-function(n,c){
    if(combinTests(c(n=n,c=c))) gamma(n+1)/gamma(c+1)
    else NA}

combination<-function(n,c){
    if(combinTests(c(n=n,c=c))) gamma(n+1)/(gamma(c+1)*gamma(n-c+1))
    else NA}

binom<-function(n,p){
    ret<-c(0,0)
    for(i in sapply(seq(n+1)-1,function(x){choose(n,x)*(1-p)^x*p^(n-x)}))
        {
            ret<-rbind(ret,c(pdf=i,cdf=ret[length(ret)][1]+i))
        }
    return(ret[2:length(ret[,1]),])
    }

pois<-function(n,p){
    ret<-c(0,0)
    for(i in sapply(seq(n+1)-1,function(x){(p*n)^x/factorial(x)*exp(-p*n)}))
        {
            ret<-rbind(ret,c(pdf=i,cdf=ret[length(ret)][1]+i))
        }
    return(ret[2:length(ret[,1]),])
    }


dataList<- rbind(c(10,12,1),
                 c(15,16,2),
                 c(11,23,1),
                 c(15,15,1))

getDistance<-function(x,p1,p2,n=2){
    sum(abs((x[p1,]-x[p2,])^n))^(1/n)}
    
flatten<-function(x){
    x[upper.tri(x,diag=TRUE)]}

getMatrixSize<-function(l){
    i=1
    dimn<-NA
    while(TRUE){    
        if(l <= sum(seq(i))){
            if(l == sum(seq(i)))dimn<- i
            break}
        i<-i+1}
    return(dimn)}

unflatten<-function(x){
    dimn<-getMatrixSize(length(x))
    A<-matrix(0,dimn,dimn)
    A[upper.tri(A,diag=TRUE)]<-x
    A<-A+t(A)
    diag(A)<-diag(A)/2
    A}

generateCovariance<-function(std,cor){
    rep<-c()
    for(i in seq(length(std)))
        for(j in seq(length(std)))
                rep<-rbind(rep,std[i]*std[j]*cor[i,j])
    return(matrix(rep,length(std),length(std) ))}

Rxyz<-unflatten(c(1,-0.5,1,0.5,0.5,1))

Cxyz<-generateCovariance(c(10,1,0.5),Rxyz)
#require("MASS")

A<-mvrnorm(10000,c(0,0,0),Cxyz)       
plot(A[,1:2])
