lzip<-function(...){
    apply(mapply(function(...)list(...),...),2,as.list)
}

#' a<-list(list(c("a","b","c")),list(c("b","c","d")),list(c("c","d","e")))
#' b<-list(list(c("a","b","c")),list(c("b","c","d")))
minter<-function(list){
    rest<-function(list,n){
        lapply(n:length(list),function(x) list[[x]])
    }
    inter<-function(list){
        intersect(unlist(list[[1]]),unlist(list[[2]]))
    }
    if(length(list)<2)
        return(list)
    else if(length(list)==2)
        return(inter(list))
    else
         minter(append(list(inter(list)),rest(list,3)))
}

sigG<-function(mix,n=1000){
    s<-seq(min(mix$x),max(mix$x),length=n)
    cbind(s,apply(cbind(mapply(function(mu,sigma,lambda)dnorm(s,mu,sigma)*lambda,mix$mu,mix$sigma,mix$lambda),0),1,sum ))
}

top<-function(test1,a1,n,k)
    names(unlist((test1[test1>a1$mu[n]+a1$sigma[n]*k])))


fn1<-paste0("cufflinks2/",c(
    "utx-Dox.1/genes.fpkm_tracking",
    "utx-Dox.2/genes.fpkm_tracking",
    "utx-ctl.1/genes.fpkm_tracking",
    "utx-ctl.2/genes.fpkm_tracking",
    "tal1-Dox.1/genes.fpkm_tracking",
    "tal1-Dox.2/genes.fpkm_tracking",
    "tal1-ctl.1/genes.fpkm_tracking",
    "tal1-ctl.2/genes.fpkm_tracking"))

J41<-lapply(fn1,function(fn) read.table(fn,header=T))


shared<-minter(lapply(J41,function(x) x$gene_id))

J41T<-lapply(J41,function(x) addNames(x$FPKM,x$gene_id,list=TRUE)[shared])

test<-do.call(cbind,lapply(J41T,function(x) x/sum(x)))

testLog<-apply(test,2,function(test) log(test[test!=0]))

mix<-lapply(testLog,function(test1) normalmixEM(test1,k=2))

sharedG<-minter(lapply(lzip(testLog,mix),function(testa)top(testa[[1]],testa[[2]],which.max(testa[[2]]$mu),-5) ))

fm<-do.call(cbind,lapply(testLog,"[",sharedG))


       
       


#cdn<-read.table("cuffnorm2/genes.count_table",header=T)

#temp<-addRownames(cdn[apply(cdn[,6:13],1,sum)!=0,6:13],cdn[apply(cdn[,6:13],1,sum)!=0,1])
#test<-temp/apply(temp,2,sum)

library(mixtools)

a<-normalmixEM(log(test[test[,3]!=0,3]),k=3)

all<-apply(test,1,function(x) sum(x==0))==0




               
test2<-log(test[test[,2]!=0,2])
test3<-log(test[test[,3]!=0,3])
test4<-log(test[test[,4]!=0,4])


a2<-normalmixEM(test2,k=3)
a3<-normalmixEM(test3,k=3)
a4<-normalmixEM(test4,k=3)

plot(a,which=3)


     
plot(sigG(a1),type="l",col="red")
lines(sigG(a2))

lines(density(test1))

n<-which.max(a1$mu)


plot(density(test1[test1>a1$mu[3]-a1$sigma[3]]))


                     
a<-intersect(
    intersect(top(test1,a1,3),top(test2,a2,3)),
    intersect(top(test3,a3,3),top(test4,a4,3) ))

top<-names(unlist(((test1-a1$mu[n])/a1$sigma[n])>0))

plot(density(log(test[a,2])))

plot(log(test[a,3]),log(test[a,4]))

plot(hclust(as.dist(1-cor(test[,]))),hang=-1)

plot(density(log(test[,4])))

plot(density(log(test[all,3])),col="red")
lines(density(log(test[all,1])),col="blue")
lines(density(log(test[all,2])),col="blue")

cor(test)

u=1/length(test[,1])

a<-sapply(seq(1,100)/10, function(E1)
    sum((u+(test[,2]-u)/var(test[,2]) - u+test[,1]*E1)^2))

b<-sapply(seq(1,100)/10, function(E1)
    sum((mean(test[,3])-test[,3]*E1 - mean(test[,1])-test[,1]*E1)^2))

c<-sapply(seq(1,100)/10, function(E1)
    sum((mean(test[,3])-test[,3]*E1 - mean(test[,4])-test[,4]*E1)^2))

plot(a,type="l")
lines(c,col="red")
lines(b,col="blue")

mydensity<-function(x,n=10000,c=(n-1)){
    mx<-max(x)
    mn<-min(x)
    step<-(mx-mn)/n
    steps<-seq(mn,mx,by=step)[1:(c+1)]
    out<-unlist(mapply(function(u,d)length(x[x>=u & x<d]),steps[1:(c)],steps[2:(c+1)]))[2:c]
    out/sum(out)
}


dbetaWrap<-function(x,n=10000,c=n,a,b){
    mx<-max(x)
    mn<-min(x)
    step<-(mx-mn)/n
    steps<-seq(mn,mx,by=step)[1:c]
    out<-dbeta(steps,a,b)[2:c]
    out/sum(out)    
}



source("~/r-workspace/nov-functions.r")


mean(test1)

var(test1)

E=0.01
mean(test1)+test1*E

bvar<-function(v,a,m){
    b<-function(a,m){
        (a-m*a)/m
    }
    v-(a*b(a,m)/((a+b(a,m))^2*(a+b(a,m)+1)))
}



a<-dbetaWrap(test2,n=1000000,c=200,0.5,1)
c<-dbetaWrap(test2,n=1000000,c=200,0.01,300)
b<-mydensity(test2,n=1000000,c=200)

plot(b,type="l")
lines(a,col="red")
lines(c,col="blue")




library(MASS)



mynorm<-rnorm(n=100,mean=10,sd=10)

nn<-(mynorm-min(mynorm)+1)/(max(mynorm)-min(mynorm))*0.9

fitdistr(test1,"beta",list(shape1=mean(test1),shape2=400))

fitdistr(test1[test1!=0],dnorm,list(mean=0,sd=1))

plot(dbeta(seq(0,1,length=100),0.0001,400))

a<-dbeta(seq(0,1,length=100),0.0001,400)[2:100]



fitdistr(a,dbeta,list(shape1=mean(test1),shape2=400))

set.seed(3)
x <- rgamma(1e5, 2, .2)
plot(density(x))
# normalize the gamma so it's between 0 & 1
# .0001 added because having exactly 1 causes fail
xt <- x / ( max( x ) + .0001 )
# fit a beta distribution to xt
library( MASS )
library( ggplot2)
fit.beta <- fitdistr( xt, "beta", start = list( shape1=2, shape2=5 ) )
x.beta <- rbeta(1e5,fit.beta$estimate[[1]],fit.beta$estimate[[2]])
## plot the pdfs on top of each other
plot(density(xt))
lines(density(x.beta), col="red" )
 
## plot the qqplots
qqplot(xt, x.beta)
