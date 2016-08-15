library(CCCA)
library(Biostrings)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#library(GO.db)
source("~/r-workspace/project.r")
source("~/r-workspace/project-variables.r")
source("~/r-workspace/ccca.r")

geneList<-genGeneTSS()

regions<-genomicRegions(geneList$chr,geneList$tss,geneList$strand,
                        1000,5000,1000000)


env<-getPRC20(2)

PCAgenes<-geneMatrix(env$over,env$reg[,1:4],regions,geneList,id="name")

selectGenes(PCAgenes,"ECFC",TRUE)

x<-as.list(org.Hs.egGO2EG)

id<-unlist(sapply(x,function(a) unique(names(a))))

unique(as.numeric(unlist(x)))


terms<-matrix(do.call(c,mapply(rep, names(x),sapply(x,function(a) length(unique(names(a))) ))),ncol=1)

rownames(terms)<-id

rm<-as.list(org.Hs.egGO[mappedkeys(org.Hs.egGO)])

names(rm)<-select(org.Hs.eg.db,names(rm), "SYMBOL","ENTREZID")$SYMBOL           
asl<-lapply(rm,function(y) do.call(rbind,lapply(y,function(x)cbind(x[["GOID"]],x[["Ontology"]],x[["Evidence"]]))))


do.call(rbind,lapply(rm[[3]],function(x)cbind(x[["GOID"]],x[["Ontology"]],x[["Evidence"]])))

do.call(rbind,lapply(rm[["10"]],function(x)cbind(x[["GOID"]],x[["Ontology"]],x[["Evidence"]])))



agenes<-selectGenes(PCAgenes,"Leukemia",TRUE)
lgenes<-length(agenes)

gt<-as.data.frame(do.call(rbind,lapply(agenes,function(i)
    do.call(rbind,lapply(rm[[i]],function(x)cbind(GOID=x[["GOID"]],Ontology=x[["Ontology"]],Evidence=x[["Evidence"]],gene=i))))))


##lapply(list(gt$Ontology=="MF", gt$Ontology=="BP", gt$Ontology=="CC"),sum)



temp<-sapply(as.character(unique(gt[gt$Ontology=="BP"&gt$Evidence=="TAS","GOID"])),function(sec)gt[gt$GOID==sec,"gene"])

lapply(seq(length(temp)),function(i){
    i<<-i
    out<-t(sapply(which(ret[,i]),function(x) cbind(length(temp[[x]]),x)))
    out[order(out[,1],decreasing=TRUE)[1],2]
})



#cs<-makeForkCluster(4)
#ret<-parLapply(cs,combn(temp,2),Vectorize(function(a) 1==length(intersect(a[1],a[2]))/min(length(a[1]),length(a[2])) ))
#stopCluster(cs)
library(compiler)

inner<-Vectorize(function(a,b) 1==length(intersect(unique(a),unique(b)))/min(length(unique(a)),length(unique(b))))

inc<-cmpfun(inner,options=3)
                 
ret<-outer(temp,temp,inc)

temp

system.time(outer(seq(100),seq(1000),inc))
system.time(outer(seq(100),seq(1000),inner))

temp<-as.data.frame(t(sapply(as.character(unique(gt[gt$Ontology=="BP"&gt$Evidence=="TAS","GOID"])),function(sec) c(length(gt[gt$GOID==sec,"gene"]),gosize[sec],numGenes,lgenes))))




numGenes<-length(rm)
gosize<-sapply(as.list(x),length)



p<-dbinom(temp[,1],temp[,4],temp[,2]/temp[,3])
l<-sum(sort(p<1e-2))

rownames(temp[order(p)[1:l],])
select(GO.db,rownames(temp[order(p)[1:l],]),"TERM","GOID")[,2]

select(GO.db,keys(GO.db)[1:10],"DEFINITION","GOID")[,2]



.Internal(unique(x, FALSE, FALSE, min(length(x), nlevels(x) + 1L)))

x<-as.factor(c(1,1,2,1,1))

system.time(lapply(seq(100000),function(x) .Internal(unique(x, FALSE, FALSE, min(length(x), nlevels(x) + 1L)))))

system.time(lapply(seq(100000),function(x) unique(x)))

choices<-c(1,2)

system.time(lapply(seq(100000),function(x) .Internal(pmatch(as.character(x), as.character(choices), 0, TRUE))))

system.time(lapply(seq(100000),function(x) match(x, choices,0L)))

system.time(lapply(seq(100000),function(x){
y<-.Internal(pmatch(as.character(x), as.character(choices), 0, TRUE))
.Internal(unique(y, FALSE, FALSE, min(length(y), length(attr(y,"levels")) + 1L)))
}))


system.time(lapply(seq(100000),function(x){
y<-intersect(x,choices)
}))


system.time(lapply(seq(100000),function(x) length(attr(x,"levels"))))
system.time(lapply(seq(100000),function(x)nlevels(x)))


inner<-Vectorize(function(a,b)
    1==length(intersect(unique(a),unique(b)))/min(length(unique(a)),length(unique(b))))

innerO<-cmpfun(Vectorize(function(a,b){
    y<-.Internal(pmatch(as.character(a), as.character(b), 0, TRUE))
    n<-.Internal(unique(y, FALSE, FALSE, min(length(y), length(attr(y,"levels")) + 1L)))
    la<-length(.Internal(unique(a, FALSE, FALSE, min(length(a), length(attr(a,"levels")) + 1L))))
    lb<-length(.Internal(unique(b, FALSE, FALSE, min(length(b), length(attr(a,"levels")) + 1L))))
    mm<-min(la,lb)
    rem<-length(.Internal(which(n==0)))
    1==(length(n)-rem)/mm#/min(length(unique(a)),length(unique(b)))
    }),options=list(optimize=3,supressALL=TRUE))


innerOb<-cmpfun(Vectorize(function(a,b){
    ca<-as.character(a)
    cb<-as.character(b)
    y<-.Internal(pmatch(ca,cb , 0, TRUE))
    n<-.Internal(unique(y, FALSE, FALSE, length(y)))
    rem<-length(.Internal(which(n==0)))
    d<-unique(c(ca,cb), FALSE, FALSE, length(a)+length(b))
    (length(n)-rem)/length(d)
    }),options=list(optimize=3,supressALL=TRUE))




x<-c(10,20,30,40,50)

times2<-lapply(x, function(i) {print(i);system.time(outer(temp[1:i],temp[1:i],inner))})

times1<-lapply(x, function(i) {print(i);system.time(outer(temp[1:i],temp[1:i],innerO))})


png("~/Desktop/out_ops_comp.png")
y1<-log(do.call(rbind,times2)[,1],10)
y2<-log(do.call(rbind,times1)[,1],10)
#y1<-do.call(rbind,times2)[,1]
#y2<-do.call(rbind,times1)[,1]
y3<-y1
plot(x,y1,col="blue",ylab="LOG10 Time(s)",xlab="Outer Operands",type="l",main="Compare Outer 50X Speedup",ylim=c(min(min(y1),min(y2),min(y3)),max(max(y1),max(y2),max(y3))))
lines(x,y2,col="red")
#lines(x,y3,col="black")
##legend(-1,1.9,c("a","b"),col=c("blue","red"),merge=TRUE)
legend(9, 1, c("old", "new"), col = c(4,2),
       text.col = "black", lty = c(1,  1),
       merge = TRUE, bg = "gray90")
dev.off()

outer(temp[1:10],temp[1:10],inner)

res<-outer(temp,temp,innerO)

subt<-unique(apply(res,1,function(x) names(sort(sapply(which(x),function(i)length(temp[[i]])),decreasing=FALSE))[1]))

temp2<-as.data.frame(t(sapply(subt,function(sec) c(length(gt[gt$GOID==sec,"gene"]),gosize[sec],numGenes,lgenes))))




numGenes<-length(rm)
gosize<-sapply(as.list(x),length)



p<-dbinom(temp2[,1],temp2[,4],temp2[,2]/temp2[,3])
l<-sum(sort(p<1e-1))


nms<-do.call(rbind,lapply(rownames(temp2[order(p)[1:l],]),function(x) paste(temp[[x]],collapse=",")))

terse<-data.frame(GOID=rownames(temp2[order(p)[1:l],]),
                  term=select(GO.db,rownames(temp2[order(p)[1:l],]),"TERM","GOID")[,2],
           genes=nms,
           addColnames(temp2[rownames(temp2[order(p)[1:l],]),],c("a","b","c","d")))



### full run
load("~/Dropbox/UTX-Alex/Paper/genes.RData")

agenes<-selectGenes(PCAgenes,"Erythroid",TRUE)
lgenes<-length(agenes)
numGenes<-length(rm)
x<-as.list(org.Hs.egGO2EG)
gosize<-sapply(as.list(x),length)
gt<-as.data.frame(do.call(rbind,lapply(agenes,function(i)
    do.call(rbind,lapply(rm[[i]],function(x)cbind(GOID=x[["GOID"]],Ontology=x[["Ontology"]],Evidence=x[["Evidence"]],gene=i))))))
##gt$Evidence=="TAS"

temp<-sapply(as.character(unique(gt[gt$Ontology=="BP","GOID"])),function(sec)gt[gt$GOID==sec,"gene"])
temp<-Filter(function(x) length(x)>2,temp)
subt<-names(temp)
temp2<-as.data.frame(t(sapply(subt,function(sec) c(length(gt[gt$GOID==sec,"gene"]),gosize[sec],numGenes,lgenes))))
p<-dbinom(temp2[,1],temp2[,4],temp2[,2]/temp2[,3])
l<-sum(sort(p<0.01))

cn<-names(temp)[order(p)[1:l]]

temp<-lapply(order(p)[1:l],function(x) temp[[x]])
names(temp )<-cn
    
print(length(temp)^2*y2[5]/50^2)

res<-outer(temp,temp,innerO)
subt<-unique(apply(res,1,function(x) names(sort(sapply(which(x),function(i)length(temp[[i]])),decreasing=TRUE))[1]))



temp2<-as.data.frame(t(sapply(subt,function(sec) c(length(gt[gt$GOID==sec,"gene"]),gosize[sec],numGenes,lgenes))))
p<-dbinom(temp2[,1],temp2[,4],temp2[,2]/temp2[,3])
l<-sum(sort(p<0.01))

nms<-do.call(rbind,lapply(rownames(temp2[order(p)[1:l],]),function(x) paste(temp[[x]],collapse=",")))
terse<-data.frame(GOID=rownames(temp2[order(p)[1:l],]),
                  term=select(GO.db,rownames(temp2[order(p)[1:l],]),"TERM","GOID")[,2],
           p=p[order(p)[1:l]],genes=nms,
           addColnames(temp2[rownames(temp2[order(p)[1:l],]),],c("a","b","c","d")))
terse[,2]


tsep<-lapply(terse[,1],function(x) temp[[x]])

res2<-outer(tsep,tsep,innerOb)

plot(((1-res2)%*%prcomp(1-res2)$rotation)[,c(1,2)])

png("~/Desktop/clust_leuk_pca.png",1024,1024)
plot(hclust(dist(1-addNames(CCCA::qn(res2),paste0(substr(terse[,2],1,100)," : 1E-",round(-1*log(terse[,3],10),1),"(",terse[,5],")")),method="manh"),method="ward.D2"), hang=-1,xlab="Ery")
#plot(hclust(dist(1-addNames(CCCA::qn(res2),substr(terse[,2],1,100)),method="manh"),method="ward.D2"), hang=-1)
#plot(hclust(as.dist(1-addNames(res2,substr(terse[,2],1,100))),method="ward.D2"), hang=-1)
dev.off()


match(temp[[1]],)

##match(temp[[1]],geneList$name

listToMatrix<-function(list){
    rn<-unique(unlist(list))
    cn<-names(list)
    lenc<-length(cn)
    lenr<-length(rn)
    out<-matrix(FALSE,lenr,lenc)
    rownames(out)<-rn
    colnames(out)<-cn
    for(i in seq(lenc))
        out[as.character(list[[i]]),i]=TRUE
    out
}

BP<-unique(unlist(do.call(function(x){x$GOID[x$Ontology=="BP"]},list(as.data.frame(t(sapply(rm,function(x) x[[1]])))))))

y<-lapply(BP,function(i) x[[i]])
names(y)<-BP

test<-listToMatrix(y)


d<-dist(test)

library(fastcluster)

plot(hclust(dist(test)),hang=-1)



library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(fastcluster)
x<-as.list(org.Hs.egGO2EG)
rm<-as.list(org.Hs.egGO[mappedkeys(org.Hs.egGO)])
names(rm)<-select(org.Hs.eg.db,names(rm), "SYMBOL","ENTREZID")$SYMBOL
BP<-unique(unlist(do.call(function(x){x$GOID[x$Ontology=="BP"]},list(as.data.frame(t(sapply(rm,function(x) x[[1]])))))))
y<-lapply(BP,function(i) x[[i]])
names(y)<-BP

test<-listToMatrix(y)



td<-t(test)
out<-td %*% t(td)
mag<-rowSums(td)
sim<-out/outer(mag,mag,"*")


D<-diag(colSums(sim))
D1<-diag(1/colSums(sim))

P<-D1%*%sim

sv<-svd(P)

X<-sv$v[,1:25]


e<-eigen(P)



##ret2<-lapply(seq(71,100),function(n){

n<-93
X<-e$vectors[,1:n]
Y<-X/apply(X,1,function(a) sqrt(sum(a^2)))
##d<-dist(Y,"manh")





h<-hclust(d,method="single")
cut<-function(sim,a,b,t)sum(sim[t==a,t==b])
t<-cutree(h,n)
V<-sapply(seq(n),function(a) sum(diag(D)[t==a]))    
C<-outer(seq(n),seq(n),Vectorize(function(a,b)cut(sim,a,b,t)))
sum((C/V)[upper.tri(C/V)])  + sum(e$values[1:n]) - n
##})

#cgap<-cbind(seq(45,100),c(unlist(ret),unlist(ret2)))

#plot(seq(45,100),c(unlist(ret),unlist(ret2)),type="l")

sum((C/V)[upper.tri(C/V,diag=TRUE)])-n  


plot(sim[1:10,1:10])

200-sum(sv$d[1:200]/sv$d[1])

100-sum(e$value[1:100])

sum(sim[t==1,t==2])/sum(diag(D)[t==1])


plot(h,hang=-1)

tree
plot(sort(tmat[1,]))

d<-apply(a,1,function(b) rowSums(abs(a-b)))
        
plot(hclust(as.dist(d),method="ward.D2"),hang=-1)

plot(hclust(as.dist(d),method="average"),hang=-1)

100

out<-matrix(0,10,10)
out[lower.tri(out)]<-apply(combn(10,2),2,function(a)sum(td[a[1],]&td[a[2],]))
diag(out)<-apply(combn(10,1),2,function(a)sum(td[a,]))
out

plot(hclust(as.dist(out),method="ward.D2"),hang=-1)

outer(seq(10),seq(10),)


plot(hclust(dist(test)),hang=-1)


apply(combn(100,2),2,function(a)sum(td[a[1],]&td[a[2],]))

apply(combn(100,2),2,function(a))

library(data.table)

tdb<-as.data.table(td)

apply(combn(100,2),2,function(a)sum(tdb[a[1],]&tdb[a[2],]))



heatmap(sim[1:300,1:300])
