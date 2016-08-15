source("~/r-workspace/dec-functions.r")
source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-variables.r")
source("~/Masters/CCCA/inst/scipts/eboxFrequency.r")
library(msa)

# Functions

freq2SDFull<-function(f1,lm,n=1000){
    randomLocs<-lapply(sapply(lm,length),function(x)lapply(seq(n),function(y)sample(seq(length(fasta5)),x)))
    eboxStatsP<-do.call(abind,append(lapply(randomLocs,function(randomReg) apply(motifAllFrequency(genEboxCombs(),seq(n),motifLoc,randomReg),1,function(x) cbind(mean(x),sqrt(var(x))))),list(along=3) ))
    (f1-eboxStatsP[1,,])/eboxStatsP[2,,]
}

freq2SDQ<-function(f1,eboxStats){
    (f1-eboxStats[1,])/eboxStats[2,]
}


# Init Variabls

over5=read.table(mfn22(pvalues[3]),header=T)
reg5p<-cbind(reg5,HSC5,Other5)
motifLoc<-sapply(genEboxCombs(),grepMotifs,fasta5)
peakLoc<-lapply(seq(4,dim(over5)[2]),Compose(sel2(over5),as.logical,which))
uniqueLocs<-apply(over5[,4:22],1,sum)==1
regLoc5<-lapply(seq(dim(reg5p)[2]),Compose(sel2(reg5p),which))



randomReg<-lapply(seq(1000),function(x)sample(seq(length(fasta5)),2000))
eboxStats<-apply(motifAllFrequency(genEboxCombs(),seq(1000),motifLoc,randomReg),1,function(x) cbind(mean(x),sqrt(var(x))))


allData<-data.frame(count=do.call(rbind,lapply(motifLoc,length)))
write.table(allData,"eboxAll-pvalue=5.matrix", quote=FALSE,row.names=FALSE,col.names = FALSE)

bedEboxUniqueFreq<-motifUniqueFrequency(genEboxCombs(),colnames(over5)[4:dim(over5)[2]],motifLoc,peakLoc,uniqueLocs)
write.table(allData,"ebox-bed-unique-freq-pvalue=5.matrix", quote=FALSE,row.names=FALSE,col.names = FALSE)

bedEboxSD<-freq2SDQ(bedEboxUniqueFreq,eboxStats)
write.table(bedEboxSD,"ebox-bed-unique-sd-1000-pvalue=5.matrix", quote=FALSE,row.names=FALSE,col.names = FALSE)

bedEboxFreq<-motifAllFrequency(genEboxCombs(),colnames(over5)[4:dim(over5)[2]],motifLoc,peakLoc)
write.table(allData,"ebox-bed-freq-pvalue=5.matrix", quote=FALSE,row.names=FALSE,col.names = FALSE)

bedEboxSD<-freq2SDQ(bedEboxFreq,eboxStats)
write.table(bedEboxSD,"ebox-bed-sd-1000-pvalue=5.matrix", quote=FALSE,row.names=FALSE,col.names = FALSE)

pcaEboxFreq<-motifAllFrequency(genEboxCombs(),colnames(reg5p),motifLoc,regLoc)
write.table(allData,"ebox-pca-freq-pvalue=5.matrix", quote=FALSE,row.names=FALSE,col.names = FALSE)

pcaEboxSD<-freq2SDQ(pcaEboxFreq,eboxStats)
write.table(pcaEboxSD,"ebox-pca-sd-1000-pvalue=5.matrix", quote=FALSE,row.names=FALSE,col.names = FALSE)

over20=read.table(mfn22(pvalues[9]),header=T)
motifLoc<-sapply(genEboxCombs(),grepMotifs,fasta20)
peakLoc<-lapply(seq(4,dim(over20)[2]),Compose(sel2(over20),as.logical,which))
uniqueLocs<-apply(over20[,4:22],1,sum)==1
regLoc20<-lapply(seq(dim(reg20)[2]),Compose(sel2(reg20),which))

randomReg<-lapply(seq(1000),function(x)sample(seq(length(fasta20)),2000))
eboxStats<-apply(motifAllFrequency(genEboxCombs(),seq(1000),motifLoc,randomReg),1,function(x) cbind(mean(x),sqrt(var(x))))


allData<-data.frame(count=do.call(rbind,lapply(motifLoc,length)))
write.table(allData,"eboxAll-pvalue=20.matrix", quote=FALSE,row.names=FALSE,col.names = FALSE)

bedEboxUniqueFreq<-motifUniqueFrequency(genEboxCombs(),colnames(over20)[4:dim(over20)[2]],motifLoc,peakLoc,uniqueLocs)
write.table(allData,"ebox-bed-unique-freq-pvalue=20.matrix", quote=FALSE,row.names=FALSE,col.names = FALSE)

bedEboxSD<-freq2SDQ(bedEboxUniqueFreq,eboxStats)
write.table(bedEboxSD,"ebox-bed-unique-sd-1000-pvalue=20.matrix", quote=FALSE,row.names=FALSE,col.names = FALSE)

bedEboxFreq<-motifAllFrequency(genEboxCombs(),colnames(over20)[4:dim(over20)[2]],motifLoc,peakLoc)
write.table(allData,"ebox-bed-freq-pvalue=20.matrix", quote=FALSE,row.names=FALSE,col.names = FALSE)

bedEboxSD<-freq2SDQ(bedEboxFreq,eboxStats)
write.table(bedEboxSD,"ebox-bed-sd-1000-pvalue=20.matrix", quote=FALSE,row.names=FALSE,col.names = FALSE)

pcaEboxFreq<-motifAllFrequency(genEboxCombs(),colnames(reg20),motifLoc,regLoc20)
write.table(allData,"ebox-pca-freq-pvalue=20.matrix", quote=FALSE,row.names=FALSE,col.names = FALSE)

pcaEboxSD<-freq2SDQ(pcaEboxFreq,eboxStats)
write.table(pcaEboxSD,"ebox-pca-sd-1000-pvalue=20.matrix", quote=FALSE,row.names=FALSE,col.names = FALSE)


w=12
dev.new(width=w, height=w*2/5)
par(mfrow=c(2,5),mar = c(.5,.5,2,.5))
for (i in substr(genEboxCombs(),3,4))
    pcaEboxView(pcsO[4],i ,motifLoc,prc20,eboxStats)
#dev.off()

pcaEboxView<-function(pc,eboxc,motifLoc,prc20,eboxStats){
    ebox<-paste("CA",eboxc,"TG",sep="")
    n=50
    l<-floor(length(fasta20)/n)
    s<-floor((length(fasta20)/n-l)*n/2)
    a<-((seq(floor(length(fasta20)/l))-1)*l+1+s)
    tr<-cbind(a[1:length(a)-1],a[2:length(a)]-1)
    or<-apply(prc20$eigenVectors[,],2,order)
    subs<-lapply(seq(dim(tr)[1]),function(x) or[tr[x,1]:tr[x,2],pc])
    subsEboxFreq<-motifAllFrequency(genEboxCombs(),seq(dim(tr)[1]),motifLoc,subs)
    subsEboxSD<-freq2SDQ(subsEboxFreq,eboxStats)
    x=a[1:dim(subsEboxSD)[2]]
    plot(x,subsEboxSD[ebox,],xlab="",ylab="",main=ebox,axes=FALSE, frame.plot=TRUE)
    axis(2,tick=FALSE,labels=FALSE)
    axis(1,tick=FALSE,line=FALSE,labels=FALSE)
    lines(x,rep(0,length(x)),type="l",col="red",lwd=2)
    lines(x,rep(1,length(x)),type="l",col="red",lwd=2,lty=2)
    lines(x,rep(-1,length(x)),type="l",col="red",lwd=2,lty=2)
    lines(x,subsEboxSD[ebox,],type="p",pch=20)
}


interweave<-function(a,b){
    out<-rep("",length(a))
    for (i in seq(length(a))){
        out[(i-1)*2+1]=a[i]
        out[(i-1)*2+2]=b[i]
    }
    return(out)
}

write.table(interweave(paste(">",bed20[,1],":",bed20[,2],"-",bed20[,3],sep=""),as.character(fasta20)),"~/thesis-december/fasta20.matix",col.names=FALSE,row.names=FALSE,quote=FALSE)

writeXStringSet(fasta5,"~/thesis-december/fasta5.fasta")






majority<-function(x){
    opts<-unique(x)
    opts[which.max(sapply(opts,function(y) length(which(y==x))))]    
}

accum<-function(x){
    last<-function(x) x[length(x)]
    rest<-function(x) x[2:lenth(x)]
    first<-function(x) x[1]
    ret<-c(first(x))
    sapply(rest(x),function(x) ret<<-append(ret,last(ret)+x))
    ret
}


olgn<-8
mers<-countNMers("~/thesis-december/fasta5.fasta",prc5$eigenVectors[,1],olgn)
normMers<-countNMers("~/thesis-december/fasta5.fasta",rep(1,length(prc5$eigenVectors[,1])),olgn)[,2]
bounds<-cbind(accum(c(0,normMers)[1:length(normMers)])+1,accum(normMers))
tround<-t(outer(seq(10),normMers,Vectorize(function(x,y)sum(sample(prc5$eigenVectors[,1],x,replace=TRUE)))))
m2<-(mers[,2]-apply(tround,1,mean))/apply(tround,1,Compose(var,sqrt))


size=32
tmer<-as.character(mers[order(m2,decreasing=TRUE)[1:size],1])
top<-do.call(abind,append(lapply(tmer,function(x)rbind(score(pairwiseAlignment(tmer,x, gapOpening = 0, gapExtension = 1,type="global-local")),score(pairwiseAlignment(sapply(tmer,compliment),x, gapOpening = 0, gapExtension = 1,type="global-local")))),list(along=3)))
sel<-(top[1,,]-top[2,,]<0)+1
otop<-addNames(outer(seq(size),seq(size),Vectorize(function(x,y)top[sel[x,y],x,y])),tmer)
hc=hclust(dist(otop))

n=5
subs<-cutree(hc,n)#-cutree(hc,k=n)
mtif<-lapply(seq(n),function(n){    
    a<-apply(outer(which(subs==n),which(subs==n),Vectorize(function(x,y)sel[x,y])),1,majority )    
    if(length(a)>1){
        if(length(which(a==2)>0))
            names(a[a==2])<-compliment(names(a[a==2]))
        msaResults<-msa(names(a),"ClustalW",type="dna")        
        c<-consensusString(msaResults)
        b<-consensusMatrix(msaResults)
        bo<-(b/apply(b,2,sum))[1:4,]        
    }
    else{
    c<-names(a)
    bo<-t(do.call(rbind,lapply(strsplit("GATTAGC","")[[1]],function(x,dna)  as.numeric(x==dna),c("A","C","G","T"))))
}
        rbind(c(paste(">",gsub("-","",c),sep=""),c,"",""),round(t(bo),4))
})

write.table(do.call(rbind,mtif),"~/thesis-december/test.pwm",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
stampT<-stampWrapper("~/thesis-december/test.pwm","TRANSFAC_Fams","test")
stampJ<-stampWrapper("~/thesis-december/test.pwm","JASPAR_Fams","test")

cbind(stampT[,c("motif","genes")],stampJ[,c("genes")])

names(subs[])

A2Rplot(hc,k=4,boxes=FALSE,main=name)

t<-sapply(seq(25),function(n){
subs<-cutree(hc,n)
mean(sapply(seq(n),function(i) sum(dist(otop[subs==i,subs==i])^2)))/
          mean(sapply(seq(n),function(i) sum(dist(otop[subs==i,subs!=i])^2)))})

plot(log(abs(t[2:length(t)]-t[1:(length(t)-1)])))
