library(CCCA)
library(functional)
library(parallel)

#####
createZip<-function(x)
    Map(function(i,j)cbind(x[i],x[j]),seq(1,length(x),2),seq(2,length(x),2))
getSwapCats<-function(cats,names){    
    swCat<-cats
    names(swCat)<-names
    function(x){as.character(unlist(lapply(x,function(x) swCat[x])))}
}
stringToSwap<-function(x)do.call(getSwapCats,splitZip(createZip(strsplit(x," ")[[1]])))
splitZip<-function(inList)
    list(sapply(inList,"[",1),sapply(inList,"[",2))
getSwapCats<-function(cats,names){    
    swCat<-cats
    names(swCat)<-names
    function(x){as.character(unlist(lapply(x,function(x) swCat[x])))}
}
simpleSwapFun<-function(char){
stringToSwap(paste (rev(strsplit(char, " ")[[1]]),collapse=" "))
}
swapFunB<-stringToSwap(paste (rev(strsplit("rpmi Leukemia tall_p3 Leukemia tall_p2 Leukemia tall_p1 Leukemia tall_p2_1 Leukemia tall_p2_2 Leukemia tall_p3_1 Leukemia tall_p3_2 Leukemia jurk_sandar_1 Leukemia jurk_sandar Leukemia jurk Leukemia rpmi_1 Leukemia rpmi_2 Leukemia cem_1 Leukemia cem_2 Leukemia cem Leukemia ecfc-tsa ECFC ecfc.tsa ECFC ecfc ECFC meka MEKA cd133 HSC cd34 HSC cd34_new HSC eryt Erythroid eryt_f Erythroid eryt_a Erythroid k562 Erythroid k562_1 Erythroid k562_2 Erythroid", " ")[[1]]),collapse=" "))
swapFun<-stringToSwap(paste (rev(strsplit("rpmi RPMI tall_p3 Prima5 tall_p2 Prima2 tall_p1 Prima5 tall_p2_1 Prima2 tall_p2_2 Prima2 tall_p3_1 Prima5 tall_p3_2 Prima5 jurk_sandar_1 Jurakt jurk_sandar Jurkat jurk Jurkat rpmi_1 RPMI rpmi_2 RPMI cem_1 CEM cem_2 CEM cem CEM ecfc-tsa ECFC ecfc.tsa ECFC ecfc ECFC meka MEKA cd133 HSC cd34 HSC cd34_new HSC eryt Erythroid eryt_f Erythroid eryt_a Erythroid k562 K562 k562_1 K562 k562_2 K562", " ")[[1]]),collapse=" "))
swapFunC<-stringToSwap(paste (rev(strsplit("Leukemia blue Erythroid red HSC orange ECFC green MEKA brown", " ")[[1]]),collapse=" "))
#####


catFile<-"~/Dropbox/UTX-Alex/jan/catagories"
#catFile<-"/mnt/brand01-00/mbrand_analysis/peaks/jan/catagories"
readCategories<-Compose(read.table,unlist,as.character)
categories<-readCategories(catFile)
seql<-Compose(length,seq)
rawDataFiles<-function(categories)paste("/mnt/brand01-00/mbrand_analysis/data_sets/",categories,"/",categories,"_unique_nodupes.bed",sep="")
#makePeakFiles<-function(categories)paste("~/.peaktemp/",categories,"~combined_mock_peaks.xls",sep="")
makePeakFiles<-function(categories)paste("~/Dropbox/Data/august_peaks/",categories,"~combined_mock_peaks.xls",sep="")
all(sapply(makePeakFiles(categories),file.exists))
#AFS_22x22_0<-makeAFS(makePeakFiles(categories),categories)
wrapAFS<-function(makePeakFiles,categories){
    function(pvalue=0)
        makeAFS(makePeakFiles(categories),categories,pValue=pvalue)
}
renameColumn<-function(x,colName,newColName)
    {
        UseMethod("renameColumn",x)
    }
renameColumn.data.default<-function(x,colName,newColName){
    warning("renameColumn only works for data.frame currently.")
    x
}
renameColumn.data.frame<-function(x,colName,newColName){
    whichName<-(colnames(x)==colName)
    if(! any(whichName))
        warning(colName," was not found in colnames(x)!")
    else{
    newCols<-replace(colnames(x),whichName,newColName)
    colnames(x)<-newCols}
    x
}
categoriesB<-c("cd34_new","eryt","jurk","cem_1")
fullAFS<-wrapAFS(makePeakFiles,categories)
pilotAFS<-wrapAFS(makePeakFiles,categoriesB)
#test<-fullAFS()
#sdata<-hg19Sort(renameColumn(test,"chr","chro")[,1:3])

#save.image("~/AFS.RData")
#install_github("alexjgriffith/CCCA")
#load("~/AFS.RData")
#library(CCCA)

#score<-pileUp(sdata,rawDataFiles(categories),22)
#colnames(score)<-categories

#cor_22x22_0<-cor(score)

n=30
cs<-makeForkCluster(n,renice=0)
pvalues<-seq(from=0,by=2.5,length=n)
ret<-parallel::parLapply(cs,pvalues,fullAFS)
stopCluster(cs)

for(i in seql(ret)){
   fileName<-paste("~/thesis-november/","22x22-pvalue=",pvalues[i],"_heights.matrix",sep="")
   matrixName<-paste("~/thesis-november/","22x22-pvalue=",pvalues[i],".matrix",sep="")   
   print(fileName)
   sdata<-hg19Sort(renameColumn(ret[[i]],"chr","chro"))
   score<-pileUp(sdata[,1:3],rawDataFiles(categories),20)   
   write.table(score,fileName,quote=FALSE,row.names=FALSE)
   write.table(sdata,matrixName,quote=FALSE,row.names=FALSE)
   cmd<-paste("scp",fileName,paste("griffita@ogic.ca:/data/websites/dropbox.ogic.ca/mbrandlab/alex_analysis/thesis/exdata/","22x22-pvalue=",pvalues[i],"_heights.matrix",sep=""))
   system(cmd)
   cmd<-paste("scp",matrixName,paste("griffita@ogic.ca:/data/websites/dropbox.ogic.ca/mbrandlab/alex_analysis/thesis/exdata/","22x22-pvalue=",pvalues[i],".matrix",sep=""))
   system(cmd)
}


n=30
pvalues<-seq(from=0,by=2.5,length=n)
getCor<-function(i){
    file<-paste("~/thesis-november/22x22-pvalue=",i,"_heights.matrix",sep="")
test<-read.table(file,header=T)
colnames(test)<-categories
cor(test)
}

getOve<-function(i,cats=NULL,fun=overlap){
    file<-paste("~/thesis-november/22x22-pvalue=",i,".matrix",sep="")
    test<-read.table(file,header=T)
    test<-test[4:length(test)]
    depth<-dim(test)[1]
    if(is.null(cats))
        cats<-colnames(test)
    out<-matrix(numeric(depth*length(cats)),ncol=length(cats))
    colnames(out)<-cats
    for(i in cats){
        if(i %in% colnames(test))
            out[,i]=test[,i]
    }
    fun(out)
    #out
}


overlap<-function (f) 
{
    permutations<-function(n)
        do.call(rbind,lapply(seq(n),function(x) cbind(seq(n),x)))

    overlapNorm <- function(x, f) {
        a <- f[,x[1]]
        b <- f[,x[2]]
        nu<-sum(a & b)
        de<-sum(a | b)
        if(de==0)
            return (1)
        else
            return (nu/de)
    }    
    width <- dim(f)[2]
    m <- matrix(apply( apply(permutations(seq(width)),1, as.numeric), 2, overlapNorm, f ), ncol = width)
    colnames(m) <- colnames(f)
    rownames(m) <- colnames(f)
    m
}



correl<-getCor(pvalues[1])
for(i in seq(2,n))
    correl<-abind(correl,getCor(pvalues[i]),along=3)

overl<-getOve(pvalues[1])
cats<-colnames(overl)
for(i in seq(2,n))
    overl<-abind(overl,getOve(pvalues[i],cats),along=3)

c2<-sapply(seq(22),function(i)lapply(seq(22),function(j)cor(pvalues,correl[i,j,])))
diag(c2)<-1
colnames(c2)<-swapFun(categories)
rownames(c2)<-swapFun(categories)
c3<-matrix(as.numeric(c2),ncol=22)
colnames(c3)<-swapFun(categories)
rownames(c3)<-swapFun(categories)
heatmap(c3)

plot(hclust(1-as.dist(c2),"complete"),hang=-1)

library(grid)
library(ggplot2)

#p1<-heightHist(t0,which.max(t0)+min(h)-1)
#p2<-locHist(h,paste(m1,m2,sep="-"))
plotV<-function(a,b,v,ylab,title){    
    df<-data.frame(pvalues=pvalues,corelation=v[a,b,])
    best<-coef(lm(corelation ~pvalues , data=df))
    print(best)
    ggplot(df,aes(x=pvalues,y=corelation))+xlab("Cut Off (10^(-pvalue))")+ylab(ylab)+geom_point(size=3,shape=1)+ggtitle(title)+
        geom_abline(intercept=best[1],slope=best[2])+
        theme_bw() +
            theme(axis.line = element_line(colour = "black"),
                  axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title=element_text(size=14,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
}


p1<-plotV("eryt","jurk",correl,"Correlation","Jurkat and Erythroid")
p2<-plotV("eryt","jurk",overl,"Overlap","Jurkat and Erythroid")
p3<-plotV("cem_1","jurk",correl,"Correlation","Jurkat and CEM")
p4<-plotV("cem_1","jurk",overl,"Overlap","Jurkat and CEM")
#svg("~/Dropbox/overlap-correlation.svg")
pushViewport(viewport(layout=grid.layout(2,2)))    
print(p1,vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(p2,vp=viewport(layout.pos.row=1,layout.pos.col=2))
print(p3,vp=viewport(layout.pos.row=2,layout.pos.col=1))
print(p4,vp=viewport(layout.pos.row=2,layout.pos.col=2))
#dev.off()

data<-list()
for(i in pvalues){
    file<-paste("~/Dropbox/Data/22x22-pvalue=",i,"_heights.matrix",sep="")
    test<-read.table(file,header=T)
    colnames(test)<-categories
    data<-append(data,list(test))
}

matr<-array(numeric(length(categories)^2*length(pvalues)),c(length(categories),length(categories),length(pvalues)))
            
for(i in seql(data)){
    pcs<-pca(data[[i]])
    temp<-t(pcs$normData)%*%pcs$eigenVectors
    matr[,,i]<-temp
}

rownames(matr)<-categories
colnames(matr)<-paste("PC",seql(categories),sep="")

swapFunC<-stringToSwap(paste (rev(strsplit("Leukemia blue Erythroid red HSC orange ECFC green MEKA brown", " ")[[1]]),collapse=" "))

cols<-swapFunB(categories)




shiftdf<-function(category,x,y ,df=NULL){
    dft<-data.frame(categories=category, x=x, y=y)
    if(is.null(df))
        return(dft)
    else
        return(rbind(df,dft))
}
shiftCols<-function(x,y,categories,shiftdf){
    nx<-diff(range(x))/2
    ny<-diff(range(y))/2
    cats<-as.character(unlist(categories))
    cat<-as.character(unlist(shiftdf$categories))
    for(i in seq(dim(shiftdf)[1])){
        if(cat[i] %in% cats){
            print(shiftdf$y*ny)
           x[which(cat[i]==cats)]= x[which(cat[i]==cats)]+shiftdf$x[i]*nx
           y[which(cat[i]==cats)]= y[which(cat[i]==cats)]+shiftdf$y[i]*ny
       }
    }
    data.frame(x=x,y=y,categories=cats)
}

### shift for 22x22-pca
sdf<-shiftdf("ecfc-tsa",0,-0.05)
sdf<-shiftdf("cd34",-0.12,0,sdf)
sdf<-shiftdf("cd133",-0.12,0,sdf)
sdf<-shiftdf("cd34_new",0,-.05,sdf)
sdf<-shiftdf("meka",0,-.05,sdf)
sdf<-shiftdf("eryt_f",0,.05,sdf)
sdf<-shiftdf("eryt",0,.05,sdf)
sdf<-shiftdf("eryt_a",-.17,0,sdf)
sdf<-shiftdf("k562_1",-.15,0,sdf)
sdf<-shiftdf("k562_2",-.05,-.05,sdf)
sdf<-shiftdf("jurk_sandar_1",.1,0.05,sdf)
sdf<-shiftdf("rpmi_1",.1,0.05,sdf)
sdf<-shiftdf("rpmi_2",.15,.02,sdf)
sdf<-shiftdf("tall_p1",.01,-.03,sdf)
sdf<-shiftdf("tall_p2_1",.15,.02,sdf)
sdf<-shiftdf("tall_p2_2",.2,0,sdf)
sdf<-shiftdf("jurk",-.02,.03,sdf)
sdf<-shiftdf("tall_p3_1",-0.1,-0.05,sdf)
sdf<-shiftdf("tall_p3_2",0.15,-0.05,sdf)
sdf<-shiftdf("cem_1",-0.11,.04,sdf)
sdf<-shiftdf("cem_2",.15,-.01,sdf)
sdf<-shiftdf("jurk_sandar",-.02,-.03,sdf)
df<-data.frame(x=matr[,"PC1",1],y=matr[,"PC3",1],categories=swapFun(categories),Conditions=cols)
postext<-shiftCols(df$x,df$y,categories,sdf)
p<-ggplot(df,aes(x=x,y=y,col=Conditions,label=categories))+ylab("PC3")+xlab("PC1")+geom_point(size=3,shape=1)+theme_bw()+geom_text(x=postext$x,y=postext$y,show_guide=F,size=4)+scale_color_manual(values=c("green","red" , "orange","blue","brown"))
                                        #scale_color_manual(values=c("green","red" , "orange","blue","brown"))
ggsave("~/Dropbox/pca-1-3-labled.png")
#plot(df[,c("x","y")])


##### looking at how clusters change as the p-value is shifted
a=1
b=4
c=8
labels<-paste("PC",seq(22),sep="")
df<-data.frame(x=matr[,labels[a],c],y=matr[,labels[b],c],categories=swapFun(categories),Conditions=cols)
ggplot(df,aes(x=x,y=y,col=Conditions,label=categories))+ylab(labels[b])+xlab(labels[a])+geom_point(size=3,shape=1)+theme_bw()+ scale_color_manual(values=c("green","red" , "orange","blue","brown"))+geom_text(show_guide=F)

ggsave(paste("~/Dropbox/pca-",a,"-",b,"-",pvalues[c],".png"))


n=30
test<-applyPeakPartitions(pca(data[[30]])$eigenVector,list(list(1,"lt",1),list(1,"gt",1),list(4,"lt",1),list(4,"gt",1)))
apply(test,2,sum)

plot(pca(data[[30]])$eigenVector[,c(1,4)])

#### Contribution Plots
tags<-list()
tags<-append(tags,list(getOve(0,fun=function(x) x)))
ctemp<-colnames(tags[[1]])
colnames(tags[[1]])<-categories
for(i in pvalues[2:30]){
    test<-getOve(i,cats=ctemp,fun=function(x) x)
    print(dim(test))
    colnames(test)<-categories
    tags<-append(tags,list(test))
}



test<-applyPeakPartitions(pca(data[[n]])$eigenVector,list(list(1,"lt",1),list(1,"gt",1),list(4,"lt",1),list(4,"gt",1)))
head(tags[[n]][test[,1],])



#stackedContrib(pcs$eigenVectors[,2],"contrib2",n=3,tags=tags,swapFun=swapFun)


#apply(mergeFun(tags[[n]],swapFun),2,sum)
#simpleSwapFun("green red orange blue brown")
#temp<-mergeFun(tags[[n]],swapFun)
n=10
pcs<-pca(data[[n]])$eigenVector[,4]
stackedContrib(pcs,command="contrib2",tags=mergeFun(tags[[n]],swapFunB),sum=TRUE,n=6,colourOveride=swapFunC)

mergeFun<-function(ma,swapFun){
    newCols<-swapFunB(colnames(ma))
    #print(newCols)
    unc<-unique(newCols)
    outL<-c()
    for(i in unc){
        pos<-which(i==newCols)
        #print(pos)
        #print(data.frame(c=colnames(ma)[pos],v=do.call(rbind,lapply(pos,function(x) sum(ma[,x])))))
        if(length(pos)>1)
            temp<-orM(ma[,pos])
        else{
            temp<-ma[,pos]
            print(sum(temp))
        }
        outL<-cbind(outL,temp)}
    
    colnames(outL)<-unc
    outL
}

#### Coloured histograms
library(ggplot2)
library(zoo)
library(ggdendro)
library(plyr)
cut<-4
hc<-hclust(as.dist(1-cor(qn(data[[1]]))),"ave")
dendr<-dendro_data()
clust<-cutree(hc,k=cut)
clust.df<-data.frame(label=names(clust),cluster=clust)

height<-


for(i in seql(ret)){
   fileName<-paste("~/thesis-november/","4x4-pvalue=",pvalues[i],".matrix",sep="")   
   print(fileName)
   write.table(dim(ret[[i]]),fileName,quote=FALSE,row.names=FALSE)
}

for(i in seql(ret)){
   fileName<-paste("~/thesis-november/","22x22-pvalue=",pvalues[i],"_heights.matrix",sep="")
   print(fileName)
   pups<-pileUp(ret[[i]],rawDataFiles(categories),22)
   write.table(dim(ret[[i]]),fileName,quote=FALSE,row.names=FALSE)
}
   

# off cluster
pvalues<-seq(from=0,by=2.5,length=n)
for(i in seql(ret)){
   fileName<-paste("~/thesis-november/","22x22-pvalue=",pvalues[i],".matrix",sep="")   
   print(fileName)
   write.table(dim(ret[[i]]),fileName,quote=FALSE,row.names=FALSE)
}

sdata<-hg19Sort(renameColumn(test,"chr","chro")[,1:3])

score<-pileUp(sdata,rawDataFiles(categories[1:2]))

data.frame.cha
replace(colnames(test),colnames(test)=="chr","chro")
data<-read.table("~/Dropbox/Data/thesis-november/4x4-0-overlap-matrix.data",header=TRUE)
