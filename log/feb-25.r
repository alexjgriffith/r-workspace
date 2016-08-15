# AGGCGG SCACTG (look at prefered distances between erythroid and hsc
# ACNNAC SMNANCSN (ALL)

source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/jan-functions.r")
source("~/r-workspace/feb-functions.r")
source("~/Masters/CCCA/inst/scipts/eboxFrequency.r")
source("~/r-workspace/preferedDistance.r")

smad="CCAGACA"
nr="AGGTCA"
fox="TTGTT[AC]"
ebox=IUPACtoBase("CANNTG")
gata="GATAA"
mzf1<-IUPACtoBase("ACNNAC")

jaspar<-CCCA::loadPWM("~/Masters/mulcal/inst/data/jaspar_motifs.txt",version="jaspar")

# jasparMotifs[which(jasparMotifs[,2]=="ACNNAC"),]
# mzf1

jasparMotifs<-cbind(unlist(jaspar[,2]),unlist(lapply(jaspar[,3],function(x) PWMtoCons(x))))

env<-getPRC20(2)

mList<-c(smad,fox,ebox,gata,nr,mzf1,"AGGCCG",IUPACtoBase("SCACTG"),IUPACtoBase("SMNANCSN") )

motifs<-c(genEboxCombs(),mList)

cList<-sapply(motifs,compliment)
locationsM<-lapply(motifs,grep,env$fasta)
locationsC<-lapply(cList,grep,env$fasta)


name<-plotFastMot(env,motifs,cList,locationsM,locationsC,13,17)


png("~/Dropbox/UTX-Alex/br-data/ebox/CANNTG-GATA.png")
name<-plotFastMot(env,motifs,cList,locationsM,locationsC,13,14)
dev.off()

for(reg in c("Leukemia","Erythroid","HSC","ECFC","ALL","NONE")){
    png(paste0("~/Dropbox/UTX-Alex/br-data/ebox/Eboxs-GATA-",reg,".png"))
    mplotFastMot(env,motifs,cList,locationsM,locationsC,1:10,14,reg)
    dev.off()
}
    
a<-getMinDistance(env$fasta,consenusIUPAC(motifs[13]),consenusIUPAC(motifs[17]),env$ALL)

png("~/Dropbox/UTX-Alex/br-data/ebox/CANNTG-AGGCCG.png")
qstem(a,paste0(consenusIUPAC(motifs[13]),"_",consenusIUPAC(motifs[17])))
dev.off()

db<-read.table("~/feb-homer/motifAnnotations_20_2-alt.table",header=T)

fm<-unique(unlist(sapply(
    as.character(db$filename[intersect(which(db$size==6), grep("NONE",db$dataset))])
     , function(x) sapply(loadPWM(x)[,3],PWMtoCons))))


## number of peaks in each conditnio
apply(env$reg,2, Compose(which,length))


tts(env$bed[env$Leukemia.alt,][order(abs(env$prc$eigenVectors[env$Leukemia.alt,1]),decreasing=TRUE),])[1:10]

## find loc

leukemic<-as.numeric(rownames(env$bed[env$Leukemia.alt,][order(abs(env$prc$eigenVectors[env$Leukemia.alt,1]),decreasing=FALSE),])[1:20])

erythroid<-as.numeric(rownames(env$bed[env$Erythroid.alt,][order(abs(env$prc$eigenVectors[env$Erythroid.alt,1]),decreasing=FALSE),])[1:20])

hsc<-as.numeric(rownames(env$bed[env$HSC.alt,][order(abs(env$prc$eigenVectors[env$HSC.alt,1]),decreasing=FALSE),])[1:20])

ecfc<-as.numeric(rownames(env$bed[env$ECFC.alt,][order(abs(env$prc$eigenVectors[env$ECFC.alt,1]),decreasing=FALSE),])[1:20])

leukemic.top<-as.numeric(rownames(env$bed[env$Leukemia.alt,][order(abs(env$prc$eigenVectors[env$Leukemia.alt,1]),decreasing=TRUE),])[1:20])

erythroid.top<-as.numeric(rownames(env$bed[env$Erythroid.alt,][order(abs(env$prc$eigenVectors[env$Erythroid.alt,1]),decreasing=TRUE),])[1:20])

hsc.top<-as.numeric(rownames(env$bed[env$HSC.alt,][order(abs(env$prc$eigenVectors[env$HSC.alt,1]),decreasing=TRUE),])[1:20])

ecfc.top<-as.numeric(rownames(env$bed[env$ECFC.alt,][order(abs(env$prc$eigenVectors[env$ECFC.alt,1]),decreasing=TRUE),])[1:20])


leukemic.all<-as.numeric(rownames(env$bed[env$Leukemia.alt,][order(abs(env$prc$eigenVectors[env$Leukemia.alt,1]),decreasing=TRUE),]))

erythroid.all<-as.numeric(rownames(env$bed[env$Erythroid.alt,][order(abs(env$prc$eigenVectors[env$Erythroid.alt,1]),decreasing=TRUE),]))

hsc.all<-as.numeric(rownames(env$bed[env$HSC.alt,][order(abs(env$prc$eigenVectors[env$HSC.alt,1]),decreasing=TRUE),]))

ecfc.all<-as.numeric(rownames(env$bed[env$ECFC.alt,][order(abs(env$prc$eigenVectors[env$ECFC.alt,1]),decreasing=TRUE),]))


specialLocs<-strsplit(c(
    "HDAC7 chr12 48225000 Leukemia",
    "Sox1 chr13 112715000 Leukemia",
    "erg chr21 39846000 Leukemia",
    "gata2 chr3 128200000 Leukemia",
    "Sp1 chr11 47380000 HSC",
    "runx1 chr21 36280000 HSC",
    "gata3 chr10 8105000 ECFC")," ")


specialLocs<-strsplit(c(
    "HDAC7 chr12 48225000 Leukemia",
    "erg chr21 39846000 Leukemia",
    "runx1 chr21 36280000 HSC",
    "runx1 chr21 36280000 Leukemia",
    "cnr2 chr13 24205000 Leukemia",
    "gata3 chr11 35140000 Leukemia",    
    "gata3 chr10 8105000 ECFC")," ")


keyVals<-rev(c(sapply(specialLocs,function(x)
    closestPeak(env$bed,x[2],as.numeric(x[3]),env$reg[,x[4]]))))

keyVals.top<-c(leukemic.top,erythroid.top,hsc.top,ecfc.top)

keyVals.all<-c(leukemic.all,erythroid.all,hsc.all,ecfc.all)

keyVals.bottom<-c(leukemic,erythroid,hsc,ecfc)


data<-env$prc$normData[keyVals,]
png("~/Dropbox/Data/new-figs/heatmap-normData-pick.png")
heatmap(addColnames(data,swapFunM(env$categories)),margin=c(5,0),Rowv=NA,Colv=NA)
dev.off()


data<-env$prc$normData[keyVals.top,]
heatmap(addColnames(data,swapFunM(env$categories)),margin=c(5,0),Rowv=NA,Colv=NA,scale="none")

data<-env$prc$normData[keyVals.bottom,]
heatmap(addColnames(data,swapFunM(env$categories)),margin=c(5,0),Rowv=NA,Colv=NA)

data<-env$prc$normData[keyVals.all,]
png("~/Dropbox/Data/new-figs/heatmap-normData-all.png")
heatmap(addColnames(data,swapFunM(env$categories)),margin=c(5,0),Rowv=NA,Colv=NA)
dev.off()

setRange<-function(values,min=0,max=1){
    (values-min(values))/(max(values)-min(values))
}

par(mfrow=c(2,1),mar=c(3,3,3,3),ann=FALSE)
image(as.matrix(apply(data,1,setRange)),axes=FALSE,xlab="Data Set")
plot(seq(0,1,length.out=22),rep(0,22),ylim=c(1.9,2.1))
text((seq(22)-1)/21,2,labels=swapFunM(env$categories),cex=1,srt=90)



image(as.matrix((0:4)/4))

closestPeak<-function(peakList,chr,start,reg){
    sel<-data.frame(n=1:length(env$bed[,"start"]),start=env$bed[,"start"],diff=abs(env$bed[,"start"]-start),chr=env$bed[,"chr"]==chr)
    selInter<-sel[sel$chr&reg,]
    myStart<-selInter[order(selInter[,"diff"]),"start"][1]
    which(peakList[,"start"]==myStart & peakList[,"chr"]==chr)    
}


## bar plot


rna<-addColnames(read.table("~/Dropbox/Data/rna/TAL1-KD_DESeq2_pvalFC.rnk"),c("gene","fc"))

geneFile<-"/home/agriffith/Dropbox/UTX-Alex/jan/hg19.RefSeqGenes.csv"
geneList<-read.delim(geneFile)
chrom<-as.character(geneList$chrom)
tss<-as.numeric(geneList$txStart)
strand<-geneList$strand
# double check to make sure that levels(strand) > c("-","+")
levels(strand)<-c(-1,1)
strand<-as.numeric(strand)
#statsAndData<-readSplitTable(heightFile)

regions<-genomicRegions(chrom,tss,strand,1000,5000,1000000)

env<-getMotifInfo(addFasta(getPRC20(2)),c(genEboxCombs()))

genes<-geneMatrix(env$over,env$reg,regions,geneList)

gn<-getGenes(genes)


Le<-setdiff(gn$Leukemia,gn$Erythroid)
Ee<-setdiff(gn$Erythroid,gn$Leukemia)



## get cc-ebox locations
### from feb-24.r

gm<-sapply(Le,function(x) which(x==as.character(rna$gene)))

gme<-sapply(Ee,function(x) which(x==as.character(rna$gene)))

boxplot(log(unlist(gm),10),log(unlist(gme),10))

boxplot(x=rna$fc[unlist(gm)])


LeukemiaBed<-as.logical(mergeFun(env$over[,4:dim(env$over)[2]],swapFunB)[,"Leukemia"])

gnnbed<-findMotifDist(LeukemiaBed,"CAGATG","GATAA",LeukemiaBed,-14,-15,-16)
gnnpca<-findMotifDist(env$Leukemia.alt,"CAGATG","GATAA",env$Leukemia.alt,-14,-15,-16)

ccbed<-intersect(union(grep("CACCTG",env$fasta,ignore.case=TRUE),grep("CAGGTG",env$fasta,ignore.case=TRUE)),which(LeukemiaBed))
ccpca<-intersect(union(grep("CACCTG",env$fasta,ignore.case=TRUE),grep("CAGGTG",env$fasta,ignore.case=TRUE)),which(env$Leukemia.alt))



reg<-do.call(cbind,lapply(list(gnnbed=gnnbed,gnnpca=gnnpca,ccbed=ccbed,ccpca=ccpca),function(x)makeLogic(x,length(LeukemiaBed))))



genes2<-geneMatrix(env$over,reg,regions,geneList)

gn<-getGenes(genes2)

gm<-sapply(gn$ccpca,function(x) which(x==as.character(rna$gene)))
gme<-sapply(gn$gnnpca,function(x) which(x==as.character(rna$gene)))

boxplot(unlist(gm),unlist(gme))


plot(sort(unlist(gm)),type="l")

plot(sort(unlist(gme)),type="l")

t.test(unlist(gm),unlist(gme))

addNames(lapply(genEboxCombs(),function(i)
    intersect(eenv$eboxLoc[[i]], which(LeukemiaBed))),genEboxCombs(),list=TRUE)

a<-motifDistances(env$eboxLoc,env$fasta)
neb<-which(apply(do.call(cbind,lapply(a,makeLogic,length(env$reg[,"Leukemia"]))),1,sum)==0)
b<-append(a,addNames(list(neb),c("NO"),list=TRUE))
locEboxPCA<-do.call(cbind,lapply(b,function(x) makeLogic(intersect(x,which(env$reg[,"Leukemia"])),length(env$reg[,"Leukemia"]))))
genes2<-geneMatrix(env$over,locEbox,regions,geneList)

#gn<-getGenes(genes2)


uniquePCA<-rep(TRUE,length(env$reg[,"Leukemia"]))#unlist(apply(genes2,1,sum)==1)

#greatGeneAssoc(env$over[locEbox[,1],c(1,2,3)],regions,geneList)$name2

gme<-sapply(rownames(genes2),function(x) {a=which(x==as.character(rna$gene))
                                          if(length(a)>0)
                                              return( rna$fc[a])
                                          else
                                              return( NA)
                                          }
                                                 )

myrmna<-function(x){Filter(function(x)!is.na(x), x)}

gl<-addNames(lapply(seq(11),function(i) myrmna(gme[ intersect(which(uniquePCA),which(genes2[,i]==1))])),colnames(locEbox),list=TRUE)

which(genes[,2]==1)

#gm<-lapply(gn,function(a)sapply(a,function(x) which(x==as.character(rna$gene))))



boxplot(gl[1],gl[[2]],gl[[3]],gl[[4]],gl[[5]],gl[[6]],gl[[7]],gl[[8]],gl[[9]],gl[[10]],gl[[11]])

boxplot(gl,las=2)

outer(1:11,1:11,Vectorize(function(i,j) t.test(gl[[i]],gl[[j]])$p.value))

t2<-lapply(gl,function(x)log(abs(x),2))

outer(1:11,1:11,Vectorize(function(i,j) t.test(t2[[i]],t2[[j]])$p.value))

 outer(1:5,1:5,Vectorize(function(i,j) t.test(tt[[i]],tt[[j]])$p.value))

### bed version
LeukemiaBed<-as.logical(mergeFun(env$over[,4:dim(env$over)[2]],swapFunB)[,"Leukemia"])

t.test(unlist(lapply(seq(10),function(i) t2[[1]])),t2[[11]])

a<-motifDistances(env$eboxLoc,env$fasta)
neb<-which(apply(do.call(cbind,lapply(a,makeLogic,length(env$reg[,"Leukemia"]))),1,sum)==0)
b<-append(a,addNames(list(neb),c("NO"),list=TRUE))

locEboxBED<-do.call(cbind,lapply(b,function(x) makeLogic(intersect(x,which(LeukemiaBed)),length(env$reg[,"Leukemia"]))))
genes3<-geneMatrix(env$over,locEboxBED,regions,geneList)

#gn<-getGenes(genes2)


uniqueBED<-rep(TRUE,length(env$reg[,"Leukemia"]))#unlist(apply(genes3,1,sum)==1)

#greatGeneAssoc(env$over[locEbox[,1],c(1,2,3)],regions,geneList)$name2

gmeBED<-sapply(rownames(genes3),function(x) {a=which(x==as.character(rna$gene))
                                          if(length(a)>0)
                                              return( rna$fc[a])
                                          else
                                              return( NA)
                                          }
                                                 )

glBED<-addNames(lapply(seq(11),function(i) myrmna(unlist(gmeBED[ intersect(which(uniqueBED),which(genes3[,i]==1))]))),colnames(locEboxBED),list=TRUE)

boxplot(glBED[[1]],gl[[1]],las=2)

boxplot(glBED[[2]],gl[[2]],las=2)

boxplot(glBED[[3]],gl[[3]],las=2)

boxplot(glBED[[4]],gl[[4]],las=2)

boxplot(glBED[[5]],gl[[5]],las=2)

boxplot(glBED[[6]],gl[[6]],las=2)

boxplot(unlist(glBED),unlist(gl))

mapply(function(i,j) t.test(log(abs(i),2),log(abs(j),2))$p.value,gl,glBED)


length(which(env$reg[,"Leukemia"] &LeukemiaBed))/length(which(env$reg[,"Leukemia"]))

length(intersect(names(glBED[[2]]),names(gl[[2]])))/length(names(gl[[2]]))




locEboxALL<-do.call(cbind,lapply(b,function(x) makeLogic(intersect(x,which(LeukemiaBed)),length(env$reg[,"Leukemia"]))))

genes<-geneMatrix(env$over,env$reg,regions,geneList)

#gn<-getGenes(genes2)


uniqueALL<-rep(TRUE,length(env$reg[,"Leukemia"]))#unlist(apply(genes3,1,sum)==1)

#greatGeneAssoc(env$over[locEbox[,1],c(1,2,3)],regions,geneList)$name2

gmeALL<-sapply(rownames(genes3),function(x) {a=which(x==as.character(rna$gene))
                                          if(length(a)>0)
                                              return( rna$fc[a])
                                          else
                                              return( NA)
                                          }
                                                 )

glALL<-addNames(lapply(seq(11),function(i) myrmna(unlist(gmeALL[ intersect(which(uniqueALL),which(genes3[,i]==1))]))),colnames(env$reg),list=TRUE)

boxplot(glALL$ALL,glALL$Leukemia)

boxplot(unlist(rna$fc),glALL$ALL,glALL$Leukemia)

boxplot(lapply(list(glALL$Leukemia,unlist(glBED)), function(x) log(abs(x),2)),las=2)

boxplot(lapply(list(gl,glBED[[11]]), function(x) log(abs(x),2)),las=2)

boxplot(lapply(list(unlist(gl),unlist(glBED),unlist(glALL)), function(x) log(abs(x),2)),las=2)



boxplot(t,las=2)

t<-lapply(list(PCA=unlist(lapply(seq(10),function(i) gl[[i]])),BED=unlist(lapply(seq(10),function(i) glBED[[i]]))), function(x) log(abs(x),2))




locGNN<-cbind(makeLogic(gnnbed,dim(env$reg)[1]),makeLogic(gnnpca,dim(env$reg)[1]),makeLogic(ccbed,dim(env$reg)[1]),makeLogic(ccpca,dim(env$reg)[1]))

genesGNN<-geneMatrix(env$over,locGNN,regions,geneList)

uniqueGNN<-rep(TRUE,length(env$reg[,"Leukemia"]))

gmeGNN<-sapply(rownames(genesGNN),function(x) {a=which(x==as.character(rna$gene))
                                          if(length(a)>0)
                                              return( rna$fc[a])
                                          else
                                              return( NA)
                                          }
                                                 )

glGNN<-addNames(lapply(seq(4),function(i) myrmna(unlist(gmeGNN[ intersect(which(uniqueGNN),which(genesGNN[,i]==1))]))),colnames(locGNN),list=TRUE)

boxplot(addNames(lapply(glGNN,function(x)log(abs(x),2)),c(paste0("<",seq(-7,0),"sd"),"ALL"),list=TRUE))

boxplot(glGNN,las=2)

t1<-lapply(glGNN,function(x)log(abs(x),2))
tt<-addNames(append(t1,list(log(abs(rna$fc),2))),c("G-GA BED","G-GA PCA","CC BED","CC PCA","ALL Genes"),list=TRUE)
png("~/Dropbox/UTX-Alex/br-data/box-plot-gata-ebox.png")
boxplot(tt,las=2)
dev.off()

glGNN<-slowP(env,cbind(
    normalize(env$prc$eigenVectors[,1])>6,
    normalize(env$prc$eigenVectors[,1])>5,
    normalize(env$prc$eigenVectors[,1])>4,
    normalize(env$prc$eigenVectors[,1])>3,
    normalize(env$prc$eigenVectors[,1])>2,
    normalize(env$prc$eigenVectors[,1])>1,
    normalize(env$prc$eigenVectors[,1])>0,
    normalize(env$prc$eigenVectors[,1])>-1,
    normalize(env$prc$eigenVectors[,1])>-2))


glGNN<-slowP(env,cbind(
    normalize(env$prc$eigenVectors[,1])<(-7),
    normalize(env$prc$eigenVectors[,1])<(-6),
    normalize(env$prc$eigenVectors[,1])<(-5),
    normalize(env$prc$eigenVectors[,1])<(-4),
    normalize(env$prc$eigenVectors[,1])<(-3),
    normalize(env$prc$eigenVectors[,1])<(-2),
    normalize(env$prc$eigenVectors[,1])<(-1),
    normalize(env$prc$eigenVectors[,1])<(-0),
    TRUE
    )             
             )

png("~/Dropbox/UTX-Alex/br-data/boxplot-PC1-sds.png")
boxplot(addNames(lapply(glGNN,function(x)log(abs(x),2)),c(paste0("<",seq(-7,0),"sd"),"ALL"),list=TRUE))
dev.off()

png("~/Dropbox/UTX-Alex/br-data/boxplot-eboxs-bed-sds.png")
boxplot(lapply(glBED,function(x)log(abs(x),2)),las=2)
dev.off()

slowP<-function(env,locGNN){
#locGNN<-cbind(makeLogic(gnnbed,dim(env$reg)[1]),makeLogic(gnnpca,dim(env$reg)[1]),makeLogic(ccbed,dim(env$reg)[1]),makeLogic(ccpca,dim(env$reg)[1]))
genesGNN<-geneMatrix(env$over,locGNN,regions,geneList)
uniqueGNN<-rep(TRUE,length(env$reg[,"Leukemia"]))
gmeGNN<-sapply(rownames(genesGNN),function(x) {a=which(x==as.character(rna$gene))
                                          if(length(a)>0)
                                              return( rna$fc[a])
                                          else
                                              return( NA)
                                          }
                                                 )
addNames(lapply(1:dim(locGNN)[2],function(i) myrmna(unlist(gmeGNN[ intersect(which(uniqueGNN),which(genesGNN[,i]==1))]))),colnames(locGNN),list=TRUE)
}



t<-lapply(list(glALL$ALL,unlist(rna$fc)), function(x) log(abs(x),2))
boxplot(t,las=2)

png("~/Dropbox/UTX-Alex/br-data/pvalue-PC1-sds.png")
plot(seq(-7,0),sapply(seq(8),function(i) t.test(t1[[i]],t[[2]])$p.value),cex=1.25)
lines(c(-8,1),c(0.05,0.05))
dev.off()

boxplot(append(t1,t[2]))
