source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/jan-functions.r")
source("~/r-workspace/feb-functions.r")
source("~/r-workspace/feb-variables.r")


tts<-function(x){
    paste0(x[,1],":",x[,2]-2000,"-",x[,3]+2000)
}


env20<-getPRC20(2);

write.table(tts(env20$bed[env20$reg[,"Erythroid"],]),"~/Dropbox/UTX-Alex/br-data/eryt_sd_2.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(tts(env20$bed[env20$reg[,"Leukemia"],]),"~/Dropbox/UTX-Alex/br-data/jurk_sd_2.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(tts(env20$bed[env20$reg[,"ECFC"],]),"~/Dropbox/UTX-Alex/br-data/ecfc_sd_2.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(tts(env20$bed[env20$reg[,"HSC"],]),"~/Dropbox/UTX-Alex/br-data/stem_sd_2.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

env20<-addFasta(env20)

smad="CCAGACA"
nr="AGGTCA"
fox="TTGTT[AC]"
ebox=IUPACtoBase("CANNTG")

getMotifInfo

length(intersect(which(env20$reg[,"ECFC"])))

grep(TTGTTAC,env20$fasta[env20$reg[,"ECFC"],])

grep("TTGTT",env20$fasta[env20$reg[,"ECFC"],])


selmin<-function(x,y){
    a<-x-y
    i<-which.min(abs(a))
    a[i]
}


m1<-fox
m2<-ebox

slowMot<-function(m1,m2,dr)

stem <- function(x,y,pch=16,linecol=1,clinecol=1,...){
if (missing(y)){
    y = x
    x = 1:length(x) }
    plot(x,y,pch=pch,...)
    for (i in 1:length(x)){
       lines(c(x[i],x[i]), c(0,y[i]),col=linecol)
    }
    lines(c(x[1]-2,x[length(x)]+2), c(0,0),col=clinecol)
}

dr<-env20$fasta[locs,]

## this stuff is broken don't trust it (slowMot and plotSlowMot)
slowMot<-function(m1,m2,dr,reg,title){
    locs<-intersect(intersect(
        union(grep(m1,dr,ignore.case=TRUE),
              grep(compliment(m1),dr,ignore.case=TRUE)),
        union(grep(m2,dr,ignore.case=TRUE),
              grep(compliment(m2),dr,ignore.case=TRUE)
            )),
                    which(reg)
                    )
a<-mapply(function(x,y) c(x,300-y), gregexpr(m1, dr[locs,],ignore.case=TRUE),gregexpr(compliment(m1), dr[locs,],ignore.case=TRUE))
b<-mapply(function(x,y) c(x,300-y), gregexpr(m2, dr[locs,],ignore.case=TRUE),gregexpr(compliment(m2), dr[locs,],ignore.case=TRUE))
dists<-mapply(function(x,y){selmin(c(outer(x,y,Vectorize(function(x,y) x-y))),0)},a,b)
x<-do.call(seq,as.list(range(dists)))
stem(x,getHeights(dists),xlim=c(-32,32),main=title)
}

plotSlowMot<-function(m1,m2,env20,cats=c("ECFC","Erythroid","HSC","Leukemia")){
par(mfrow=c(2,2),mar=c(2,2,3,2), oma=c(0, 0, 2, 0))
slowMot(m1,m2,env20$fasta,env20$reg[,"ECFC"],"ECFC")
slowMot(m1,m2,env20$fasta,env20$reg[,"Erythroid"],"Erythroid")
slowMot(m1,m2,env20$fasta,env20$reg[,"HSC"],"HSC")
slowMot(m1,m2,env20$fasta,env20$reg[,"Leukemia"],"Leukemia")
title(paste(consenusIUPAC(m1),"-",consenusIUPAC(m2)),outer=TRUE)
}

plotSlowMot("GATAA",ebox,env20)

plot(getHeights(apply(cbind(sapply(gregexpr(m1, env20$fasta[locs,],ignore.case=TRUE),selmin,150),
      sapply(gregexpr(compliment(m1), env20$fasta[locs,],ignore.case=TRUE),selmin,150  )),1,selmin,0)))


env20<-addFasta(env20)

#' values<-loadHeightFile(heightFile)
#' data<-values$data
#' reg<-ascore(data,1,"top",3)
#' test<-readDNAStringSet(fastaFile,use.names = TRUE)
#' motifs<-as.matrix(read.table("data/normal_not_abnormal_motifs"))
#' mList<-unlist(lapply(c(motifs,addmotifs),IUPACtoBase))
#' cList<-unlist(lapply(lapply(c(motifs,"CGNNGC"),IUPACtoBase),compliment))
#' locationsM<-lapply(mList,grep,test)
#' locationsC<-lapply(cList,grep,test)
#' l<-length(cList)
#' motifHist(mList,cList,locationsM,locationsC,4,l,reg)


mList<-c(fox,ebox,"GATAA",nr)
cList<-sapply(mList,compliment)
locationsM<-lapply(mList,grep,env20$fasta)
locationsC<-lapply(cList,grep,env20$fasta)

plotFastMot(env20,mList,cList,locationsM,locationsC,2,3)

plotFastMot<-function(env20,mList,cList,locationsM,locationsC,n1,n2,reg=c("ECFC","Erythroid","HSC","Leukemia")){
par(mfrow=c(2,2),mar=c(2,2,3,2), oma=c(0, 0, 2, 0))
makeStem(motifHist(env20$fasta,mList,cList,locationsM,locationsC,n1,n2,env20$reg[,reg[1]]),main=reg[1])
makeStem(motifHist(env20$fasta,mList,cList,locationsM,locationsC,n1,n2,env20$reg[,reg[2]]),main=reg[2])
makeStem(motifHist(env20$fasta,mList,cList,locationsM,locationsC,n1,n2,env20$reg[,reg[3]]),main=reg[3])
makeStem(motifHist(env20$fasta,mList,cList,locationsM,locationsC,n1,n2,env20$reg[,reg[4]]),main=reg[4])
title(paste(consenusIUPAC(mList[n1]),"-",consenusIUPAC(mList[n2])),outer=TRUE)
paste0(consenusIUPAC(mList[n1]),"-",consenusIUPAC(mList[n2]))
}


makeStem(h)

makeStem<-function(y,xlim=c(-32,32),...){
    yn<-y[y>xlim[1]&y<xlim[2]]
    yt<-getHeights(yn);
    
    stem(do.call(seq,as.list(range(yn))),yt,xlim=xlim,...)
    #stem(xlim,yt,xlim=xlim,...)
}


source("~/Masters/CCCA/inst/scipts/eboxFrequency.r")

motifs<-c(genEboxCombs(),mList)
cList<-sapply(motifs,compliment)
locationsM<-lapply(motifs,grep,env$fasta)
locationsC<-lapply(cList,grep,env$fasta)

i=11
name<-plotFastMot(env,motifs,cList,locationsM,locationsC,7,13)

png(paste0("~/Dropbox/UTX-Alex/br-data/ebox/",name,".png"))
name<-plotFastMot(env20,motifs,cList,locationsM,locationsC,i,13)
dev.off()

chr9:98205264-98270831
a<-env20$prc$normData
ls<-which(env20$bed$chr=="chr9" &
env20$bed$start>98205264 &
    env20$bed$end<98270831)


env20$prc$normData[ls[1],]

norm<-qn(log(env20$heights+1,2))

cat(t(norm)[,ls[1]])

cor(norm[ls[1],],env20$prc$normData[ls[1],])


CCL<-intersect(union(locationsC$CACCTG,locationsM$CACCTG),which(env20$reg[,"Leukemia"]))

GAL<-intersect(union(locationsC$CAGATG,locationsM$CAGATG),which(env20$reg[,"Leukemia"]))

GAE<-intersect(union(locationsC$CAGATG,locationsM$CAGATG),which(env20$reg[,"Erythroid"]))

NN<-locationsC$`CA[ACGT][ACGT]TG`

CC<-union(locationsC$CACCTG,locationsM$CACCTG)




G<-union(locationsC$GATAA,locationsM$GATAA)

GAEG<-intersect(GAE,G)

GALG<-intersect(GAL,G)

geneGALG<-unique(greatGeneAssoc(env20$bed[GALG,],regions,geneList)$name2)
write.table(geneGALG,"~/Dropbox/UTX-Alex/br-data/genes/feb13_GALG.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

geneGAEG<-unique(greatGeneAssoc(env20$bed[GAEG,],regions,geneList)$name2)
write.table(geneGAEG,"~/Dropbox/UTX-Alex/br-data/genes/feb13_GAEG.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)



length(intersect(geneGALG,geneGAEG))

geneGAL<-unique(greatGeneAssoc(env20$bed[GAL,],regions,geneList)$name2)
write.table(geneGAL,"~/Dropbox/UTX-Alex/br-data/genes/feb13_GAL.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

geneGAE<-unique(greatGeneAssoc(env20$bed[GAE,],regions,geneList)$name2)
write.table(geneGAE,"~/Dropbox/UTX-Alex/br-data/genes/feb13_GAE.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)


length(intersect(geneGAL,geneGAE))

geneCCL<-unique(greatGeneAssoc(env20$bed[CCL,],regions,geneList)$name2)
write.table(geneCCL,"~/Dropbox/UTX-Alex/br-data/genes/feb13_CCL.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)


write.table(tts(env20$bed[CCL,]),"~/Dropbox/UTX-Alex/br-data/CC_leukemia_sd_2_2.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)

printBB(env20$bed[CCL,],"~/Dropbox/CC_leukemia_sd_2_new.bb")

gata3<-addColnames(loadBedFile("~/Dropbox/UTX-Alex/br-data/jurk-gata/combined~jurk_gata3_peaks.bed"),c("chr","start","end"))


bedGAL<-env20$bed[GAL,]

bedCCL<-env20$bed[CCL,]

bedCC<-env20$bed[CC,]

bedNN<-env20$bed[NN,]

bedNNE<-env20$bed[intersect(NN,which(env20$reg[,"Leukemia"])),]

bedOverlaps<-function(bed1,bed2){
    #bed1$chr==bed2$chr & bed1&start>bed2$end & bed1&end<bed2$start}
    chrs=as.character(unique(bed1$chr))
    lapply(chrs, function(c){length(which(bed1$chr==))},bed2)
}





bedOverlaps<-function(a,b){
    pu<-sortDataFrame(rbind(cbind(a,name="a",summit=(a$end+a$start)/2),
                        cbind(b,name="b",summit=(b$end+b$start)/2)),"chr","summit")
afs<-unifyBedFile(pu,
             700
)
afs$chr=as.character(unique(pu$chr))[afs$chr];
afs<-afs[afs$a==afs$b,]
return( data.frame(chr=afs$chr,start=afs$summit-350,end=afs$summit+350,a=afs$a,b=afs$b))
#print(dim(afs)[1])
}

val<-function(a,b){
    un<-dim(bedOverlaps(a,b))[1]
    un#(dim(a)[1]+dim(b)[1]-un)
}

gataTal1<-bedOverlaps(env$bed,gata3)

dim(bedOverlaps(bedNNE,gata3))[1]/(dim(bedNNE)[1]+dim(gata3)[1]-dim(bedOverlaps(bedNNE,gata3))[1])

dim(bedOverlaps(bedNN,gata3))[1]/(dim(bedNN)[1]+dim(gata3)[1])

dim(bedOverlaps(bedCC,gata3))[1]/(dim(bedCC)[1]+dim(gata3)[1])

dim(bedOverlaps(bedCCL,gata3))[1]/(dim(bedCCL)[1]+dim(gata3)[1])

val(bedOverlaps(bedCC,gata3))

dim(bedOverlaps(bedGAL,gata3))

mots<-mapply(function(x,y) union(x,y),locationsC,locationsM)

ou<-addColnames(outer(mots,colnames(env20$reg),Vectorize(function(a,b){
    val(env20$bed[intersect(a, which(env20$reg[,b])),],gata3)
})),colnames(env20$reg))

t(t(ou)/sapply(mots[1:10],length))

ou2<-outer(mots[1:10],apply(env20$reg,2,which),Vectorize(function(x,y) length(intersect(x,y))))


ou2<-outer(gata3,apply(env20$reg,2,which),Vectorize(function(x,y) length(intersect(x,y))))

apply(env20$reg,2,function(x) dim(bedOverlaps(env20$bed[which(x),],gata3)))


write.table(bedCC,"~/Dropbox/UTX-Alex/br-data/bed/bedCC.bed",quote=FALSE,row.names=FALSE)
write.table(bedCCL,"~/Dropbox/UTX-Alex/br-data/bed/bedCCL.bed",quote=FALSE,row.names=FALSE)
write.table(bedNN,"~/Dropbox/UTX-Alex/br-data/bed/bedNN.bed",quote=FALSE,row.names=FALSE)
write.table(bedNNE,"~/Dropbox/UTX-Alex/br-data/bed/bedNNE.bed",quote=FALSE,row.names=FALSE)
write.table(bedGAL,"~/Dropbox/UTX-Alex/br-data/bed/bedGAL.bed",quote=FALSE,row.names=FALSE)

gata3<-addColnames(loadBedFile("~/Dropbox/UTX-Alex/br-data/jurk-gata/single~jurk_gata3_peaks.bed"),c("chr","start","end"))


outer(seq(10),seq(6),Vectorize(function(i,j) ou[i,j]/ou2[i,j]))


bedOverlaps<-function(a,b,l=700){
    pu<-sortDataFrame(rbind(cbind(a,name="a",summit=(a$end+a$start)/2),
                        cbind(b,name="b",summit=(b$end+b$start)/2)),"chr","summit")
    afs<-unifyBedFile(pu,l)
#    afs[afs$a==afs$b,]
    afs
}

ou2<-outer(mots[1:10],apply(env20$reg,2,which),Vectorize(function(x,y) length(intersect(x,y))))


ou2<-outer(gata3,apply(env20$reg,2,function(i) env20$bed[i,]),Vectorize(function(x,y) length(bedOverlaps(x,y,300))))




#outer(mots[1:10],apply(env20$reg,2,function(x) bedOverlaps(env20$bed[which(x),],gata3))





apply(env20$reg,2,function(x) dim(bedOverlaps(env20$bed[which(x),],gata3)))


addFasta2<-function(env,width=300){
    getSeq(BSgenome.Hsapiens.UCSC.hg19,env$chr,start=shiftFromZero(env$summit-150),width=width)
}




a<-addFasta2(bedOverlaps(env$bed[env$reg[,"ECFC"],],gata3))

motifs<-genEboxCombs()
cList<-sapply(motifs,compliment)
locationsM<-lapply(motifs,grep,env20$fasta)
