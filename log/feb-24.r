source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/jan-functions.r")
source("~/r-workspace/feb-functions.r")

source("~/Masters/CCCA/inst/scipts/eboxFrequency.r")

env<-getMotifInfo(addFasta(getPRC20(2)),c(genEboxCombs(),"CANNTG","GATAA"))

                 
mbed<-apply(mergeFun(env$over[,4:dim(env$over)[2]],swapFunB),2,as.logical)

uniqueLoc<-apply(mbed,1,sum)==1

LeukemiaBed<-as.logical(mergeFun(env$over[,4:dim(env$over)[2]],swapFunB)[,"Leukemia"])

env$bed[intersect(env$eboxLoc$CACCTG,which(LeukemiaBed)),]

env$bed[intersect(intersect(env$eboxLoc$CACCTG,which(LeukemiaBed)),which(uniqueLoc)),]

env$bed[intersect(env$eboxLoc$CACCTG,which(env$Leukemia.alt)),]



getMinDistance<-function(fasta,motifa,motifb,loc){
    a<-cbind(which(loc),selectRegionsByDistance(IUPACtoBase(motifa),IUPACtoBase(motifb),loc,env$fasta))
    b<-a[!is.na(a[,2]),]
    b    
}

qstem(a)

qstem<-function(b){
    stem(min(b[,2]):max(b[,2]),getHeights(b[,2]),xlim=c(-32,32))
}

getLoc2Motif<-function(minDistance){
    t<-getHeights(minDistance[,2])
    b<-which.max(t)
    tm<-t
    tm[b]<-0
    c<-which.max(tm)
    c(b+min(minDistance)-1,c+min(minDistance)-1)    
}

selMotifReg<-function(minDistance,...){
    getMin<-function(loc,minDistance)
        minDistance[which(minDistance[,2]==loc),1]
    sort(unique(unlist(lapply(list(...),getMin,minDistance))))
}

findMotifDist<-function(fasta,motifa,motifb,loc,...){
    a<-getMinDistance(fasta,motifa,motifb,loc)
    #getLoc2Motif(a))
    selMotifReg(a,...)
}

gnnbed<-env$bed[findMotifDist(LeukemiaBed,"CANNTG","GATAA",LeukemiaBed,-14,-15,-16),]
gnnpca<-env$bed[findMotifDist(env$Leukemia.alt,"CANNTG","GATAA",env$Leukemia.alt,-14,-15,-16),]
ccbed<-env$bed[intersect(env$eboxLoc$CACCTG,which(LeukemiaBed)),]
ccpca<-env$bed[intersect(env$eboxLoc$CACCTG,which(env$Leukemia.alt)),]


write.table(tts(orderBed(gnnbed)),
            "~/Dropbox/UTX-Alex/br-data/GATA_NN_leukemia_sd_2_bed.txt",
            quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(tts(orderBed(gnnpca)),
            "~/Dropbox/UTX-Alex/br-data/GATA_NN_leukemia_sd_2_pca.txt",
            quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(tts(orderBed(ccbed)),
            "~/Dropbox/UTX-Alex/br-data/CC_leukemia_sd_2_bed.txt",
            quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(tts(orderBed(ccpca)),
            "~/Dropbox/UTX-Alex/br-data/CC_leukemia_sd_2_pca.txt",
            quote=FALSE,col.names=FALSE,row.names=FALSE)


printBB(orderBed(gnnbed),"~/Dropbox/GATA_NN_leukemia_sd_2_bed.bb")
printBB(orderBed(gnnpca),"~/Dropbox/GATA_NN_leukemia_sd_2_pca.bb")
printBB(orderBed(ccbed),"~/Dropbox/CC_leukemia_sd_2_bed.bb")
printBB(orderBed(ccpca),"~/Dropbox/CC_leukemia_sd_2_pca.bb")

a<-getMinDistance(env.fasta,"CANNTG","GATAA",LeukemiaBed)
qstem(a)
minDistance<-a
    t<-getHeights(minDistance[,2])
    b<-which.max(t)
    tm<-t
    tm[b]<-0
    c<-which.max(tm)
tm2<-tm
tm2[c]<-0
d<-which.max(tm2)

c(b+min(minDistance)-1,c+min(minDistance)-1,d+min(minDistance)-1)

t[c(b,c,d)]
