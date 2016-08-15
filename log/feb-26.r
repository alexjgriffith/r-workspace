source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/jan-functions.r")
source("~/r-workspace/feb-functions.r")

source("~/Masters/CCCA/inst/scipts/eboxFrequency.r")


library(seqLogo)

db<-read.table("~/feb-homer/motifAnnotations_20_2-alt.table",header=T)


fm<-sapply(
    as.character(unique(db$filename[ grep("NONE",db$dataset)]))
     , function(x) loadPWM(x))


psl<-function(data){
    seqLogo(makePWM(data),xaxis=FALSE,yaxis=FALSE,ic.scale=FALSE)
}

l=c(6,6,6)
motif<-c("ACCACA","CATCSG","CACCTG")
name=c("Leukemia","Leukemia","Leukemia")

makeLogos<-function(l,motif,name)
    mapply(function(l,motif,name){
        fi<-paste0("~/thesis-feb/",name,"_NONE_pvalue=20_len=",l,"_sd=2_alt.motif")
        loc<-grep(motif,fm[[fi]][,1])
        png(paste0("~/Dropbox/UTX-Alex/br-data/logos/",name,"_",l,"_",motif,".png"))
        psl(as.data.frame(fm[[fi]][loc,3]))
        dev.off()
    }, l,motif,name)


makeLogos(c(7,7,7),c("CGGAARC","TGTTTTC","GATAACA"),c("ECFC","ECFC","ECFC"))

makeLogos(c(6,6,6,7,8),c("CCACAG","GGGGGC","KGCGTG","CTGGCTG","CSCYCCTC"),rep("HSC",5))

makeLogos(c(7,7,8,8),c("CTTATCT","CCAGCTG","GCCTTGTC","GGCTGTCA"),rep("Erythroid",4))


### search for motif with prefered distances
source("~/Masters/CCCA/inst/scipts/eboxFrequency.r")
source("~/r-workspace/preferedDistance.r")

jaspar<-CCCA::loadPWM("~/Masters/mulcal/inst/data/jaspar_motifs.txt",version="jaspar")
jasparMotifs<-cbind(unlist(jaspar[,2]),unlist(lapply(jaspar[,3],function(x) PWMtoCons(x))))


env<-addFasta(getPRC20(2))

ebox<-IUPACtoBase("CANNTG")

mList<-sapply(jasparMotifs[,2],IUPACtoBase)
motifs<-c(mList,ebox)
cList<-sapply(motifs,compliment)
locationsM<-lapply(motifs,grep,env$fasta)
locationsC<-lapply(cList,grep,env$fasta)


## h is backwards need to flip ebox
h<-lapply(1:129,function(i){
    list(num=i,mota=motifs[i],motb=motifs[length(motifs)],distribution=motifHist(env$fasta,motifs,cList,locationsM,locationsC,i,length(motifs),env$reg[,"ALL"]))})
stats<-do.call(rbind,lapply(h,heightStat))
write.table(stats,"~/thesis-feb/jaspar_ebox_stats.txt")

stats[28:dim(stats)[1],]

33,82

n<-6
qstem(h[[n]][[4]],xlim=c(-128,128))
stats[n,]
jasparMotifs[n,]

plotFastMot(env,motifs,cList,locationsM,locationsC,129,length(motifs),xlim=c(-64,64))


plotFastMot(env,motifs,mList,locationsM,locationsC,82,length(motifs),xlim=c(-128,128))


comEboxs<-function(env,motif,xlim=c(-32,32),reg="ALL",q=FALSE){
    motifs2<-c(genEboxCombs(),motif)
    cList2<-sapply(motifs2,compliment)
    lM2<-lapply(motifs2,grep,env$fasta)
    lC2<-lapply(cList2,grep,env$fasta)
    p<-lapply(1:10,function(i) 
    motifHist(env$fasta,motifs2,cList2,lM2,lC2,i,11,env$reg[,reg]))
    if(! q){
    par(mfrow=c(5,2))
    mapply(function(k,n) qstem(k,title=paste0(consenusIUPAC(motifs2[11]),"_",n),xlim),p,as.list(motifs2[1:10]))
}
    p
}


PCSpread<-function(env,motifa,motifb,PC,xlim=c(-32,32)){
    prc<-env$prc$eigenVectors[,PC]
    norm<-normalize(prc)
    motifs2<-c(IUPACtoBase(motifa),IUPACtoBase(motifb))
    cList2<-sapply(motifs2,compliment)
    lM2<-lapply(motifs2,grep,env$fasta)
    lC2<-lapply(cList2,grep,env$fasta)
    x<-seq(-1,1,length=7)
    p<-lapply(seq(6),function(i) {
       reg<-norm>x[i]&norm<x[i+1]
       motifHist(env$fasta,motifs2,cList2,lM2,lC2,1,2,reg)
   })              
    par(mfrow=c(3,2))
    mapply(function(k,n) if(! is.na(k[[1]]))qstem(k,title=paste0(motifa,"_",motifb,"_",n),xlim),p,seq(6))
    p
}

PCSpreadE<-function(env,motifa,motifb,PC,xlim=c(-32,32),q=FALSE){
    prc<-env$prc$eigenVectors[,PC]
    norm<-normalize(prc)
    snorm<-sort(norm)
    step<-floor(length(prc)/10)
    motifs2<-c(IUPACtoBase(motifa),IUPACtoBase(motifb))
    cList2<-sapply(motifs2,compliment)
    lM2<-lapply(motifs2,grep,env$fasta)
    lC2<-lapply(cList2,grep,env$fasta)
    x<-seq(10)
    p<-lapply(x,function(i) {
        mrange<-c(snorm[step*(i-1)+1], snorm[step*i])
        reg<-norm>mrange[1] & norm < mrange[2]       
        motifHist(env$fasta,motifs2,cList2,lM2,lC2,1,2,reg)
   })
    if(! q){
    par(mfrow=c(5,2))
    mapply(function(k,i) if(! is.na(k[[1]]))qstem(k,title=paste0(motifa,"_",motifb,"_",round(snorm[(step*(i-1)+1)],2) ,":",round(snorm[(step*i)],2)),xlim),p,x)
}
    p
}


#png("~/Dropbox/UTX-Alex/br-data/ebox/PC2_CANNTG_H1F1A_ARNT.png")
p<-PCSpreadE(env,consenusIUPAC(motifs[129]),"CANNTG",2,xlim=c(-62,62))
                                        #dev.off()

imBody<-function(p,mrange){
    par(mar=c(0,0,0,0))
    image(t(as.matrix(do.call(rbind,lapply(lapply(p,qstem,xlim=mrange,q=TRUE),getHeights,mrange)))), xaxt='n', yaxt='n', ann=FALSE)
}

makeImage<-function(env,motifa,motifb,pc,mrange=c(-32,32)){
    p<-PCSpreadE(env,motifa,motifb,pc,xlim=mrange,q=TRUE)
    imBody(p,mrange)
}

png("~/Dropbox/UTX-Alex/br-data/ebox/PC1_CANNTG-GATAA.png")
makeImage(env,"CANNTG","GATAA",1)
dev.off()
png("~/Dropbox/UTX-Alex/br-data/ebox/PC2_CANNTG-GATAA.png")
makeImage(env,"CANNTG","GATAA",2)
dev.off()
png("~/Dropbox/UTX-Alex/br-data/ebox/PC4_CANNTG-GATAA.png")
makeImage(env,"CANNTG","GATAA",4)
dev.off()
png("~/Dropbox/UTX-Alex/br-data/ebox/PC1_CANNTG-ARNT.png")
makeImage(env,"CANNTG",consenusIUPAC(motifs[129]),1,c(-62,62))
dev.off()
png("~/Dropbox/UTX-Alex/br-data/ebox/PC2_CANNTG-ARNT.png")
makeImage(env,"CANNTG",consenusIUPAC(motifs[129]),2,c(-62,62))
dev.off()
png("~/Dropbox/UTX-Alex/br-data/ebox/PC4_CANNTG-ARNT.png")
makeImage(env,"CANNTG",consenusIUPAC(motifs[129]),4,c(-62,62))
dev.off()
png("~/Dropbox/UTX-Alex/br-data/ebox/PC1_CANNTG-TEB.png")
makeImage(env,"CANNTG",consenusIUPAC(motifs[69]),1,c(-62,62))
dev.off()
png("~/Dropbox/UTX-Alex/br-data/ebox/PC2_CANNTG-TEB.png")
makeImage(env,"CANNTG",consenusIUPAC(motifs[69]),2,c(-62,62))
dev.off()
png("~/Dropbox/UTX-Alex/br-data/ebox/PC4_CANNTG-TEB.png")
makeImage(env,"CANNTG",consenusIUPAC(motifs[69]),4,c(-62,62))
dev.off()
png("~/Dropbox/UTX-Alex/br-data/ebox/PC1_CANNTG-PAX5.png")
makeImage(env,"CANNTG",consenusIUPAC(motifs[6]),1,c(-62,62))
dev.off()
png("~/Dropbox/UTX-Alex/br-data/ebox/PC2_CANNTG-PAX5.png")
makeImage(env,"CANNTG",consenusIUPAC(motifs[6]),2,c(-62,62))
dev.off()
png("~/Dropbox/UTX-Alex/br-data/ebox/PC4_CANNTG-PAX5.png")
makeImage(env,"CANNTG",consenusIUPAC(motifs[6]),4,c(-62,62))
dev.off()
png("~/Dropbox/UTX-Alex/br-data/ebox/PC1_CANNTG-AGGCCG.png")
makeImage(env,"CANNTG","AGGCCG",1,c(-32,32))
dev.off()
png("~/Dropbox/UTX-Alex/br-data/ebox/PC2_CANNTG-AGGCCG.png")
makeImage(env,"CANNTG","AGGCCG",2,c(-32,32))
dev.off()
png("~/Dropbox/UTX-Alex/br-data/ebox/PC4_CANNTG-AGGCCG.png")
makeImage(env,"CANNTG","AGGCCG",4,c(-32,32))
dev.off()



png("~/Dropbox/UTX-Alex/br-data/ebox/PC1_CACCTG-TEB.png")
makeImage(env,"CACCTG",consenusIUPAC(motifs[69]),1,c(-62,62))
dev.off()
png("~/Dropbox/UTX-Alex/br-data/ebox/PC2_CACCTG-TEB.png")
makeImage(env,"CACCTG",consenusIUPAC(motifs[69]),2,c(-62,62))
dev.off()
png("~/Dropbox/UTX-Alex/br-data/ebox/PC4_CACCTG-TEB.png")
makeImage(env,"CACCTG",consenusIUPAC(motifs[69]),4,c(-62,62))
dev.off()


png("~/Dropbox/UTX-Alex/br-data/ebox/PC1_CAGCTG-TEB.png")
makeImage(env,"CAGCTG",consenusIUPAC(motifs[69]),1,c(-62,62))
dev.off()
png("~/Dropbox/UTX-Alex/br-data/ebox/PC2_CAGCTG-TEB.png")
makeImage(env,"CAGCTG",consenusIUPAC(motifs[69]),2,c(-62,62))
dev.off()
png("~/Dropbox/UTX-Alex/br-data/ebox/PC4_CAGCTG-TEB.png")
makeImage(env,"CAGCTG",consenusIUPAC(motifs[69]),4,c(-62,62))
dev.off()


png("~/Dropbox/UTX-Alex/br-data/ebox/CANNTG_H1F1A_ARNT.png",height=240,width=240)
qstem(h[[n]][[4]]*-1,xlim=c(-64,64))
dev.off()




PCSpread2<-function(env,motifa,motifb,PC,xlim=c(-32,32)){
    prc<-env$prc$eigenVectors[,PC]
    norm<-normalize(prc)
    motifs2<-c(IUPACtoBase(motifa),IUPACtoBase(motifb))
    cList2<-sapply(motifs2,compliment)
    lM2<-lapply(motifs2,grep,env$fasta)
    lC2<-lapply(cList2,grep,env$fasta)
    p<-list(motifHist(env$fasta,motifs2,cList2,lM2,lC2,1,2,prc<(-0.5)),
         motifHist(env$fasta,motifs2,cList2,lM2,lC2,1,2,prc>(0.5))
        )
    par(mfrow=c(1,2))
    mapply(function(k,n) if(! is.na(k))qstem(k,title=paste0(motifa,"_",motifb,"_",n),xlim),p,seq(2))
    p
}


qstem(hist2Motifs(env,"CANNTG",consenusIUPAC( motifs[129]),"ALL"),xlim=c(-62,62))




plotFastMot(env,motifs,cList,locationsM,locationsC,129,length(motifs),xlim=c(-64,64))

### Fix the motifHist function

motifs2<-c(genEboxCombs(),"GATAA",ebox,"AGGCCG")
cList2<-sapply(motifs2,compliment)
lM2<-lapply(motifs2,grep,env$fasta)
lC2<-lapply(cList2,grep,env$fasta)


par(mfrow=c(1,2))
qstem(cbind(0,motifHist(env$fasta,motifs2,cList2,lM2,lC2,12,11,env$reg[,"Erythroid"]))
     ,"ALL")
qstem(cbind(0,motifHist(env$fasta,motifs2,cList2,lM2,lC2,12,11,env$reg[,"Erythroid"]))
     ,"ALL")

plotFastMot(env,motifs2,cList2,lM2,lC2,12,13)

mplotFastMot(env,motifs2,cList2,lM2,lC2,13,1:10,"Leukemia")

reg<-"ALL"
p<-lapply(1:10,function(i) 
    motifHist(env$fasta,motifs2,cList2,lM2,lC2,i,13,env$reg[,reg]))

#png("~/Dropbox/UTX-Alex/br-data/ebox/Eboxs-AGGCCG.png",width=240,height=240)
par(mfrow=c(5,2))
mapply(function(k,n) qstem(k,title=paste0(motifs2[13],"_",n),xlim=c(-32,32)),p,as.list(motifs2[1:10]))
#dev.off()

#png("~/Dropbox/UTX-Alex/br-data/ebox/Leukemic_CANNTG-AGGCCG.png",width=240)
par(mfrow=c(2,1))
qstem(motifHist(env$fasta,motifs2,cList2,lM2,lC2,6,13,env$reg[,"Leukemia"]),"Leukemic")
qstem(motifHist(env$fasta,motifs2,cList2,lM2,lC2,6,13,!env$reg[,"Leukemia"]),"Non Leukemic")
#dev.off()

png("~/Dropbox/UTX-Alex/br-data/ebox/CACCTG-TBP.png",height=240,width=240)
qstem(hist2Motifs(env,"CACCTG",consenusIUPAC(motifs[69]),"ALL"),xlim=c(-64,64))
dev.off()


png("~/Dropbox/UTX-Alex/br-data/ebox/CAGCTG-TBP.png",height=240,width=240)
p<-qstem(hist2Motifs(env,"CAGCTG",consenusIUPAC(motifs[69]),"ALL"),xlim=c(-64,64))
dev.off()


png("~/Dropbox/UTX-Alex/br-data/ebox/CANNTG-MIZF.png",height=240,width=240)
qstem(hist2Motifs(env,"CANNTG",consenusIUPAC(motifs[82]),"ALL"),xlim=c(-32,32))
dev.off()




png("~/Dropbox/UTX-Alex/br-data/ebox/CANNTG-MIZF.png",height=240,width=240)
qstem(hist2Motifs(env,"CANNTG",consenusIUPAC(motifs[82]),"ALL"),xlim=c(-32,32))
dev.off()


munit<-function(x) (x+min(x))/(max(x)-min(x))

eboxImage<-function(env,motif,mrange=c(-32,32),reg="ALL",norm=pass){
    imBody<-function(p,mrange,norm=pass){
        par(mar=c(0,0,0,0))
        makeData<-function(x)
            norm(getHeights(qstem(x,xlim=mrange,q=TRUE),mrange))
        data<-lapply(p, makeData)        
        image(t(as.matrix(do.call(rbind,data))), xaxt='n', yaxt='n', ann=FALSE)
    }    
    p<-comEboxs(env,motif,xlim,reg=reg,q=TRUE)
    imBody(p,mrange,norm)
}

generateHML<-function(){
png("~/Dropbox/UTX-Alex/br-data/ebox/HM_leukemia_CANNTG-GATA.png")
eboxImage(env,"GATAA",norm=function(x){log(x+1,10)},reg="Leukemia")
dev.off()
png("~/Dropbox/UTX-Alex/br-data/ebox/HM_erythroid_CANNTG-GATA.png")
eboxImage(env,"GATAA",norm=function(x){log(x+1,10)},reg="Erythroid")
dev.off()
png("~/Dropbox/UTX-Alex/br-data/ebox/HM_ecfc_CANNTG-GATA.png")
eboxImage(env,"GATAA",norm=function(x){log(x+1,10)},reg="ECFC")
dev.off()
png("~/Dropbox/UTX-Alex/br-data/ebox/HM_hsc_CANNTG-GATA.png")
eboxImage(env,"GATAA",norm=function(x){log(x+1,10)},reg="HSC")
dev.off()
}


eboxImage(env,motifs[69],norm=function(x){log(x+1,10)},reg="ALL",mrange=c(-62,62))


eboxImage(env,motifs[69],norm=pass,mrange=c(-62,62))


png("~/Dropbox/UTX-Alex/br-data/ebox/eboxs-tbp.png")
comEboxs(env,motifs[69],c(-62,62))
dev.off()



comEboxs(env,"TAGTTA",c(-32,32))

png("~/Dropbox/UTX-Alex/br-data/ebox/AGGCGG-SCACTG.png",height=240,width=240)
qstem(hist2Motifs(env,"AGGCGG","SCACTG","ALL"),xlim=c(-64,64))
dev.off()

png("~/Dropbox/UTX-Alex/br-data/ebox/PC1_AGGCGG-SCACTG.png")
makeImage(env,"AGGCGG","SCACTG",1,c(-64,64))
dev.off()
png("~/Dropbox/UTX-Alex/br-data/ebox/PC2_AGGCGG-SCACTG.png")
makeImage(env,"AGGCGG","SCACTG",2,c(-64,64))
dev.off()
png("~/Dropbox/UTX-Alex/br-data/ebox/PC4_AGGCGG-SCACTG.png")
makeImage(env,"AGGCGG","SCACTG",4,c(-64,64))
dev.off()


fastMot2(env,"AGGCGG","SCACTG")
