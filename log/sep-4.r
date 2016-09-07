library(CCCA)

library('BSgenome.Hsapiens.UCSC.hg19')

source("~/r-workspace/project/project.r")
source("~/r-workspace/project/project-variables.r")
source("~/r-workspace/project/ccca.r")

env<-getPRC20(2) ;; redone with ECFC

apply(env$reg,2,sum)

env<-addFasta(env)
## Find motifs ECFC

#CCCA::homerWrapper(env$fasta,env$reg[,"ECFC"],)

## Call Homer
i="ECFC"
homerWrapper(env$fasta,env$reg[,i],env$reg[,"NONE"],"~/binaries/homer/bin/homer2",motifsFile=paste0("~/Dropbox/Data/homer_",i,"_sqrt_PCA_6.txt"),opts="-S 25 -len 6")
homerWrapper(env$fasta,env$reg[,i],env$reg[,"NONE"],"~/binaries/homer/bin/homer2",motifsFile=paste0("~/Dropbox/Data/homer_",i,"_sqrt_PCA_7.txt"),opts="-S 25 -len 7")
homerWrapper(env$fasta,env$reg[,i],env$reg[,"NONE"],"~/binaries/homer/bin/homer2",motifsFile=paste0("~/Dropbox/Data/homer_",i,"_sqrt_PCA_8.txt"),opts="-S 25 -len 8")

## Rebuild stamp db
stampDF<-function(combs,type,pvalue=20,sd=2,dir="~/Dropbox/Data/homer-paper/homer_")    
    do.call(rbind,apply(combs,1,function(x) data.frame(file=paste(dir,x[1],"_",x[2],".txt",sep=""),compatator=x[3],name=x[1],size=x[2],pvalue=pvalue,sd=sd)))

combs<-stampDF(cbind(strsplit("Erythroid_PCA Erythroid_BED Leukemia_PCA Leukemia_BED HSC_PCA HSC_BED ECFC_sqrt_PCA ECFC_BED"," ")[[1]],c(rep(6,8),rep(7,8),rep(8,8)),c(rep("TRANSFAC_Fams",24),rep("JASPAR_Fams",24))))

library(httr)

db<-multiStamp(combs)

write.table(db,"~/Dropbox/Data/homer-paper/motifAnnotations_ecfc_sqrt.txt",quote=FALSE,row.names=FALSE)

## cluster motifs

db<-read.table("~/Dropbox/Data/homer-paper/motifAnnotations_ecfc_sqrt.txt",header=T)

buildFM("ECFC_sqrt","PCA",6,20)


pwms<-apply(combs,1,function(x) loadPWM(x[1]))

names<-as.character(combs$name)
names(pwms)<-names

library(Biostrings)

mot<-function(n,j,l=8)sapply(seq(n)-1,function(i) l*i+j)

motifs<-gsub(">","",unlist(lapply(pwms[mot(4,7)],function(x) x[,1])))

ecfc<-apply(combn(motifs,2),2,function(x) Biostrings::pairwiseAlignment(x[1],x[2],scoreOnly=TRUE))

## relies heavily on functions from aug-1.r
       
names<-unique(as.character(combs$name))

as<-lapply(names,function(name){
    motifs<-unique(gsub(">","",unlist(lapply(pwms[combs$name==name],function(x) x[,1]))))
    com<-t(combn(motifs,2))
    a<-apply(com,1,function(x) Biostrings::pairwiseAlignment(x[1],x[2],scoreOnly=TRUE))
})


motifs<-lapply(names,function(name){
    motifs<-unique(gsub(">","",unlist(lapply(pwms[combs$name==name],function(x) x[,1]))))
})


names(as)<-names

stopCluster(cs)

Xs<-lapply(motifs,function(motifs) matrix(0,ncol=length(motifs),nrow=length(motifs)))

for(i in seq(length(Xs))){
    Xs[[i]][lower.tri(Xs[[i]])]<-as[[i]]
    colnames(Xs[[i]])<-motifs[[i]]
    rownames(Xs[[i]])<-motifs[[i]]
    t
}

ds<-lapply(Xs,function(X) as.dist(X))


save(as,Xs,ds,file="~/Dropbox/clust_sqrt.RData")


clust<-hclust(ds[[2]],method="ward.D2")

regs<-cutree(clust,k=30)

subdb<-db[meq(motifs[[2]][regs==1],db$motif) & db$dataset==names[[2]],]


#lapply(pwms,function(x) x[,2])

subdbOrder<-subdb[order(subdb$pvalue,subdb$escore),]

mergeBase()


mergeBase<-function(motifs){
    sapply(motifs,CCCA::IUPACtoBase)
    Reduce("==",sapply(motifs,nchar))
}

mots<-motifs[[2]][regs==2]

                                        #CCCA::consenusIUPAC()

devtools::install_github("Bioconductor-mirror/msa")

#library(msa)

msa:
CCCA::IUPACtoBase("CANNTG",TRUE)

