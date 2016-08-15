source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/jan-functions.r")
source("~/r-workspace/feb-functions.r")
source("~/r-workspace/feb-variables.r")

source("~/Dropbox/thesis/R/contrib_new.R")



qhwF<-function(fasta,r,p1,p2,pvalue,len=8,sd=0){
    print(paste("~/thesis-feb/",p1,"_",p2,"_pvalue=",pvalue,"_len=",len,"_sd=",sd,".motif",sep=""))
    homerWrapper(fasta,r[,p1],r[,p2],"~/Masters/mulcal/inst/lib/homer-4.7/bin/homer2",paste("~/thesis-feb/",p1,"_",p2,"_pvalue=",pvalue,"_len=",len,"_sd=",sd,".motif",sep=""),opts=paste("-S 25 -len ",len," > /tmp/homerTrash 2>&1 ",sep=""))
} 

env<-getPRC20(2)
env$fasta<-getSeq(BSgenome.Hsapiens.UCSC.hg19,env$bed$chr,start=env$bed$start+150,width=300)

    
matrix<-with(env,{
    sapply(list(Erythroid.alt,ECFC.alt,Leukemia.alt,HSC.alt),function(x)length(which(x)))
})

callHomer<-function(env,pvalue,sd){
matrix<-with(env,{
    r<-cbind(ECFC=ECFC.alt,NONE=NONE,ALL=ALL,
             Erythroid=Erythroid.alt,Leukemia=Leukemia.alt,HSC=HSC.alt)
    for(b in c("NONE","ALL")){
        for(a in c("ECFC","HSC","Erythroid","Leukemia")){
            for(len in c(6,7,8))
                qhwF(fasta,r             
                    ,a,b,pvalue,len,sd)
        }
    }
})
}

callHomerALL<-function(env,pvalue,sd){
matrix<-with(env,{
    r<-cbind(ECFC=ECFC,NONE=NONE,ALL=ALL,
             Erythroid=Erythroid,Leukemia=Leukemia,HSC=HSC)
    for(b in c("NONE","ALL")){
        for(a in c("ECFC","HSC","Erythroid","Leukemia")){
            for(len in c(6,7,8))
                qhwF(fasta,r             
                    ,a,b,pvalue,len,paste0(sd,"_alt"))
        }
    }
})
}


stampDF<-function(combs,pvalue,sd,dir="~/thesis-feb/")
    do.call(rbind,apply(combs,1,function(x,pvalue) data.frame(file=paste(dir,x[1],"_pvalue=",pvalue,"_len=",x[2],"_sd=",sd,".motif",sep=""),compatator=x[3],name=x[1],size=x[2],pvalue=pvalue),pvalue ))


sd=3
combs20<-stampDF(cbind(strsplit("Erythroid_NONE Erythroid_ALL Leukemia_NONE Leukemia_ALL HSC_NONE HSC_ALL ECFC_NONE ECFC_ALL"," ")[[1]],c(rep(6,8),rep(7,8),rep(8,8)),c(rep("TRANSFAC_Fams",24),rep("JASPAR_Fams",24))),pvalue=20,sd=sd)

combs5<-stampDF(cbind(strsplit("Erythroid_NONE Erythroid_ALL Leukemia_NONE Leukemia_ALL HSC_NONE HSC_ALL ECFC_NONE ECFC_ALL"," ")[[1]],c(rep(6,8),rep(7,8),rep(8,8)),c(rep("TRANSFAC_Fams",24),rep("JASPAR_Fams",24))),pvalue=5,sd=sd,"~/thesis-feb/")

allCombs<-rbind(combs5,combs20)
df<-allCombs

i=1
x<-stampWrapper(as.character(df$file[i]),df[i,2],df[i,3])

db<-multiStamp(allCombs)

write.table(db,"~/thesis-feb/motifAnnotations_sd_3.table",quote=FALSE,row.names=FALSE)


ahome<-function(){
    env<-getPRC5(3);
    env$fasta<-getSeq(BSgenome.Hsapiens.UCSC.hg19,env$bed$chr,start=env$bed$start+150,width=300)
    callHomer(env,5,3)
    callHomerALL(env,5,3)
    rm(env)
    env<-getPRC20(3);
    env$fasta<-getSeq(BSgenome.Hsapiens.UCSC.hg19,env$bed$chr,start=env$bed$start+150,width=300)
    callHomer(env,20,3)
    callHomerALL(env,20,3)
    rm(env)
    env<-getPRC5(1);
    env$fasta<-getSeq(BSgenome.Hsapiens.UCSC.hg19,env$bed$chr,start=env$bed$start+150,width=300)
    callHomerALL(env,5,1)
    rm(env)
    env<-getPRC20(1);
    env$fasta<-getSeq(BSgenome.Hsapiens.UCSC.hg19,env$bed$chr,start=env$bed$start+150,width=300)
    callHomerALL(env,20,1)
    rm(env)

                                        #env<-getPRC5(2);
#env$fasta<-getSeq(BSgenome.Hsapiens.UCSC.hg19,env$bed$chr,start=env$bed$start+150,width=300)
#callHomer(env,5,2)
#env<-getPRC5(1);
#env$fasta<-getSeq(BSgenome.Hsapiens.UCSC.hg19,env$bed$chr,start=env$bed$start+150,width=300)
#callHomer(env,5,1)
#env<-getPRC20(2);
#env$fasta<-getSeq(BSgenome.Hsapiens.UCSC.hg19,env$bed$chr,start=env$bed$start+150,width=300)
#callHomer(env,20,2)
#rm(env)
#env<-getPRC20(1);
#env$fasta<-getSeq(BSgenome.Hsapiens.UCSC.hg19,env$bed$chr,start=env$bed$start+150,width=300)
#callHomer(env,20,1)
}

ahome<-function(){
#env<-getPRC5(2);
#env$fasta<-getSeq(BSgenome.Hsapiens.UCSC.hg19,env$bed$chr,start=env$bed$start+150,width=300)
#callHomer(env,5,2)
#env<-getPRC5(1);
#env$fasta<-getSeq(BSgenome.Hsapiens.UCSC.hg19,env$bed$chr,start=env$bed$start+150,width=300)
#callHomer(env,5,1)
env<-getPRC20(2);
env$fasta<-getSeq(BSgenome.Hsapiens.UCSC.hg19,env$bed$chr,start=env$bed$start+150,width=300)
callHomer(env,20,2)
rm(env)
env<-getPRC20(1);
env$fasta<-getSeq(BSgenome.Hsapiens.UCSC.hg19,env$bed$chr,start=env$bed$start+150,width=300)
callHomer(env,20,1)
}

with(env,{
    n=1
    par(mfrow=c(2,2))
    x<-apply(cbind(qn(heights)[Erythroid.alt,swapFunD(categories)=="Erythroid"],0),1,max)
    hist(x,xlab="",ylab="",main="Erythroid",xlim=c(0,200),breaks=500)
    x<-apply(cbind(qn(heights)[Leukemia.alt,swapFunD(categories)=="Leukemia"],0),1,max)
    hist(x,xlab="",ylab="",main="Leukemia",xlim=c(0,200),breaks=500)
    x<-apply(cbind(qn(heights)[HSC.alt,swapFunD(categories)=="HSC"],0),1,max)
    hist(x,xlab="",ylab="",main="HSC",xlim=c(0,200),breaks=500)
    x<-apply(cbind(qn(heights)[ECFC.alt,swapFunD(categories)=="ECFC"],0),1,max)
    hist(x,xlab="",ylab="",main="ECFC",xlim=c(0,200),breaks=500)     
})


env5<-getPRC(5,"combined")

with(env5,{
    #pc1=1
    #pc2=3
    #plotPCMat2D(pca2Matr(prc),c(pc1,pc2),categories,swapFun,swapFunD,swapFunC)
    plotPCMat2D(pca2Matr(prc),c(3,5),categories,swapFun,swapFunD,swapFunC)
})
