# combine multipe peak files based on their similarity and shared annotations
# run the pariwisealignment on the cluster for significant time savings

source("~/r-workspace/CCCA.r")
source("~/r-workspace/project.r")

stampDF<-function(combs,type,pvalue=20,sd=2,dir="~/Dropbox/Data/homer-paper/homer_")
    do.call(rbind,apply(combs,1,function(x) data.frame(file=paste(dir,x[1],"_",x[2],".txt",sep=""),compatator=x[3],name=x[1],size=x[2],pvalue=pvalue,sd=sd)))


combs<-stampDF(cbind(strsplit("Erythroid_PCA Erythroid_BED Leukemia_PCA Leukemia_BED HSC_PCA HSC_BED ECFC_PCA ECFC_BED"," ")[[1]],c(rep(6,8),rep(7,8),rep(8,8)),c(rep("TRANSFAC_Fams",24),rep("JASPAR_Fams",24))))


pwms<-apply(combs,1,function(x) loadPWM(x[1]))

names<-as.character(combs$name)
names(pwms)<-names

library(Biostrings)



Biostrings::pairwiseAlignment(motifs[1],motifs[1])

cs<-makeForkCluster(100)

names<-unique(as.character(combs$name))
as<-lapply(names,function(name){
    motifs<-gsub(">","",unlist(lapply(pwms[combs$name==name],function(x) x[,1])))
com<-t(combn(motifs,2))
a<-parApply(cs,com,1,function(x) Biostrings::pairwiseAlignment(x[1],x[2],scoreOnly=TRUE))
})
names(as,names)

stopCluster(cs)

Xs<-lapply(as,function(a) matrix(0,ncol=length(motifs),nrow=length(motifs)))

for(X in Xs)
    X[lower.tri(X)]<-a


ds<-lapply(Xs,function(X) as.dist(X))


save(as,Xs,ds,file="~/Dropbox/clust.RData")

svg("~/Dropbox/clust.svg")
plot(clust,hang=-1)
dev.off()


load("~/Dropbox/clust.RData")

png("~/Dropbox/clust.png",width=1020)
plot(hclust(as.dist(d[2]),method="ward.D2"),hang=-1)
dev.off()


db<-read.table("~/Dropbox/Data/homer-paper/motifAnnotations.txt",header=T)


meq<-function(a,b){
    Reduce("|",lapply(a,"==",b))
}

names(Xs)<-unique(names)
    
tfun<-function(set){
    clust<-hclust(as.dist(Xs[[set]]),method="ward.D2")
    regs<-cutree(clust,k=30)
    motifs<-colnames(Xs[[set]])
    sel<-lapply(do.call(seq,as.list(range(regs))),function(r){
        tm<-motifs[regs==r]
        subt<-db[meq(tm,db$motif)  &db$dataset==set, ]
        ms<-lapply(tm,function(x)as.character(subt$genes[subt$motif==x]))
        cons<-do.call(intersectN,ms)
        if(length(cons)==0){        
            cbind(subt[order(subt$escore)[1:3],c("tmotif","genes","escore","pvalue")],amotif=paste(tm,collapse=","))
        }
        else{
            mm<-as.character(subt$tmotif[which.max(subt$pvalue)])
            ftab<-cbind(subt[meq(cons,subt$genes)&subt$tmotif==mm,c("tmotif","genes","escore","pvalue")],amotif=paste(tm,collapse=","))
            pv<-subt[subt$tmotif==mm,"pvalue"]
            ftab<-rbind(ftab,data.frame(tmotif=rep(mm,3),genes="NA",escore=1,pvalue=rep(pv,3),amotif=paste(tm,collapse=",")))
            ftab[order(ftab$escore)[1:3],]
        }
    })

    subt<-do.call(rbind,sel)

    osubt<-subt[order(subt$pvalue,decreasing=TRUE),]

    mots<-unique(as.character(unlist(subt$tmotif)))

    findPWM<-function(dataset,motif,pwms){
        a<-do.call(rbind,lapply(which(dataset==names),function(x) pwms[[x]]))
        f<-paste0(">",motif)
        print(f)
        a[a[,1]==f,3][[1]]
    }

    flatenPWM<-function(pwm)
        paste0("[",paste(c("A","C","G","T"),
                         apply(a[a[,1]==f,3][[1]],1,
                               function(x)
                                   paste(pwm,collapse=",")),collapse=";",sep=","),"]")
    print(set)
    print(str(as.character(osubt$tmotif)))
    out<-cbind(osubt,as.character(unlist(sapply(as.character(osubt$tmotif),function(x) flatenPWM(findPWM(set,x,pwms))))))
    
    do.call(rbind,lapply( seq(dim(out)[1]/3),function(i){
        print(i)
        p=((i-1)*3)+1
        data.frame(motif=out[p,1],motifs=out[p,5],pvalue=out[p,4],
                   annot1=out[p,2],evalue1=out[p,3],
                   annot2=out[p+1,2],evalue2=out[p+1,3],
                   annot3=out[p+2,2],evalue3=out[p+2,3],
                   pwm=out[p,6])
    }))
}

res<-lapply(as.character(unique(names)),tfun)

names(res)<-as.character(unique(names))
    
    
categories<-as.character(unlist(read.table("~/Dropbox/Data/categories/22-categories.txt")))

contexts<-unique(swapFunD(categories))

file<-path.expand("~/Dropbox/Data/test.xlsx")

library(xlsx)

wb<-createWorkbook()

#cols<-c("ann","motif","dataset","order","rank","escore","MACS","pvalue")

for(treat in contexts)
    for(control in c("PCA","BED")){
        sn<-paste0(treat,"_",control)
        ts<-createSheet(wb,sn)
        addDataFrame(res[[sn]],sheet=ts,row.names = FALSE)
        print(sn)
    }


saveWorkbook(wb,file)

