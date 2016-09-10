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

subdb<-db[meq(motifs[[1]][regs==1],db$motif) & db$dataset==names[[2]],]


#lapply(pwms,function(x) x[,2])

subdbOrder<-subdb[order(subdb$pvalue,subdb$escore),]

mergeBase()


mergeBase<-function(motifs){
    sapply(motifs,CCCA::IUPACtoBase)
    Reduce("==",sapply(motifs,nchar))
}

mots<-motifs[[2]][regs==2]





combineClusterdMotifs<-function(name,names,ds,pwms){    
    require(msa)
    clust<-hclust(ds[[name]],method="ward.D2")
    regs<-cutree(clust,k=30)
    ns<-unique(names)
    ts<-do.call(rbind,pwms[names==name])
    sapply(do.call(seq,as.list(range(regs))),
           function(i){
               tn<-names(which(regs==i))
               if(length(tn)>1){
                   log<-capture.output({samp<-consensusMatrix(msa(tn,type="dna"))})
                   mergeConsensus(samp)
                   
               }
               else
                   tn
               })
}

mergeBase<-function(string){
    temp<-paste0(sort(unique(Filter(function(x) !x %in% c("[","]",".","+","-"),strsplit(string,"")[[1]]))),collapse="")
    if(nchar(temp)==1)
        temp
    else
        paste0("[",temp,"]")
}

mergeConsensus<-function(samp){
    require(msa)
    sampNames<-rownames(samp)
    consenusIUPAC(paste(sapply(seq(dim(samp)[2]),function(i) mergeBase(IUPACtoBase(paste0(Filter(function(x) !x %in% c("[","]",".","+","-"),sampNames[samp[,i]>0]),collapse="")))),collapse=""))
}



consensus2pwm<-function(string){
    chars<-strsplit(string,"")[[1]]
    n<-length(chars)
    do.call(rbind,lapply(chars,function(x){
        a<-matrix(rep(0.001,4),ncol=4)
        colnames(a)<-c("A","C","T","G")
        i<-Filter(function(i) !i %in% c("]","["), strsplit(IUPACtoBase(x),"")[[1]])
        a[,i]<-1/length(i)-(4-length(i))*0.001
        a
    }))                      
}

addMotifHeader<-function(motif,pwm){
    paste0(">",motif,"\t","0_",motif,"\t",100,"\t",-20,"\n",paste(apply(pwm,1,paste,collapse=" "),collapse="\n"),"\n")
}

motifList2Homer<-function(motif){
    paste(sapply(motif,function(m)
        addMotifHeader(m,consensus2pwm(m))),collapse="")
}

scorePWMHomer<-function(){
    spl<-strsplit(motif,"")[[1]]
    sum(sapply(spl,function(x) log(0.997/(0.25*length(IUPACtoBase(x,rl=TRUE))))))
    
}


scoreMotifHomer<-function(motif){
    spl<-strsplit(motif,"")[[1]]
    sum(sapply(spl,function(x) log(0.997/(0.25*length(IUPACtoBase(x,rl=TRUE))))))
    
}

lapply(seq(3),function(n) log(apply(combn(tmp[,3][[1]][,1],n),2,sum)/(0.250)))

log(.997/0.25)*5+log(.001/0.25)
mots<-combineClusterdMotifs("Leukemia_PCA",combs$name,ds,pwms)


consenusIUPAC(motifString(tmp[,3][[1]]))

cat(motifList2Homer(mots),file="~/Desktop/motifs.motifs")

cat(motifList2Homer("BBB"),file="~/Desktop/motifs.motifs")


homerWrapperKnown(ccca$fasta,ccca$reg[,"Leukemia"],ccca$reg[,"NONE"],"~/binaries/homer/bin/homer2",mots)

lapply(mots,function(mot){
sequences<-c(ccca$fasta[ccca$reg[,"Leukemia"]],ccca$fasta[ccca$reg[,"NONE"]])
n<-length(ccca$fasta[ccca$reg[,"Leukemia"]])
a<-grep(paste0(IUPACtoBase(compliment(mot)),"\\\\|",IUPACtoBase(compliment(mot))),sequences)
yes<-c(sum(a<=n),sum(a>n))
tot<-c(n,length(sequences)-n)
no<-tot-yes
print(yes/no)
test<-chisq.test(rbind(yes,no))
test$p.value
})

fg<-ccca$fasta[ccca$reg[,"Leukemia"]]
bg<-ccca$fasta[ccca$reg[,"NONE"]]

baseMots<-
revMots<-compliment(IUPACtoBase(mots[1]))

count<-function(list,ops=sort(unique(list))){    
    ret<-sapply(ops,function(x) sum(list==x))
    names(ret)<-ops
    ret
}

countMotifs<-function(sequence,motifs){
    ops<-IUPACtoBase(motifs)
    oprs<-IUPACtoBase(motifs,TRUE)
    oprsrev<-sapply(oprs,compliment)
    ret<-rbind(count(unlist(lapply(Filter(function(x) length(x)>0, str_match_all(sequence,ops)),unique)),oprs),
               count(unlist(lapply(Filter(function(x) length(x)>0, str_match_all(sequence,compliment(ops))),unique)),oprsrev))
    rownames(ret)<-c("motif","complement")
    ret
}

countMotifsOne<-function(sequence,motifs){
    oprsfor<-IUPACtoBase(motifs,TRUE)
    oprsrev<-sapply(oprsfor,compliment)
    ops<-paste(oprsrev,oprsfor,collapse="|",sep="|")
    sum(grepl(ops,sequence))
    
}


countMotifsWG<-function(sequence,motifs,wiggle){
    oprs<-wiggle(motifs)
    ops<-paste0(oprs,collapse="|")
    oprsrev<-sapply(oprs,compliment)
    opsrev<-paste0(oprsrev,collapse="|")
    ret<-rbind(count(unlist(lapply(Filter(function(x) length(x)>0, str_match_all(sequence,ops)),unique)),oprs),
               count(unlist(lapply(Filter(function(x) length(x)>0, str_match_all(sequence,opsrev)),unique)),oprsrev))
    rownames(ret)<-c("motif","complement")
    ret
}

count2freq(counts2PWM(countMotifsWG(fg,"GATAA",wiggle)))

counts2PWM(countMotifs(fg,"GATAA"))

wiggle<-function(motif,n=1){
    spl<-strsplit(motif,"")[[1]]
    unique(unlist(sapply(seq(length(spl)),function(i) {ret<-spl; ret[i]<-"N"; IUPACtoBase(paste0(ret,collapse=""),TRUE)})))
    
}

consensus2pwm<-function(string){
    chars<-strsplit(string,"")[[1]]
    n<-length(chars)
    do.call(rbind,lapply(chars,function(x){
        a<-matrix(rep(0.001,4),ncol=4)
        colnames(a)<-c("A","C","T","G")
        i<-Filter(function(i) !i %in% c("]","["), strsplit(IUPACtoBase(x),"")[[1]])
        a[,i]<-1/length(i)-(4-length(i))*0.001
        a
    }))                      
}

counts2PWM<-function(counts){
    weights<-colSums(counts)
    mat<-matrix(0,ncol=4,nrow=nchar(names(weights)[1]))
    colnames(mat)<-c("A","C","G","T")
    mapply(function(name,weight){
        cs<-strsplit(name,"")[[1]]

        for(i in seq(length(cs))){
            mat[i,cs[i]]<<-mat[i,cs[i]]+weight
                                 }
    },
           names(weights),weights)
    mat
}

count2freq<-function(mat){
    round(t(apply(mat,1,function(x){
        freq<-x/sum(x)
        nil<-freq==0
        freq[nil]<-0.001
        freq[which.max(freq)]<-freq[which.max(freq)]-0.001*sum(nil)
        freq
    })),3)
}

homerProbs<-function(a,b,na,nb,a1,b1){
    fgmcr<-do.call("/",as.list(rowSums(a)))
    bgmcr<-do.call("/",as.list(rowSums(b)))
    fgc<-sum(a1)
    bgc<-sum(b1)
    fgr<-fgc/na
    bgr<-bgc/nb
    
    list(fgc=fgc,fgr=fgr,fgt=na,bgc=bgc,bgr=bgr,bgt=nb,fgmcr=fgmcr,bgmcr=bgmcr)
        
}

homerFind<-function(fg,bg,motif){
    a<-countMotifs(fg,motif)
    b<-countMotifs(bg,motif)
    a1<-countMotifsOne(fg,motif)
    b1<-countMotifsOne(bg,motif)
    na<-length(fg) # mot and compl
    nb<-length(bg)
    probs<-homerProbs(a,b,na,nb,a1,b1)
    pwm<-count2freq(counts2PWM(a))
    square<-rbind(c(probs$fgc,probs$fgt),c(probs$bgc,probs$bgt))
    score<-scoreMotifHomer(motif)
    test<-chisq.test(square)
    p.value<-test$p.value
    pstring=paste0("T:",probs$fgc,".0(",round(probs$fgr*100,2),"%)",
                   "B:",probs$bgc,".0(",round(probs$bgr*100,2),"%)",
                   "P:1e",floor(log(p.value,10)))
    list(motif=motif,name=motif,score=score,pvalue=log(p.value),pstring=pstring,pwm=pwm)
    
}

addMotifHeader.homer<-function(x,...){
    ord<-order(sapply(x,function(obj) obj$rank))
    lapply(x[ord],function(obj)with(obj,{
    paste0(head,paste(apply(pwm,1,paste,collapse="\t"),collapse="\n"),"\n")}))
}


motifs2Homer<-function(fg,bg,mots){
    summary<-lapply(mots,function(mot)homerFind(fg,bg,mot))
    names(summary)<-mots
    ranks<-rank(sapply(summary,function(x) x$pvalue),ties.method="first")
    names(ranks)<-mots
    headString<-lapply(mots,function(mot){
        summ<-summary[[mot]]

        paste0(">",summ$motif,"\t",
               ranks[mot],"-",summ$motif,"\t",
               round(summ$score,6),"\t",
               round(summ$pvalue,6),"\t",
               0,"\t",
               summ$pstring,"\n"
               )
    })
    names(headString)<-mots
    ret<-lapply(mots,function(mot){
        append(summary[[mot]],list(rank=ranks[mot],head=headString[mot]))
    })
    class(ret)<-"homer"
    ret
}

write.homer<-function(x,file){
    cat(paste0(addMotifHeader.homer(x)),file=file,sep="")
}

print.homer<-function(x,...){

        cat(paste0(addMotifHeader.homer(x)),sep="")
}


count2freq(counts2PWM(countMotifs(fg,"CANNTG")))


homerFind(fg,bg,mots[20])




writeXStringSet(ccca$fasta[ccca$reg[,"Leukemia"]],"~/Desktop/fg.fasta")
writeXStringSet(ccca$fasta[ccca$reg[,"NONE"]],"~/Desktop/bg.fasta")


###
library(kmerhmm)
library(Biostrings)
library(msa)

mers<-countNMers("~/Desktop/fixed.fasta",-1*ccca$prc$eigenVectors[ccca$reg[,"Leukemia"],1],8)

reg<-mers[,2]>median(mers[,2])+4*mad(mers[,2])


aligned<-msa(as.vector(mers[reg,1]),cluster="upgma",type="dna")

alignedChar<-apply(as.matrix(aligned),1,paste0,collapse="")

library(abind)

sequences<-as.matrix(aligned)
levels<-c("A","C","G","T","-")

kmerhmm::buildHMM(sequences,levels,50)

