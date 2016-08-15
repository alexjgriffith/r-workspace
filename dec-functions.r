library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(httr)


#quick homer wrapper
qhw<-function(fasta,r,p1,p2,pvalue,len=8){
    homerWrapper(fasta,r[,p1],r[,p2],"~/Masters/mulcal/inst/lib/homer-4.7/bin/homer2",paste("~/thesis-november/",p1,"_",p2,"_pvalue=",pvalue,"_len=",len,".motif",sep=""),opts=paste("-S 25 -len ",len,sep=""))
}


qPeakP<-function(eigen,commands)
    addColnames(applyPeakPartitions(eigen,gPObj(commands)),strsplit(commands,":")[[1]])

#' generate partition objext
#' @export
#' @examples
#' gPObj("PC1-1:PC2+1:PC2-7")
#' x<-seq(-10,10)/2
#' reg<-as.list(paste("PC",1,"+",seq(-10,10)/2,sep=""))
#' gPObj(do.call(paste,append(reg,list(sep=":")))))
#' y<-apply(applyPeakPartitions(prc20$eigenVectors,obj,2,sum))
#' plot(x,y/dim(prc20$eigenVectors)[1],xlab="SD",ylab="#peaks")
gPObj<-function(string){
    gObjAux<-function(string,ingnore){
        if(length(ingnore)<2)
            return(gsub(ingnore[1],"",string))    
        else
            return(gObjAux(gsub(ingnore[1],"",string),ingnore[2:length(ingnore)]))
    }
    gObjAux2<-function(gtsp,operations,swapFun){
        l<-sapply(strsplit(gtsp,"")[[1]],function(x) x %in% operations)
        a<-unique(c(0,which(l),(nchar(gtsp)+1)))
        n<-do.call(rbind,createZip(Filter(function(x) x, unlist(apply(cbind(a[1:length(a)-1]+1,a[2:length(a)]),1,function(x) if(x[1]!=x[2])x else FALSE)))))
        num<-apply(n,1,function(x) substr(gtsp,x[1],x[2]-1))
        ope<-sapply(which(l),function(x) substr(gtsp,x,x))
        t<-""
        if(length(ope)==2){
            t<-ope[2]
        }
        list(as.numeric(num[1]),swapFun(ope[1]),as.numeric(paste(t,num[2],sep="")))
    }
    ingnore<-c("PC"," ",",","&")
    end<-":"
    operations<-c("+","-","w","o")
    swapFun<-simpleSwapFun("+ gt - lt w w o o")
    so<-unlist(strsplit(gObjAux(string,ingnore),":"))
    lapply(so,gObjAux2,operations,swapFun)
}

ifZeroShift<-function(data){
    locZ<-which(data$start==0)
    if(length(locZ)>0){
        data[locZ,"start"]<-data[locZ,"start"]+1
        data[locZ,"end"]<-data[locZ,"end"]+1
    }
    data
}

mrb<-function(reg)
    sample(c(rep(TRUE,length(which(reg))),rep( FALSE,length(which(! reg)))))

stampDF<-function(combs,pvalue,dir= "~/thesis-november/")
    do.call(rbind,apply(combs,1,function(x,pvalue) data.frame(file=paste(dir,x[1],"_pvalue=",pvalue,"_len=",x[2],".motif",sep=""),compatator=x[3],name=x[1],size=x[2],pvalue=pvalue),pvalue ))


logPvalueFH<-function(filename)
    as.numeric(sapply(sapply(strsplit(sapply(strsplit(do.call(rbind,sapply(loadPWM(as.character(filename))[,2],strsplit,","))[,3],":"),"[[",2),"T"),"[",1),function(x) substr(x,4,nchar(x))))

multiStamp<-function(df){
    n<-dim(df)[1]
    do.call(rbind,
            lapply(seq(n),function(i){x<-stampWrapper(as.character(df$file[i]),df[i,2],df[i,3])
                                      y<-cbind(x,MACS=df$pvalue[i],size=df$size[i],filename=df$file[i],sd=df$sd[i])
                                      pvalue<-c(sapply(logPvalueFH(df$file[i]),rep,5))
                                      motifs<-sapply(loadPWM(as.character(df$file[i]))[,1],function(x) substr(x,2,nchar(x)))
                                      order<-c(unlist(sapply(motifs,function(x,y) which(x==y),unique(y$motif))))
                                      names(order)<-NULL
                                      names(motifs)<-NULL
                                      names(pvalue)<-NULL
                                      cbind(y,pvalue=c(sapply(pvalue[order],rep,5)),tmotif=c(sapply(motifs,rep,5)),order=c(sapply(order,rep,5)))
                                  })
            )
}
      
stampWrapper<-function(motifFile,comparator,name=motifFile){
    postFormB<-function(motif,comparator){
        motif_file<-path.expand(motif)
        list(input="",  MAX_FILE_SIZE="10000000",  motiffile=upload_file(motif_file),  mtrim="on",  ICbox="0.4",  num_hits="5",  match_db=comparator,  MAX_FILE_SIZE="1000000",  metric="PCC",  align="SWU",  mulalign="IR",  tree="UPGMA",  submit="Submit",  metric_des="Pearson%27s+Correlation+Coefficient%3A%0D%0A-------------------------------%0D%0A%28Recommended%29%0D%0AMin%3A+-1%2C+Max%3A+%2B1+%28per+column%29.%0D%0A ", align_des="Ungapped+Smith-Waterman%3A%0D%0A------------------------------%0D%0A%28Recommended%29%0D%0ALocal+alignment+for+comparisons+between+short+motifs.",  mulalign_des="Iterative+Refinement%3A%0D%0A------------------------------%0D%0AUpdates+multiple+alignment+by+leave-one-out+iteration.",  tree_des="UPGMA%3A%0D%0A------------------------------%0D%0A%28Recommended%29%0D%0ADistance-based+tree+construction.+")
    }
    ua<- "Mozilla/5.0 (Windows NT 6.1; WOW64; rv:33.0) Gecko/20100101 Firefox/33.0"
    body<-postFormB(motifFile,comparator)
    a<-POST("http://www.benoslab.pitt.edu/stamp/run.php" ,body=body,user_agent(ua))
    datalink<-strsplit(strsplit(grep("parent.location",strsplit(content(a,"text"),"\n")[[1]],value=TRUE)," ")[[1]][5],"'")[[1]][2]
    b<-GET(datalink)
    motifData<-data.frame()
    for(i in (seq(100)-1)){
        c<-grep(paste("\"motif",i,"\"",sep=""),strsplit(content(b,"text"),"\n")[[1]],value=TRUE)
        if(length(c)==0)
            break
        eScores<-as.numeric(unlist(lapply(grep("width=\"75\"",strsplit(c,"</TD>")[[1]],value=TRUE),function(x) strsplit(x,">")[[1]][2]))[2:6])
        genes<-unlist(lapply(grep("width=\"175\"",strsplit(c,"</TD>")[[1]],value=TRUE),function(x) strsplit(x,">")[[1]][2]))[2:6]
    motif<-unlist(strsplit(grep("<strong>",strsplit(c,"</strong></")[[1]],value=TRUE),">"))
        t<-data.frame(motif=motif[length(motif)],rank=seq(5),genes=genes,escore=eScores,dataset=name,comparator=comparator)
        motifData<-rbind(motifData,t)
    }
    motifData
}


pairSwitch<-function(x){
    is.even<-function(x)
        (x-floor(x/2)*2)==0
    
    out<-rep("",length(x))
    for( i in seq(length(x))){
        if (is.even(i)){
            out[i-1]=x[i]
        }
        else
            out[i+1]=x[i]
        }
    return(out)
}


opt20<-strsplit("PC1+1:PC1-1:PC4+1:PC4-1:PC2+1:PC2-1",":")[[1]]
opt5<-strsplit("PC1+1:PC1-1:PC3+1:PC3-1:PC5+1:PC5-1:PC7+1:PC7-1",":")[[1]]

motifFileNames<- function(p1,p2,pvalue,len)
    paste("~/thesis-november/",p1,"_",p2,"_pvalue=",pvalue,"_len=",len,".motif",sep="")

motifFiles<-c(motifFileNames(opt20,pairSwitch(opt20),rep(20,6),rep(8,6)),
  motifFileNames(opt20,pairSwitch(opt20),rep(20,6),rep(6,6)),
  motifFileNames(opt5,pairSwitch(opt5),rep(5,8),rep(8,8)),
  motifFileNames(opt5,pairSwitch(opt5),rep(5,8),rep(6,8)),
  motifFileNames(c("HSC","Other"),c("NotHSC","NotOther"),rep(5,2),rep(6,2)),
  motifFileNames(c("HSC","Other"),c("NotHSC","NotOther"),rep(5,2),rep(8,2))
  )

sel2<-function(y){
    function(x)
        y[,x]
}
