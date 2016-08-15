source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/jan-functions.r")
###source("~/Masters/CCCA/inst/scipts/eboxFrequency.r")
source("~/Masters/mulcal/newR/rGREAT.r")
#source("~/Dropbox/R/makeLatexTable.R")


# generate variables
source("~/r-workspace/dec-variables.r")
source("~/r-workspace/jan-variables.r")

Leuk<-loadPWM("~/thesis-january/Leukemia_All_pvalue=5_len=6_sd=1_0.25.motif")
ECFC<-loadPWM("~/thesis-january/ECFC_All_pvalue=5_len=6_sd=1_0.25.motif")
Eryt<-loadPWM("~/thesis-january/Erythroid_All_pvalue=5_len=6_sd=1_0.25.motif")

ECFCMotifs<-strsplit("TTGTTT",split=" ")[[1]]
ErytMotifs<-strsplit("ATCCGA CCGTAT",split=" ")[[1]]
LeukMotifs<-strsplit("CAGACW ACCGCT AAAACA",split=" ")[[1]]




# source("http://bioconductor.org/biocLite.R")
# biocLite("seqLogo")

library(seqLogo)

psl<-function(data){
    seqLogo(makePWM(data),xaxis=FALSE,yaxis=FALSE,ic.scale=FALSE)
}

for(x in sapply(LeukMotifs,function(x) which(paste0(">",x)== Leuk[,1]))){
    png(paste0("~/Dropbox/UTX-Alex/br-data/motifs/",gsub(">","",Leuk[x,1][[1]]),"_logo.png"), bg = "transparent")
    psl(Leuk[x,3][[1]])
    dev.off()
}
for(x in sapply(ECFCMotifs,function(x) which(paste0(">",x)== ECFC[,1]))){
    png(paste0("~/Dropbox/UTX-Alex/br-data/motifs/",gsub(">","",ECFC[x,1][[1]]),"_logo.png"), bg = "transparent")
    psl(ECFC[x,3][[1]])
    dev.off()
}
for(x in sapply(ErytMotifs,function(x) which(paste0(">",x)== Eryt[,1]))){
    png(paste0("~/Dropbox/UTX-Alex/br-data/motifs/",gsub(">","",Eryt[x,1][[1]]),"_logo.png"), bg = "transparent")
    psl(Eryt[x,3][[1]])
    dev.off()
}





plot(prc5$eigenVectors[,c(1,3)])


ocup5<-read.table(mfn22(5),header=T)

ocup5Data<-ocup5[,4:dim(ocup5)[2]]

cats<-swapFunD(colnames(ocup5Data))

ords<-sapply(unique(cats), function(x) which(x == cats))

compOrds<-addColnames(do.call(cbind,lapply(ords,function(o) apply(cbind(ocup5Data[, o],0),1,function(x) {y<-sum(x) ; if(y>0) 1 else 0}))), names(ords))

unique<-apply(compOrds,1,function(x) if(sum(x)>1) FALSE else TRUE)



df<-as.data.frame(addColnames(prc5$eigenVectors[unique,c(1,3)],c("x","y")))

df$col=apply(compOrds[unique,],1,function(x,y) y[which(x==1)],colnames(compOrds))

colours=swapFunC(levels(as.factor(df[,3])))


p<-ggplot(df,aes(x=x,y=y,col=col))+ggplot2::geom_point()+theme_bw()+
    scale_color_manual(values=colours)+
    theme(axis.line = element_blank(),
          axis.text=element_blank(),
          axis.title=element_blank(),
          title=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())


##png(,bg="transparent")
ggsave("~/Dropbox/test.png",p,bg="transparent")


### find ebox logo for Leuk

Leuk<-reg5FSel(1)[,"Leukemia"]

bed5<-ifZeroShift(read.table(mfn16Make(5,heights="afs",control="combined"),header=T)[,1:3])
#fasta5<-getSeq(BSgenome.Hsapiens.UCSC.hg19,bed5$chr,start=bed5$start+150,width=300)

temp<-grep(IUPACtoBase("CANNTG"),fasta5[Leuk],value=TRUE)
v<-gregexpr(IUPACtoBase("CANNTG"),temp)

alphabet<-list(A=1,C=2,G=3,T=4,
               a=1,c=2,g=3,t=4)

tansA<-function(input)
    unlist(alphabet[input])

bfreq<-do.call(rbind,strsplit(unlist(mapply(function(string,start) sapply(start, function(x) substr(string,x,x+5)),temp,v)),""))

lr<-seq(dim(bfreq)[2])
ld<-seq(dim(bfreq)[1])

brapt<-matrix(rep(0,6*4),ncol=4,nrow=6)

for(i in ld)
    for(j in lr){
        k<-tansA(bfreq[i,j])
        brapt[j,k]= brapt[j,k]+1
    }


###call homer between leukemia and jurkat
a<-qhwJ(fasta5,reg5FSel(1),"Leukemia","Erythroid",5,6,1)


png(paste0("~/Dropbox/UTX-Alex/CACCTG-Leukemia-vs-Erythroid-logo.png"), bg = "transparent")
seqLogo(makePWM(a[5,3][[1]]),ic.scale=FALSE,xaxis=FALSE,yaxis=FALSE)
dev.off()

png(paste0("~/Dropbox/UTX-Alex/CACCTG-logo.png"), bg = "transparent")
seqLogo(makePWM(matrix(c(c(0,1,0,0),c(1,0,0,0),c(0,1,0,0),c(0,1,0,0),c(0,0,0,1),c(0,0,1,0)),ncol=6,nrow=4)),ic.scale=FALSE,xaxis=FALSE,yaxis=FALSE)
dev.off()
