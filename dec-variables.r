source("~/r-workspace/dec-functions.r")
heights20<-read.table(mfn22(pvalues[9],TRUE),header=T)
heights5<-read.table(mfn22(pvalues[3],TRUE),header=T)
prc20<-pca(heights20)
prc5<-pca(heights5)


opt20s<-"PC1+1:PC1-1:PC4+1:PC4-1:PC2+1:PC2-1"
opt5s<-"PC1+1:PC1-1:PC3+1:PC3-1:PC5+1:PC5-1:PC7+1:PC7-1"

opt20<-strsplit(opt20s,":")[[1]]
opt5<-strsplit(ot5s,":")[[1]]


reg20<-qPeakP(prc20$eigenVectors,opt20s)
reg5<-qPeakP(prc5$eigenVectors,opt5s)

HSC5<-addColnames(cbind(reg5[,"PC3+1"]&reg5[,"PC7+1"],
                  reg5[,"PC3-1"]&reg5[,"PC7-1"])|(reg5[,"PC3-1"]&reg5[,"PC7+1"])|(reg5[,"PC3+1"]&reg5[,"PC7-1"]),c("HSC","NotHSC"))

Other5<-addColnames(cbind(reg5[,"PC3+1"]&reg5[,"PC5+1"],
                  reg5[,"PC3-1"]&reg5[,"PC5-1"])|(reg5[,"PC3-1"]&reg5[,"PC5+1"])|(reg5[,"PC3+1"]&reg5[,"PC5-1"]),c("Other","NotOther"))


bed20<-ifZeroShift(read.table("~/thesis-november/22x22-pvalue=20.matrix",header=T)[,1:3])
bed5<-ifZeroShift(read.table("~/thesis-november/22x22-pvalue=5.matrix",header=T)[,1:3])

#fasta20<-getSeq(BSgenome.Hsapiens.UCSC.hg19,bed20$chro,start=bed20$start+150,width=300)

#fasta5<-getSeq(BSgenome.Hsapiens.UCSC.hg19,bed5$chro,start=bed5$start+150,width=300)
#fasta5<-read

motifFiles<-c(motifFileNames(opt20,pairSwitch(opt20),rep(20,6),rep(8,6)),
  motifFileNames(opt20,pairSwitch(opt20),rep(20,6),rep(6,6)),
  motifFileNames(opt5,pairSwitch(opt5),rep(5,8),rep(8,8)),
  motifFileNames(opt5,pairSwitch(opt5),rep(5,8),rep(6,8)),
  motifFileNames(c("HSC","Other"),c("NotHSC","NotOther"),rep(5,2),rep(6,2)),
  motifFileNames(c("HSC","Other"),c("NotHSC","NotOther"),rep(5,2),rep(8,2))
  )


