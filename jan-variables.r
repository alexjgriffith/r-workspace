over20<-addColnames(read.table(mfn16Make(20,heights="afs",control="combined"),header=T)[,1:length(categories)+3],categories)
#over5<-ifZeroShift(read.table(mfn16Make(5,heights="afs",control="combined"),header=T)[,1:3])

heights20<-addColnames(read.table(mfn16Make(20,control="combined"),header=T),categories)
heights5<-addColnames(read.table(mfn16Make(5,control="combined"),header=T),categories)
prc20<-pca(heights20)
prc5<-pca(heights5)


opt20all<-"PC1+%var%:PC1-%var%:PC4+%var%:PC4-%var%:PC2+%var%:PC2-%var%"
opt5all<-"PC1+%var%:PC1-%var%:PC3+%var%:PC3-%var%:PC5+%var%:PC5-%var%:PC7+%var%:PC7-%var%"



opt20s<-gop(opt20all,1)
opt5s<-gop(opt20all,1)


opt20<-strsplit(opt20s,":")[[1]]
opt5<-strsplit(opt5s,":")[[1]]


#reg20<-qPeakP(prc20$eigenVectors,opt20s)
#reg5<-qPeakP(prc5$eigenVectors,opt5s)
#reg5p<-genReg5pCombined3sd(reg5)

#reg20s1<-qPeakP(prc20$eigenVectors,gop(opt20all,1))
#reg5s1<-qPeakP(prc5$eigenVectors,gop(opt5all,1))
#reg5ps1<-genReg5pCombined1sd(reg5s1)

## replace these with function calls -->>



bed20<-ifZeroShift(read.table(mfn16Make(20,heights="afs",control="combined"),header=T)[,1:3])
bed5<-ifZeroShift(read.table(mfn16Make(5,heights="afs",control="combined"),header=T)[,1:3])

## these are big, load them manualy
#fasta20<-getSeq(BSgenome.Hsapiens.UCSC.hg19,bed20$chr,start=bed20$start+150,width=300)
#fasta5<-getSeq(BSgenome.Hsapiens.UCSC.hg19,bed5$chr,start=bed5$start+150,width=300)


