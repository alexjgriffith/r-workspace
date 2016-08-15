#fasta5<-getSeq(BSgenome.Hsapiens.UCSC.hg19,bed5$chr,start=bed5 $start+150,width=300)
#fasta20<-getSeq(BSgenome.Hsapiens.UCSC.hg19,bed20$chr,start=bed20$start+150,width=300)

over20<-read.table(mfn22(pvalues[9]),header=T)
eboxLoc20<-sapply(genEboxCombs(),grepMotifs,fasta20)
peakLoc20<-lapply(seq(4,dim(over20)[2]),Compose(sel2(over20),as.logical,which))
uniqueLocs20<-apply(over20[,4:22],1,sum)==1
regLoc20<-lapply(seq(dim(reg20)[2]),Compose(sel2(reg20),which))

conOver20<-do.call(cbind,lapply(lapply(unique(swapFunD(categories)),function(x) do.call(cbind,lapply(which(x==swapFunD(categories)), function(x) over20[,4:dim(over20)[2]][,x]))),function(x) apply(cbind(x,0),1,function(x) if(sum(x)>0)1 else 0 )))

colnames(conOver20)<-unique(swapFunD(categories))
conUniqueLocs20<-apply(over20[,4:22],1,sum)==1
conPeakLoc20<-lapply(unique(swapFunD(categories)),function(x) do.call(unionN,lapply(which(x==swapFunD(categories)), function(x) peakLoc20[[x]])))
names(conPeakLoc20)<-colnames(conOver20)

over5<-read.table(mfn22(pvalues[3]),header=T)
eboxLoc5<-sapply(genEboxCombs(),grepMotifs,fasta5)
peakLoc5<-lapply(seq(4,dim(over5)[2]),Compose(sel2(over5),as.logical,which))
uniqueLocs5<-apply(over5[,4:22],1,sum)==1
regLoc5<-lapply(seq(dim(reg5)[2]),Compose(sel2(reg5),which))


conOver5<-do.call(cbind,lapply(lapply(unique(swapFunD(categories)),function(x) do.call(cbind,lapply(which(x==swapFunD(categories)), function(x) over5[,4:dim(over5)[2]][,x]))),function(x) apply(cbind(x,0),1,function(x) if(sum(x)>0)1 else 0 )))
colnames(conOver5)<-unique(swapFunD(categories))
conUniqueLocs5<-apply(over5[,4:22],1,sum)==1
conPeakLoc5<-lapply(unique(swapFunD(categories)),function(x) do.call(unionN,lapply(which(x==swapFunD(categories)), function(x) peakLoc5[[x]])))
names(conPeakLoc5)<-colnames(conOver5)
