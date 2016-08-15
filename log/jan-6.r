source("~/r-workspace/jan-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-variables.r")
source("~/Masters/CCCA/inst/scipts/eboxFrequency.r")
source("~/Masters/mulcal/newR/rGREAT.r")
source("~/Dropbox/R/makeLatexTable.R")


fasta20<-getSeq(BSgenome.Hsapiens.UCSC.hg19,bed20$chro,start=bed20$start+150,width=300)
over20=read.table(mfn22(pvalues[9]),header=T)
motifLoc<-sapply(genEboxCombs(),grepMotifs,fasta20)
peakLoc<-lapply(seq(4,dim(over20)[2]),Compose(sel2(over20),as.logical,which))
uniqueLocs<-apply(over20[,4:22],1,sum)==1
regLoc20<-lapply(seq(dim(reg20)[2]),Compose(sel2(reg20),which))



temp<-swapFunB(categories)
conOver20<-do.call(cbind,lapply(lapply(unique(temp),function(x) do.call(cbind,lapply(which(x==temp), function(x) over20[,4:dim(over20)[2]][,x]))),function(x) apply(cbind(x,0),1,function(x) if(sum(x)>0)1 else 0 )))
colnames(conOver20)<-unique(temp)
conUniqueLocs20<-apply(over20[,4:22],1,sum)==1
conPeakLoc<-lapply(unique(temp),function(x) do.call(unionN,lapply(which(x==temp), function(x) peakLoc[[x]])))


allData<-data.frame(count=do.call(rbind,lapply(motifLoc,length)))

allDataFreq=matrix(round(allData$count/length(fasta20),3),nrow=10,dimnames=list(as.list(genEboxCombs()),list("Count") ))

cat(makeTable(allDataFreq,"|c|r|",label="tab:ebox-all-20",caption="The frequency of E-Box occurances in the peaks that were limited to those with a MACS score cut off of 20."),file="ebox-all-20.table")

bedEboxUniqueFreq<-motifUniqueFrequency(genEboxCombs(),colnames(conOver20),motifLoc,conPeakLoc,conUniqueLocs20)

cat(makeTable(as.matrix(round(bedEboxUniqueFreq,3)),c("|l|ccccc|"),label="tab:ebox-bed-20",caption="The frequency of E-Box occurances in the subsets of peaks priveded by each of the cellular contexts , the peaks were limited to those with a MACS score cut off of 20.")),file="ebox-bed-20.table")

bedPCAFreq<-motifAllFrequency(genEboxCombs(),colnames(reg20),motifLoc,lapply(seq(dim(reg20)[2]),function(x) which(reg20[,x])))

cat(makeTable(as.matrix(round(bedPCAFreq,3)),c("|l|cccccc|"),label="tab:ebox-pca-20",caption="The frequency of E-Box occurances in the subsets of peaks identified using the principle components of the unified density matrix , the peaks were limited to those with a MACS score cut off of 20."),file="ebox-pca-20.table")


fasta5<-getSeq(BSgenome.Hsapiens.UCSC.hg19,bed5$chro,start=bed5$start+150,width=300)
over5=read.table(mfn22(pvalues[3]),header=T)
motifLoc<-sapply(genEboxCombs(),grepMotifs,fasta5)
peakLoc<-lapply(seq(4,dim(over5)[2]),Compose(sel2(over5),as.logical,which))
uniqueLocs<-apply(over5[,4:22],1,sum)==1
regLoc5<-lapply(seq(dim(reg5)[2]),Compose(sel2(reg5),which))

temp<-swapFunB(categories)
conOver5<-do.call(cbind,lapply(lapply(unique(temp),function(x) do.call(cbind,lapply(which(x==temp), function(x) over5[,4:dim(over5)[2]][,x]))),function(x) apply(cbind(x,0),1,function(x) if(sum(x)>0)1 else 0 )))
colnames(conOver5)<-unique(temp)
conUniqueLocs5<-apply(over5[,4:22],1,sum)==1
conPeakLoc<-lapply(unique(temp),function(x) do.call(unionN,lapply(which(x==temp), function(x) peakLoc[[x]])))


allData<-data.frame(count=do.call(rbind,lapply(motifLoc,length)))
allDataFreq=matrix(round(allData$count/length(fasta5),3),nrow=10,dimnames=list(as.list(genEboxCombs()),list("Count") ))
cat(makeTable(allDataFreq,"|c|r|",label="tab:ebox-all-5",caption="The frequency of E-Box occurances in the peaks that were limited to those with a MACS score cut off of 5."),file="ebox-all-5.table")

bedEboxUniqueFreq<-motifUniqueFrequency(genEboxCombs(),colnames(conOver5),motifLoc,conPeakLoc,conUniqueLocs5)
cat(makeTable(as.matrix(round(bedEboxUniqueFreq,3)),c("|l|ccccc|"),label="tab:ebox-bed-5",caption="The frequency of E-Box occurances in the subsets of peaks priveded by each of the cellular contexts , the peaks were limited to those with a MACS score cut off of 5."),file="ebox-bed-5.table")

bedPCAFreq<-motifAllFrequency(genEboxCombs(),colnames(reg5),motifLoc,lapply(seq(dim(reg5)[2]),function(x) which(reg5[,x])))
cat(makeTable(as.matrix(round(bedPCAFreq,3)),c("|l|cccccc|"),label="tab:ebox-pca-5",caption="The frequency of E-Box occurances in the subsets of peaks identified using the principle components of the unified density matrix , the peaks were limited to those with a MACS score cut off of 5."),file="ebox-pca-5.table")


## motifs=~/thesis-november/motifAnnotations.table
motifData="~/thesis-november/motifAnnotations.table"
motifTable<-read.table(motifData,col.names=TRUE)

# Gene Association

geneFile<-"/home/agriffith/Dropbox/UTX-Alex/jan/hg19.RefSeqGenes.csv"
geneList<-read.delim(geneFile)
chrom<-as.character(geneList$chrom)
tss<-as.numeric(geneList$txStart)
strand<-geneList$strand
# double check to make sure that levels(strand) > c("-","+")
levels(strand)<-c(-1,1)
strand<-as.numeric(strand)
statsAndData<-readSplitTable(heightFile)
regions<-genomicRegions(chrom,tss,strand,1000,5000,1000000)








plotPCMat2D(t(prc5$normData)%*%prc5$eigenVectors,c(1,5),categories,swapFun,swapFunB,swapFunC)



id="name"
genes20pca<-geneMatrix(over20,reg20,regions,geneList,id)
lapply(colnames(genes20pca),writegenes,genes20pca,paste("~/Dropbox/genes20pca-",id,"-",sep=""),conv20)

genes20bed<-geneMatrix(over20[conUniqueLocs20,c(1,2,3)],addColnames(do.call(cbind,lapply(colnames(conOver20),function(x) conOver20[conUniqueLocs20,x]==1)),colnames(conOver20) ),regions,geneList,id)
lapply(colnames(genes20bed),writegenes,genes20bed,paste("~/Dropbox/genes20bed-",id,"-",sep=""),pass)

genes5pca<-geneMatrix(over5,reg5p,regions,geneList,id)
lapply(colnames(genes5pca),writegenes,genes5pca,paste("~/Dropbox/genes5pca-",id,"-",sep=""),conv5)

genes5bed<-geneMatrix(over5[conUniqueLocs5,c(1,2,3)],addColnames(do.call(cbind,lapply(colnames(conOver5),function(x) conOver5[conUniqueLocs5,x]==1)),colnames(conOver5) ),regions,geneList,id)
lapply(colnames(genes5bed),writegenes,genes5bed,paste("~/Dropbox/genes5bed-",id,"-",sep=""),pass)
    
id="name2"
genes20pca<-geneMatrix(over20,reg20,regions,geneList,id)
lapply(colnames(genes20pca),writegenes,genes20pca,paste("~/Dropbox/genes20pca-",id,"-",sep=""),conv20)

genes20bed<-geneMatrix(over20[conUniqueLocs20,c(1,2,3)],addColnames(do.call(cbind,lapply(colnames(conOver20),function(x) conOver20[conUniqueLocs20,x]==1)),colnames(conOver20) ),regions,geneList,id)
lapply(colnames(genes20bed),writegenes,genes20bed,paste("~/Dropbox/genes20bed-",id,"-",sep=""),pass)

genes5pca<-geneMatrix(over5,reg5p,regions,geneList,id)
lapply(colnames(genes5pca),writegenes,genes5pca,paste("~/Dropbox/genes5pca-",id,"-",sep=""),conv5)

genes5bed<-geneMatrix(over5[conUniqueLocs5,c(1,2,3)],addColnames(do.call(cbind,lapply(colnames(conOver5),function(x) conOver5[conUniqueLocs5,x]==1)),colnames(conOver5) ),regions,geneList,id)
lapply(colnames(genes5bed),writegenes,genes5bed,paste("~/Dropbox/genes5bed-",id,"-",sep=""),pass)


    


