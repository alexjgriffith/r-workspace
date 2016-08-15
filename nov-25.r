source("~/r-workspace/nov-functions.r")


rawDataFiles<-function(categories)paste("/mnt/brand01-00/mbrand_analysis/data_sets/",categories,"/",categories,"_unique_nodupes.bed",sep="")

makePeakFiles<-function(categories)paste("/mnt/brand01-00/mbrand_analysis/october/",categories,"_single_mock_peaks.xls",sep="")

fullAFS<-wrapAFS(makePeakFiles,categories)
pilotAFS<-wrapAFS(makePeakFiles,categoriesB)
pilot2AFS<-wrapAFS(makePeakFiles,c("k562_1","k562_2",categoriesB))

## for pilot
cs<-makeForkCluster(getOption("mc.cores",22),timeout=(4*60*60))

ret<-parallel::parLapply(cs,pvalues,pilotAFS)


for(i in seql(ret)){
   fileName<-paste("~/thesis-november/","4x4-single-pvalue=",pvalues[i],"_heights.matrix",sep="")
   matrixName<-paste("~/thesis-november/","4x4-single-pvalue=",pvalues[i],".matrix",sep="")   
   print(fileName)
   sdata<-hg19Sort(renameColumn(ret[[i]],"chr","chro"))
   score<-pileUp(sdata[,1:3],rawDataFiles(categories),20)   
   write.table(score,fileName,quote=FALSE,row.names=FALSE)
   write.table(sdata,matrixName,quote=FALSE,row.names=FALSE)
}

## for extended pilot

ret<-parallel::parLapply(cs,pvalues,pilot2AFS)


for(i in seql(ret)){
   fileName<-paste("~/thesis-november/","6x6-single-pvalue=",pvalues[i],"_heights.matrix",sep="")
   matrixName<-paste("~/thesis-november/","6x6-single-pvalue=",pvalues[i],".matrix",sep="")   
   print(fileName)
   sdata<-hg19Sort(renameColumn(ret[[i]],"chr","chro"))
   score<-pileUp(sdata[,1:3],rawDataFiles(categories),20,cs)   
   write.table(score,fileName,quote=FALSE,row.names=FALSE)
   write.table(sdata,matrixName,quote=FALSE,row.names=FALSE)
}


## for treatment

ret<-parallel::parLapply(cs,pvalues,fullAFS)




for(i in seql(ret)){
   fileName<-paste("~/thesis-november/","22x22-single-pvalue=",pvalues[i],"_heights.matrix",sep="")
   matrixName<-paste("~/thesis-november/","22x22-single-pvalue=",pvalues[i],".matrix",sep="")   
   print(fileName)
   sdata<-hg19Sort(renameColumn(ret[[i]],"chr","chro"))   
   score<-pileUp(sdata[,1:3],rawDataFiles(categories),20,clust=cs)   
   write.table(score,fileName,quote=FALSE,row.names=FALSE)
   write.table(sdata,matrixName,quote=FALSE,row.names=FALSE)
   cmd<-paste("scp",fileName,paste("griffita@ogic.ca:/data/websites/dropbox.ogic.ca/mbrandlab/alex_analysis/thesis/exdata/","22x22-single-pvalue=",pvalues[i],"_heights.matrix",sep=""))
   #system(cmd)
   cmd<-paste("scp",matrixName,paste("griffita@ogic.ca:/data/websites/dropbox.ogic.ca/mbrandlab/alex_analysis/thesis/exdata/","22x22-single-pvalue=",pvalues[i],".matrix",sep=""))
   #system(cmd)
}

#stopCluster(cs)

## for control
for(i in 1:30){
   fileName<-paste("~/thesis-november/","cont-single-pvalue=",pvalues[i],"_heights.matrix",sep="")
   matrixName<-paste("~/thesis-november/","22x22-single-pvalue=",pvalues[i],".matrix",sep="")
   ret<-read.table(matrixName,header =T)
   print(fileName)
   print(str(ret))
   sdata<-hg19Sort(ret)
   score<-pileUp(sdata[,1:3],rawDataFiles(conts),16,clust=cs)   
   write.table(score,fileName,quote=FALSE,row.names=FALSE)
   cmd<-paste("scp",fileName,paste("griffita@ogic.ca:/data/websites/dropbox.ogic.ca/mbrandlab/alex_analysis/thesis/exdata/","cont-single-pvalue=",pvalues[i],"_heights.matrix",sep=""))
}

stopCluster(cs)
