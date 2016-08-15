source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/dec-variables.r")
source("~/r-workspace/fasta-data.r")
source("~/r-workspace/jan-functions.r")
source("~/Masters/CCCA/inst/scipts/eboxFrequency.r")
source("~/Masters/mulcal/newR/rGREAT.r")
source("~/Dropbox/R/makeLatexTable.R")


# 1st step print AFS for combined single and no

AFSFilenames<-paste("~/Dropbox/UTX-Alex/br-data/afs/","afs_","pvalue=",rep(c(20,5),3),"_",sapply(c("combined","single","no"),rep,2),"_mock.txt",sep="")

# write table wrapper
wtw<-function(x,file,...)
    write.table(x,file,row.names=FALSE,quote=FALSE,...)

wtw((wrapAFS(makePeakFiles,categories,"combined"))(20),AFSFilenames[1])
wtw((wrapAFS(makePeakFiles,categories,"combined"))(5),AFSFilenames[2])
wtw((wrapAFS(makePeakFiles,categories,"single"))(20),AFSFilenames[3])
wtw((wrapAFS(makePeakFiles,categories,"single"))(5),AFSFilenames[4])
wtw((wrapAFS(makePeakFiles,categories,"no"))(20),AFSFilenames[5])
wtw((wrapAFS(makePeakFiles,categories,"no"))(5),AFSFilenames[6])

readAFSFile<-function(pvalue,control){
    if( pvalue %in% c(5,20) & control %in% c("combined","no","single"))
        return (read.table(paste("~/Dropbox/UTX-Alex/br-data/afs/","afs_","pvalue=",pvalue,"_",control,"_mock.txt",sep=""),header=T))
    else
        warning("Either ",pvalue," is not in c(5,20), or ",control," is not in c(\"combined\",\"no\",\"single\")")    
}






# generate files with regions (just do combined for now)
pca(addColnames(read.table(mfn16Make(20),header=T),categories))

plotPCMat2D(pca2Matr(pca(addColnames(read.table(mfn16Make(20),header=T),categories))),c("PC1","PC2"),categories,swapFun,swapFunB,swapFunC)

applyPeakPartitions(pca(addColnames(read.table(mfn16Make(20),header=T),categories))$eigenVectors,gPObj(opt20s))

colnameUR(genReg5pCombined(addColnames(applyPeakPartitions(pca(addColnames(read.table(mfn16Make(5),header=T),categories))$eigenVectors,gPObj(opt5s)),opt5 )),conv5p)
