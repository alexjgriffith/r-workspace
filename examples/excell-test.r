library(xlsx)

xlsx::read.xlsx2("")

load("~/thesis-january/ebox-overlap.RData")


#xlsx::
write.xlsx2(motF$pca5A$matrix,path.expand("~/Dropbox/testfile.xlsx"),sheetName="pca5A",)

write.xlsx2(do.call(cbind,lapply(motF,function(x){as.data.frame(x$matrix/x$colSize)})),path.expand("~/Dropbox/testfile.xlsx"),sheetName="pca5A",)

opts<-lapply(motF,function(x){as.data.frame(t(t(x$matrix)/x$colSize))})

ropts<-lapply(opts,round,3)


eboxFile<-path.expand("~/Dropbox/testfile.xlsx")

lapply(names(opts),function(x) write.xlsx2(opts[[x]],eboxFile,sheetName=x,append=TRUE))

write.xlsx2("",eboxFile,sheetName="Sheet1",append=FALSE)

ework<-loadWorkbook(eboxFile)
#removeSheet(ework,"Sheet1")
removeSheet(ework,"ebox")

s<-createSheet(ework,"ebox")

j=1;
k=1;
for(i in seq(length(names(opts)))){
    if(i==5){
        j=1
        k=13
    }
    sc<-1+6*(j-1)
    n=names(opts)[i]
    addDataFrame(ropts[[n]],sheet=s,startColumn = sc,startRow = k)
    addDataFrame(as.data.frame(n),sheet=s,startColumn = sc,startRow=k,row.names=F,col.names=F)
    j=j+1
}

saveWorkbook(ework,eboxFile)

dataSetSize<-t(do.call(cbind,lapply(motF,function(x)unlist(x$colSize))))

motifSize<-t(do.call(cbind,lapply(motF,function(x)unlist(x$rowSize))))


ework<-loadWorkbook(eboxFile)
removeSheet(ework,"eboxinfo")
s<-createSheet(ework,"eboxinfo")
addDataFrame(dataSetSize,sheet=s,startColumn = 1,startRow = 1)
addDataFrame(motifSize,sheet=s,startColumn = 6,startRow = 1)
saveWorkbook(ework,eboxFile)

