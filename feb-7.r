source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/jan-functions.r")
source("~/r-workspace/feb-functions.r")

db=read.table("~/thesis-feb/motifAnnotations.table",header=TRUE)

with(db,{data.frame(motif,genes)[dataset=="Erythroid_NONE" & MACS==20 & comparator=="TRANSFAC_Fams" ,]})

env<-getConEboxInfo(getEboxInfo(addFasta(getPRC20(2))))
a<-findEboxBoth(env)
rm(env)
env<-getConEboxInfo(getEboxInfo(addFasta(getPRC5(2))))
b<-findEboxBoth(env)
rm(env)
motF<-append(a,b)

# save eboxs and motifs to excel
opts<-lapply(motF,function(x){as.data.frame(x$matrix/apply(x$matrix,2,sum))})

library(xlsx)

(Compose(unlist,unique,length))(env$eboxLoc)

length(intersect(unique(unlist(env$eboxLoc)),which(env$Erythroid.alt)))/length(which(env$Erythroid.alt))

normalizeEbox<-function(t){
    x<-sapply(seq(dim(t$matrix)[2]),
              function(i){
                  c(t$matrix[,i]/t$colSize[1,i],
                    (t$colSize[1,i]-sum(t$matrix[,i]))/t$colSize[1,i])
              })
    addNames(x,colnames(t$matrix),c(rownames(t$matrix),"No-Ebox"))    
}

t=a$pcaA
normalizeEbox(t)

outer(apply(env$reg,1,which),md,Vectorize(function(a,b){intersect(a,b)}))

sum(sapply(md,function(x) length(intersect(x,which(env$Erythroid.alt))))/length(which(env$Erythroid.alt)))

opts<-lapply(list(a$bedD,a$pcaD,b$bedD,b$pcaD),function(x){as.data.frame(round(normalizeEbox(x),3))})

library(xlsx)

with(list(),{
    eboxFile<-path.expand("~/Dropbox/testfile2.xlsx")
    wb<-xlsx::createWorkbook()
    s<-xlsx::createSheet(wb,"ebox")
    p20<-as.data.frame(round(normalizeEbox(a$pcaD),3))
    p5<-as.data.frame(round(normalizeEbox(b$pcaD),3))
    xlsx::addDataFrame(p20,sheet=s,startColumn = 1,startRow = 1)
    xlsx::addDataFrame(p5,sheet=s,startColumn = 10,startRow = 1)
    xlsx::saveWorkbook(wb,eboxFile)
})



sortBy<-function(df,c1){
    df[order(df[[c1]]),]
}

find<-function(value,lin)
    makeLogic(grep(value,lin,ignore.case = TRUE),length(lin))


buildList<-function(x,y){
    do.call(rbind,lapply(x,function(x) do.call(rbind,lapply(y,function(y) cbind(x,y)))))
}

buildFM<-function(treat,control,Msize,cutoff){
    temp<-apply(buildList(c("TRANSFAC_Fams","JASPAR_Fams"),seq(3)),1,function(x) 
    sortBy(db[with(db,{find(treat,dataset) &find(control,dataset)&comparator==x[1]  &rank==x[2] & size==Msize &MACS==cutoff} ),cols],"order"))
cbind(temp[[1]][,c("motif","order","pvalue")],do.call(cbind,lapply(temp,function(x) x[,c("ann","escore")])))
}



db$ann<-db$genes
cols<-c("ann","motif","dataset","order","rank","escore","MACS","pvalue")

wb<-createWorkbook()

for(treat in unique(swapFunD(categories)))
    for(control in c("NONE","ALL")){
        sn<-paste0(treat,"_",control)
        ts<-createSheet(wb,sn)
        srow=1
        for(size in c(6,7,8)){            
            bf<-buildFM(treat,control,size,20)            
            addDataFrame(bf,sheet=ts,row.names = FALSE,startRow = srow)
            srow=srow+dim(bf)[1]+2
            print(sn)
            print(srow)
            print(size)
        }
    }

saveWorkbook(wb,path.expand("~/Dropbox/motifs_sd_1_macs_20.xlsx"))
