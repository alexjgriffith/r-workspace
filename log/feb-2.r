source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/jan-functions.r")
source("~/r-workspace/feb-functions.r")
source("~/r-workspace/feb-variables.r")


dend(heights[,categories],n=4,colours=c("red","orange","blue","green"),swapFun=swapFunM)

options(max.print=200)


fullAFS<-wrapAFS(makePeakFiles,categories)

ret<-fullAFS(20)


source("~/r-workspace/nov-functions.r")
ret<-read.table("~/test_density.bed",header=T)
rawDataFiles<-function(categories)paste("/mnt/brand01-00/mbrand_analysis/data_sets/",categories,"/",categories,"_sorted.bed",sep="")
#getPeakDensity(rawDataFiles(categories[15]),ret,as.character(ret$chr), length(ret$chr),700)
#score1<-pileUp(ret[[1]],rawDataFiles(categories),20,TRUE)
score7<-peakDensity(ret,rawDataFiles(categories),20,700,TRUE)

categories<-readCategories("~/Dropbox/UTX-Alex/jan/catagories")
#### on cluster
#1 make new afs for pvalues from 5-75

# redefine categories
makeLotsOfAFSs<-function(control){
    
fullAFS<-wrapAFS(makePeakFiles,categories,control)
pvalues<-seq(from=5,to=75,by=2.5)
n=length(pvalues)
cs<-makeForkCluster(n,renice=0)
ret<-parallel::parLapply(cs,pvalues,fullAFS)
stopCluster(cs)

for(i in seql(pvalues)){
    fname=paste0("~/Dropbox/Data/AFS/22_pvalue_",pvalues[i],"_control_",control,".txt")
    write.table(ret[[i]],fname,row.names=FALSE,quote=FALSE)
}
}

null<-lapply(c("single","no","combined"),makeLotsOfAFSs)


readAFS<-function(pvalue=5,control="single",directory="~/Dropbox/Data/AFS/",forward=22){
    fname<-paste0(directory,forward,"_pvalue_",pvalue,"_control_",control,".txt")
    read.table(fname,header=T)
}

readUDM<-function(pvalue=5,tr="treatment",control="single",directory="~/Dropbox/Data/AFS/",forward=22){
    fname<-paste0(directory,forward,"_",tr,"_pvalue_",pvalue,"_control_",control,".txt")
    read.table(fname,header=T)
}

### make UDM for pvalues from 5-75 for no combined and single using treatment and control save it in ~/Dropbox/Data/UDM/udm_treatment_pvalue_5_control_single.txt
rawDataFiles<-function(categories)paste("/mnt/brand01-00/mbrand_analysis/data_sets/",categories,"/",categories,"_sorted.bed",sep="")



makeALLUDM<-function(cs){
copts<-list(categories)
cnames<-c("treatment")

for(cats in seql(copts)){
    for(cont in c("single","combined","no"))
        for(pvalue in pvalues){
            ret<-readAFS(pvalue,cont)
            fname=paste0("~/Dropbox/Data/UDM/22_",cnames[cats],"_pvalue_",pvalue,"_control_",cont,".txt")
            score<-pileUp(data=orderBed(ret),rawDataFiles(copts[[cats]]),n=20, clust=cs)
            write.table(addColnames(score,copts[[cats]]),fname,quote=FALSE,row.names=FALSE)
            print(fname)
        }    
}

}

#categories<-readCategories("~/Dropbox/Data/feb16-categories.txt")
pvalues<-seq(from=5,to=75,by=2.5)

cs<-makeForkCluster(20,renice=0)
#score<-pileUp(data=orderBed(ret),rawDataFiles(categories),n=20, clust=cs)
makeALLUDM(cs)
stopCluster(cs)




categories<-readCategories("~/Dropbox/Data/feb16-categories-b.txt")
fullAFS<-wrapAFS(makePeakFiles,categories,"single")
a<-fullAFS(5)
write.table(a,paste0("~/Dropbox/Data/AFS/23_pvalue_",5,"_control_","single",".txt"),row.names=FALSE,quote=FALSE)
a<-fullAFS(20)
write.table(a,paste0("~/Dropbox/Data/AFS/23_pvalue_",20,"_control_","single",".txt"),row.names=FALSE,quote=FALSE)
fullAFS<-wrapAFS(makePeakFiles,categories,"combined")
a<-fullAFS(5)
write.table(a,paste0("~/Dropbox/Data/AFS/23_pvalue_",5,"_control_","combined",".txt"),row.names=FALSE,quote=FALSE)
a<-fullAFS(20)
write.table(a,paste0("~/Dropbox/Data/AFS/23_pvalue_",20,"_control_","combined",".txt"),row.names=FALSE,quote=FALSE)
fullAFS<-wrapAFS(makePeakFiles,categories,"no")
a<-fullAFS(5)
write.table(a,paste0("~/Dropbox/Data/AFS/23_pvalue_",5,"_control_","no",".txt"),row.names=FALSE,quote=FALSE)
a<-fullAFS(20)
write.table(a,paste0("~/Dropbox/Data/AFS/23_pvalue_",20,"_control_","no",".txt"),row.names=FALSE,quote=FALSE)


makeALLUDM2<-function(cs){
copts<-list(categories)
cnames<-c("treatment")
for(cats in seql(copts)){
    for(cont in c("single","combined","no"))
        for(pvalue in c(5,20)){
            ret<-readAFS(pvalue,cont,forward="23")
            fname=paste0("~/Dropbox/Data/UDM/23_",cnames[cats],"_pvalue_",pvalue,"_control_",cont,".txt")
            score<-pileUp(data=orderBed(ret),rawDataFiles(copts[[cats]]),n=20, clust=cs)
            write.table(addColnames(score,copts[[cats]]),fname,quote=FALSE,row.names=FALSE)
            print(fname)
        }    
}

}

cs<-makeForkCluster(20,renice=0)

#score<-pileUp(data=orderBed(ret),rawDataFiles(categories),n=20, clust=cs)
makeALLUDM2(cs)
stopCluster(cs)


makeALLUDM3<-function(cs){
copts<-list(conts)
cnames<-c("control")
for(cats in seql(copts)){
    for(cont in c("combined"))
        for(pvalue in c(5,20)){
            ret<-readAFS(pvalue,cont,forward="22")
            fname=paste0("~/Dropbox/Data/UDM/22_",cnames[cats],"_pvalue_",pvalue,"_control_",cont,".txt")
            score<-pileUp(data=orderBed(ret),rawDataFiles(copts[[cats]]),n=20, clust=cs)
            write.table(addColnames(score,copts[[cats]]),fname,quote=FALSE,row.names=FALSE)
            print(fname)
        }    
}

}

cs<-makeForkCluster(16,renice=0)

#score<-pileUp(data=orderBed(ret),rawDataFiles(categories),n=20, clust=cs)
makeALLUDM3(cs)

stopCluster(cs)


makeALLUDM4<-function(cs){
cnames<-c("meka","tall_p1","rpmi","ecfc","k562_1")    
copts<-lapply(cnames,function(x) categories[categories!=x])
for(cats in seql(copts)){
    for(cont in c("combined"))
        for(pvalue in c(5,20)){
            fullAFS<-wrapAFS(makePeakFiles,categories,"no")
            ret<-fullAFS(pvalue)
            write.table(ret,paste0("~/Dropbox/Data/AFS/",cnames[cats],"_pvalue_",pvalue,"_control_","no",".txt"),row.names=FALSE,quote=FALSE)
            #ret<-readAFS(pvalue,cont,forward="22")
            fname=paste0("~/Dropbox/Data/UDM/",cnames[cats],"_",cnames[cats],"_pvalue_",pvalue,"_control_",cont,".txt")
            #score<-pileUp(data=orderBed(ret),rawDataFiles(copts[[cats]]),n=20, clust=cs)
            #write.table(addColnames(score,copts[[cats]]),fname,quote=FALSE,row.names=FALSE)
            print(fname)
        }    
}

}


cs<-makeForkCluster(20,renice=0)

#score<-pileUp(data=orderBed(ret),rawDataFiles(categories),n=20, clust=cs)
makeALLUDM4(cs)

stopCluster(cs)

