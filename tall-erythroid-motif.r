source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/jan-functions.r")
source("~/r-workspace/feb-functions.r")


qhw<-function(fasta,r,p1,p2,pvalue,len=8,sd=0){
    print(paste("~/feb-homer/",p1,"_",p2,"_pvalue=",pvalue,"_len=",len,"_sd=",sd,".motif",sep=""))
    homerWrapper(fasta,r[,p1],r[,p2],"~/Masters/mulcal/inst/lib/homer-4.7/bin/homer2",paste("~/feb-homer/",p1,"_",p2,"_pvalue=",pvalue,"_len=",len,"_sd=",sd,".motif",sep=""),opts=paste("-S 25 -len ",len," > /tmp/homerTrash 2>&1 ",sep=""))
} 

stampDF<-function(combs,pvalue,sd,dir="~/feb-homer/")
    do.call(rbind,apply(combs,1,function(x,pvalue) data.frame(file=paste(dir,x[1],"_pvalue=",pvalue,"_len=",x[2],"_sd=",sd,".motif",sep=""),compatator=x[3],name=x[1],size=x[2],pvalue=pvalue,sd=sd),pvalue ))


env<-addFasta(getPRC20(2))

LE<-qhw(env$fasta,env$reg,"Leukemia","Erythroid",20,8,"2_alt")
LE<-qhw(env$fasta,env$reg,"Leukemia","Erythroid",20,7,"2_alt")
LE<-qhw(env$fasta,env$reg,"Leukemia","Erythroid",20,6,"2_alt")
EL<-qhw(env$fasta,env$reg,"Erythroid","Leukemia",20,8,"2_alt")
EL<-qhw(env$fasta,env$reg,"Erythroid","Leukemia",20,7,"2_alt")
EL<-qhw(env$fasta,env$reg,"Erythroid","Leukemia",20,6,"2_alt")

combs<-stampDF(cbind(strsplit("Erythroid_Leukemia Leukemia_Erythroid" ," ")[[1]],c(rep(6,2),rep(7,2),rep(8,2)),c(rep("TRANSFAC_Fams",6),rep("JASPAR_Fams",6))),pvalue=20,sd="2_alt")
db<-multiStamp(combs)
write.table(db,"~/feb-homer/motifAnnotations_erythroid_jurkat.table",quote=FALSE,row.names=FALSE)




library(xlsx)
###
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


cols<-c("ann","motif","dataset","order","rank","escore","MACS","pvalue")
db<-read.table("~/feb-homer/motifAnnotations_erythroid_jurkat.table",header=TRUE)
db$ann<-db$genes

file<-path.expand("~/Dropbox/motifs_sd_Erythroid_Jurkat.xlsx")
wb<-createWorkbook()

for(tc in list(list("Erythroid","Leukemia"),list("Leukemia","Erythroid"))){
    treat=tc[[1]]
    control=tc[[2]]
    sn<-paste0(treat,"_",control)
    ts<-createSheet(wb,sn)
    srow=1
    for(size in c(6,7,8)){            
        bf<-buildFM(sn,sn,size,20)
        addDataFrame(bf,sheet=ts,row.names = FALSE,startRow = srow)
        srow=srow+dim(bf)[1]+2
        print(sn)
        print(srow)
        print(size)
    }
}

saveWorkbook(wb,file)
####

db<-read.table("~/feb-homer/motifAnnotations_20_2-alt.table",header=TRUE)
db$ann<-db$genes

file<-path.expand("~/Dropbox/motifs_sd_2_macs_20.xlsx")

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

saveWorkbook(wb,file)
