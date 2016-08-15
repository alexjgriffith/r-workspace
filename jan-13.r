source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/jan-functions.r")
###source("~/Masters/CCCA/inst/scipts/eboxFrequency.r")
#source("~/Masters/mulcal/newR/rGREAT.r")
#source("~/Dropbox/R/makeLatexTable.R")


# generate variables
source("~/r-workspace/jan-variables.r")
fasta5<-getSeq(BSgenome.Hsapiens.UCSC.hg19,bed5$chr,start=bed5 $start+150,width=300)
fasta20<-getSeq(BSgenome.Hsapiens.UCSC.hg19,bed20$chr,start=bed20$start+150,width=300)
#source("~/r-workspace/fasta-data.r")



wreg20s1<-mwr(reg20s1,conv20s1)
wreg5s1<-mwr(reg5ps1,conv5s1)

wreg20<-mwr(reg20FSel(1),reg20FSel(0.25))

wreg5<-mwr(reg20FSel(1),reg20FSel(0.25))

allWrap(fasta20,fasta5,wreg20,wreg5,"1_0.25")

combn(unique(swapFunD(categories)))

temp<-t(combn(,2))
names<-rbind(temp,t(apply(temp,1,function(x) c(x[2],x[1]))))
names[order(names[,1]),]

temp<-unique(swapFunD(categories))
names<-t(do.call(cbind,lapply(temp,function(x) sapply(temp, function(y) c(x,y)))))

namesF<-names[names[,1]!=names[,2],]


sd="1_1"
for (fasta in list(list(fasta20,wreg20,20),list(fasta5,wreg5,5)))
    for(n in apply(namesF,1,list))
        for(l in c(8,6)){            
            a<-qhwJ(fasta[[1]],fasta[[2]],n[[1]][1],n[[1]][2],fasta[[3]],l,sd)
        }


allWrap<-function(fasta20,fasta5,wreg20,wreg5,sd){
a<-qhwJ(fasta20,wreg20,"Erythroid","All",20,8,sd)
a<-qhwJ(fasta20,wreg20,"Erythroid","None",20,8,sd)
a<-qhwJ(fasta20,wreg20,"Leukemia","All",20,8,sd)
a<-qhwJ(fasta20,wreg20,"Leukemia","None",20,8,sd)
a<-qhwJ(fasta20,wreg20,"HSC","All",20,8,sd)
a<-qhwJ(fasta20,wreg20,"HSC","None",20,8,sd)
a<-qhwJ(fasta20,wreg20,"ECFC","All",20,8,sd)
a<-qhwJ(fasta20,wreg20,"ECFC","None",20,8,sd)
a<-qhwJ(fasta20,wreg20,"Erythroid","All",20,6,sd)
a<-qhwJ(fasta20,wreg20,"Erythroid","None",20,6,sd)
a<-qhwJ(fasta20,wreg20,"Leukemia","All",20,6,sd)
a<-qhwJ(fasta20,wreg20,"Leukemia","None",20,6,sd)
a<-qhwJ(fasta20,wreg20,"HSC","All",20,6,sd)
a<-qhwJ(fasta20,wreg20,"HSC","None",20,6,sd)
a<-qhwJ(fasta20,wreg20,"ECFC","All",20,6,sd)
a<-qhwJ(fasta20,wreg20,"ECFC","None",20,6,sd)
a<-qhwJ(fasta5,wreg5,"Erythroid","All",5,8,sd)
a<-qhwJ(fasta5,wreg5,"Erythroid","None",5,8,sd)
a<-qhwJ(fasta5,wreg5,"Leukemia","All",5,8,sd)
a<-qhwJ(fasta5,wreg5,"Leukemia","None",5,8,sd)
a<-qhwJ(fasta5,wreg5,"HSC","All",5,8,sd)
a<-qhwJ(fasta5,wreg5,"HSC","None",5,8,sd)
a<-qhwJ(fasta5,wreg5,"ECFC","All",5,8,sd)
a<-qhwJ(fasta5,wreg5,"ECFC","None",5,8,sd)
a<-qhwJ(fasta5,wreg5,"Erythroid","All",5,6,sd)
a<-qhwJ(fasta5,wreg5,"Erythroid","None",5,6,sd)
a<-qhwJ(fasta5,wreg5,"Leukemia","All",5,6,sd)
a<-qhwJ(fasta5,wreg5,"Leukemia","None",5,6,sd)
a<-qhwJ(fasta5,wreg5,"HSC","All",5,6,sd)
a<-qhwJ(fasta5,wreg5,"HSC","None",5,6,sd)
a<-qhwJ(fasta5,wreg5,"ECFC","All",5,6,sd)
a<-qhwJ(fasta5,wreg5,"ECFC","None",5,6,sd)
}


stampDF<-function(combs,pvalue,sd,dir= "~/thesis-november/")
    do.call(rbind,apply(combs,1,function(x,pvalue) data.frame(file=paste(dir,x[1],"_pvalue=",pvalue,"_len=",x[2],"_sd=",sd,".motif",sep=""),compatator=x[3],name=x[1],size=x[2],pvalue=pvalue),pvalue ))

sd="1_0.25"
combs20<-stampDF(cbind(strsplit("Erythroid_None Erythroid_All Leukemia_None Leukemia_All HSC_None HSC_All ECFC_None ECFC_All"," ")[[1]],c(rep(6,8),rep(8,8)),c(rep("TRANSFAC_Fams",16),rep("JASPAR_Fams",16))),pvalue=20,sd=sd,"~/thesis-january/")
combs5<-stampDF(cbind(strsplit("Erythroid_None Erythroid_All Leukemia_None Leukemia_All HSC_None HSC_All ECFC_None ECFC_All"," ")[[1]],c(rep(6,8),rep(8,8)),c(rep("TRANSFAC_Fams",16),rep("JASPAR_Fams",16))),pvalue=5,sd=sd,"~/thesis-january/")
allCombs<-rbind(combs5,combs20)


## fix multiStamp
df<-allCombs

n<-dim(df)[1]
i=1

x<-stampWrapper(as.character(df$file[i]),df[i,2],df[i,3])
                
y<-cbind(x,MACS=df$pvalue[i],size=df$size[i],filename=df$file[i])
                
pvalue<-logPvalueFH(df$file[i])


motifs<-sapply(loadPWM(as.character(df$file[i]))[,1],function(x) substr(x,2,nchar(x)))

order<-c(unlist(sapply(motifs,function(x,y) which(x==y),unique(y$motif))))

names(order)<-NULL

names(motifs)<-NULL
names(pvalue)<-NULL

cbind(y,pvalue=c(sapply(pvalue[order],rep,5)),tmotif=c(sapply(motifs,rep,5)),order=c(sapply(order,rep,5)))
                                  })
            )


## fix logPvalueFH

filename<-"~/thesis-january/ECFC_All_pvalue=20_len=6_sd=1_0.25.motif"

logPvalueFH<-function(filename)
    sapply(sapply(strsplit(sapply(strsplit(do.call(rbind,sapply(loadPWM(as.character(filename))[,2],strsplit,","))[,3],":"),"[[",2),"T"),"[",1),function(x) abs(log10(as.numeric(x))))


db<-multiStamp(allCombs)
write.table(db,"~/thesis-january/motifAnnotations-c.table",quote=FA0LSE,row.names=FALSE)



sortBy<-function(df,c1){
    df[order(df[[c1]]),]
}

find<-function(value,lin)
    makeLogic(grep(value,lin,ignore.case = TRUE),length(lin))

unique(db$dataset[grep("ca[acgt][actg]tg",db$motif,ignore.case=TRUE)],c("genes","motif","dataset","order")],"order")



                   
sortBy(db[with(db,{makeLogic(grep("None",dataset),length(dataset)) & size==6 & rank<3 & comparator=="JASPAR_Fams" }),c("genes","motif","dataset","order")],"order")


sortBy(db[ find("gata",db$genes) & find("None",db$dataset) &db$rank<3 ,c("genes","motif","dataset","order","rank","escore")],"order")


sortBy(db[ find("fox",db$genes) & find("None",db$dataset) &db$rank<3 ,c("genes","motif","dataset","order","rank","escore")],"order")

db[db$motif=="AGAWAAAA",c("genes","motif","dataset","order","rank","escore")]

cols<-c("ann","motif","dataset","order","rank","escore","MACS","pvalue")

db[with(db,{find("All",dataset)& rank==1 & size==6 &order<5}),cols]

db$ann=db$genes

sortBy(db[with(db,{find("Erythroid",dataset) &find("All",dataset) & size==6}),cols],"order")


sortBy(db[with(db,{find("runt",ann) &find("All",dataset)  &rank<4 & size==6 &MACS==5} ),cols],"order")


buildFM<-function(df){
    with(df,{
        do.call(rbind,lapply(unique(motif),function(x) as.character(genes[x==motif])))
})}

buildList<-function(x,y){
    do.call(rbind,lapply(x,function(x) do.call(rbind,lapply(y,function(y) cbind(x,y)))))
}


buildFM<-function(treat,control,Msize,cutoff){
    temp<-apply(buildList(c("TRANSFAC_Fams","JASPAR_Fams"),seq(3)),1,function(x) 
    sortBy(db[with(db,{find(treat,dataset) &find(control,dataset)&comparator==x[1]  &rank==x[2] & size==Msize &MACS==cutoff} ),cols],"order"))
cbind(temp[[1]][,c("motif","order","pvalue")],do.call(cbind,lapply(temp,function(x) x[,c("ann","escore")])))
}

library(xlsx)

# macs score cut off of 5 only for motifs in motifs.xlsxj
wb<-loadWorkbook(path.expand("~/Dropbox/motifs.xlsx"))

for(treat in unique(swapFunD(categories)))
    for(control in c("None","All")){
        sn<-paste0(treat,"_",control,"_8")
        #removeSheet(wb,sn)
        bf<-buildFM(treat,control,8,5)
        ts<-createSheet(wb,sn)
        addDataFrame(bf,sheet=ts,row.names = FALSE)
    }

saveWorkbook(wb,path.expand("~/Dropbox/motifs.xlsx"))
