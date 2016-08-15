source("~/r-workspace/CCCA.r")
source("~/r-workspace/project.r")


plotPC<-function(x,PC1,PC2,filt,...)plotPCMat2D(pca2Matr(x$prc),c(PC1,PC2),x$categories,Compose(filt,swapFun),Compose(filt,swapFunB),swapFunC,...)



env<-getPRC20(2)

colnames(env$heights)<-env$categories
colnames(env$over)<-c("chr","start","end",env$categories)

env$name<-"combined"
envSingle<-getPRC(pvalue=20,control="single",treatment="treatment")
envSingle$name<-"single"
envNo<-getPRC(pvalue=20,control="no",treatment="treatment")
envNo$name<-"no"

control<-getPRC(pvalue=20,control="combined",treatment="control")
control$categories<-gsub("\\.","-",colnames(control$heights))
control$name<-"control"

envTALL<-removeDataSets( "tall_p1",env)
envK562<-removeDataSets(c("k562_1", "k562_2"),env)
envMEKA<-removeDataSets("meka",env)
envCD34<-removeDataSets(c("cd34","cd34_new") ,env)

envECFC<-removeDataSets("ecfc-tsa" ,env)


names<-strsplit("combined single no control tall k562 meka cd34 ecfc"," ")[[1]]

data<-list(env,envSingle,envNo,control,envTALL,envK562,envMEKA,envCD34,envECFC)

test<-function(){ print("mocktonorm")
                   mockToNorm}

names<-"ecfc"
data<-list(envECFC)

mapply(function(name,e){
    print(name)
    print(e$name)
    print(e$categories)
    name14<-paste0("~/Dropbox/UTX-Alex/Paper/PC14_20_",name,".png")
    p14<-plotPC(e,1,4,ifelse(name=="control",test(),pass),label=TRUE,blank=TRUE)
    ggsave(name14,p14)
    name12<-paste0("~/Dropbox/UTX-Alex/Paper/PC12_20_",name,".png")
    p12<-plotPC(e,1,2,ifelse(name=="control",mockToNorm,pass),label=TRUE,blank=TRUE)
    ggsave(name12,p12)  
},names,data)



data.frame(names,size=sapply(data,function(x) dim(x$prc$normData)[1]))


saveToPaperPCA<-function(x){
    fname<-paste0("~/Dropbox/UTX-Alex/Paper/Analysis/PCA_Peaks_",x,".bed")
    write.table(env$bed[env$reg[,x],],fname,row.names=FALSE,col.names=TRUE,quote=FALSE)
}

lapply(colnames(env$reg),saveToPaperPCA)

saveToPaperBED<-function(x){
    fname<-paste0("~/Dropbox/UTX-Alex/Paper/Analysis/BED_Peaks_",x,".bed")
    write.table(bed[[x]],fname,row.names=FALSE,col.names=TRUE,quote=FALSE)}

saveToPaperBEDALL<-function(x){
    fname<-paste0("~/Dropbox/UTX-Alex/Paper/Analysis/ALL_BED_Peaks_",x,".bed")
    write.table(bedALL[[x]],fname,row.names=FALSE,col.names=TRUE,quote=FALSE)}

merged<-mergeFun(env$over[,4:dim(env$over)[2]],swapFunD)

unique<-apply(merged,1,sum)==1

bed<-apply(merged,2,function(x) env$bed[as.logical(x)&unique,])

bedALL<-apply(merged,2,function(x) env$bed[as.logical(x),])

lapply(names(bed),saveToPaperBED)

lapply(names(bedALL),saveToPaperBEDALL)
