# all of this relies on jan-26 being loaded
loc1<-"chr17:49408052-49408752"
loc1b<-24227
bed20[loc1b,]

heights20[loc1b,]
## #> 10122
## grep chr17 eryt/eryt_unique_nodupes.bed | awk '$3>49408052 && $2<49408752{print $0}' | wc -l
## #> 4

genCommand<-function(x,name){
paste0("grep ",x$chr," ",name,"/",name,"_unique_nodupes.bed | awk '$3>",x$start," && $2<",x$end,"{print $0}' | wc -l\n")
}

system(cat(sapply(categories,function(x) genCommand(bed20[loc1b,],x))))


apply(heights20,2,order)[loc1b,]/dim(heights20)[1]

normData[locc1b,]


rev(sort(apply(sapply(seq(10),function(n) normalize(prc20$eigenVectors[,n]) ),1,function(x) sqrt(var(x)))))


tot<-prcomp(t(normData))$sdev

sum(tot[1:6])/sum(tot)
sum(tot[1:4]^2)/sum(tot^2)
# the first three principle components contain more thean half of the varience


# on cluster
rawdata<-(function(name)paste0("/mnt/brand01-00/mbrand_analysis/data_sets/",name,"/",name,"_unique_nodupes.bed"))(categories)

score<-pileUp(bed20,rawdata,22)
write.table(score,"~/Dropbox/udm_pvalue_20.txt",quote=FALSE,col.names = FALSE,row.names=FALSE)

## once transfered to workspace
tp20<-addColnames(read.table("~/Dropbox/udm_pvalue_20.txt"),categories)


