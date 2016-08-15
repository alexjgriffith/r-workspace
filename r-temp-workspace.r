
globalScore<-function(binSize,rawdata){
    regions<-do.call(rbind, apply(read.table("/data/binaries/BEDTools/genomes/human.hg19.genome")[1:24,],1, function(x,step) {y<-seq(1,as.numeric(x[2]),step); cbind(as.character(x[1]),as.character(y),as.character(y+step))} ,binSize))
    data<-hg19Sort(data.frame(chro=regions[,1],start=as.integer(regions[,2]), end=as.integer(regions[,3])))
    score<-pileUp(data,rawdata,n=22)
    score
}



printScore<-function(score,filename){write.table(score,filename,quote=FALSE, row.names=FALSE, col.names=FALSE)}

#cats<-read.table("/home/griffita/Dropbox/UTX-Alex/jan/catagories")
#cats<-as.data.frame(strsplit("cd133_mock cd34_mock cd34_new_mock cem_mock ecfc-tsa_mock ecfc_old_mock eryt_a_mock eryt_f_mock jurk_mock jurk_sandar_mock k562_mock_1 k562_mock_2 meka_mock rpmi_mock_1 tall_p1_mock tall_p2_mock tall_p3_mock"," "))

cats<-data.frame(strsplit("eryt cd34_new cem_1 jurk", " "))
prefix<-"/mnt/brand01-00/mbrand_analysis/data_sets/"
#suffix<-"_sorted.bed"
suffix<-"_unique_nodupes.bed"
rawdata<-apply(cats,1,function(x){paste(prefix,x,"/",x,suffix,sep="")})
binSizes=as.integer(10^seq(4 ,7) )

#for (b in binSizes){
    binSize=1000
    score<-globalScore(binSize,rawdata)
    #printScore(score,paste("global_analysis_bin_",binSize,".height",sep=""))
    printScore(score,paste("global_control_bin_",binSize,".height",sep=""))
#}

#for (binSize in binSizes){
regions<-do.call(rbind, apply(read.table("/data/binaries/BEDTools/genomes/human.hg19.genome")[1:24,],1, function(x,step) {y<-seq(1,as.numeric(x[2]),step); cbind(as.character(x[1]),as.character(y),as.character(y+step))} ,binSize))
    data<-hg19Sort(data.frame(chro=regions[,1],start=as.integer(regions[,2]), end=as.integer(regions[,3])))
printScore(data,paste("global_control_bin_",binSize,".front",sep=""))

#}
loc<-which(apply(score,1,sum)>0)


genes<-paste("E",rand(10))


inVec<-loadBedFile("peaks/jan/combined_mock.bed")
inVec<-loadBedFile("~/4x4-combined_unified_700.bed")
data<-hg19Sort(inVec)
#data<-hg19Sort(data.frame(chro=regions[,1],start=as.integer(regions[,2]), end=as.integer(regions[,3])))
score<-pileUp(data,rawdata,n=2)
printScore(score,paste("combined_control.height",sep=""))
printScore(data,paste("combined_control.front",sep=""))

write.table(cats,"../peaks/jan/contcats",quote=FALSE,row.names=FALSE,col.names=FALSE)

temp<-kmeans((function(x,cats){y<-t(x)%*%prcomp(t(x))$rotation; rownames(y)<-as.character(unlist(cats)); return(y)})(qn(score),cats)[,c(1,2)], 4)



normData=qn(score)
pcs<-prcomp(t(normData))
png("~/21x21-combined-pca1-3.png")
plotPCs(pcs,c(1,3),normData,as.character(unlist(cats)))
dev.off()

tdata<-t(apply(qn(score)[ord,],1,function(x) x/mean(x)))#/sqrt(var(x))))
tdata2<-tdata#<-t(apply(qn(score)[ord,],1,function(x) x/mean(x)))
    
clust<-kmeans(tdata,200,iter.max=100)
or1<-sapply(clust$cluster,function(x,y) y[x],order(getHeights(clust$cluster)))
ord<-order(or1)


## for small datasets only
rownames<-apply(as.matrix(data),1, function(x)(gsub(" ","", paste(x[1],":",x[2],"-",x[3],collapse="",sep=""),perl=TRUE ) ))



lrclust<-hclust(as.dist(1-cor(score)))
lrclust$labels<-as.character(unlist(cats))

png("~/21x21-combined-dend.png")
plot(lrclust,hang=-1)
dev.off()

xn<-tdata2[ord,lrclust$order]#/sqrt(var(x))))
colnames(xn)<-as.character(unlist(cats))
rownames(xn)<-seq(length(ord))#as.character(rownames[ord])

myMelt<-function(inm){
    rown<-as.integer(unlist(sapply(rownames(inm),rep,dim(inm)[2])))
    tdata<-data.frame(name=rown,id=colnames(inm),value=do.call(c,lapply(seq(dim(inm)[1]),function(pos)inm[pos,])))
    #tdata$name<-factor(tdata$name, levels=unique(tdata$name))
    tdata$id<-factor(tdata$id, levels=unique(tdata$id))
    tdata
}



t<-myMelt(xn)

png("~/21x21-200-hm.png",width=1080,height=1080)
p<-ggplot(t,aes(x=id,y=name))+geom_tile(aes(fill=value))+  scale_fill_gradient(low = "red",high = "green")
p+theme(axis.text.y=element_blank())
dev.off()


p<-p + theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
p+  scale_fill_gradient(low = "green",high = "red")
dev.off()

checkChroms<-function(data){
data[which(apply(cbind(data[2:dim(data)[1],"chro"],data[1:(dim(data)[1]-1),"chro"]),1,function(x) ! x[1]==x[2])),]}

hg19Sort<-function (data) 
{
    neworder<-Filter(function(x) x %in% levels(data$chro), strsplit("chr1 chr2 chr3 chr4 chr5 chr6 chr7 chrX chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr20 chrY chr19 chr22 chr21 chrM" ," ")[[1]]

     )
                                        #levels(inVec$chro)
    #neworder <- levels(data$chro)[c(1, 12, 16, 17, 18, 19, 20, 
    #    24, 21, 22, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 23, 14, 
    #    15, 25)]
    data$chro <- factor(data$chro, neworder)
    #sorter<-data.frame(order=seq(length(neworder)))
    #rownames(sorter)<-neworder
    #nord<-with(data,order(
    #                      as.integer(sapply(as.character(chro),function(x) sorter[x,])),start))
    nord<-with(data,order(as.integer(chro),start))
    data <- data[nord, ]

    data
    #sorter
}

(function(chrom,n=0){for(i in data$chro){if(i!=chrom){print(paste(i,as.character(n),sep=" "));n=0; chrom=i};n=n+1 }})("")
