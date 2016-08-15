source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/jan-functions.r")
source("~/r-workspace/jan-variables.r")

for (sd in c(1,2,3,4))
for (i in c("Leukemia","ECFC","Erythroid","HSC")){
    filename=paste0("~/Dropbox/",i,"_pvalue_20","_sd_",sd,"_j26.bed")
    print (filename)
    write.table(bed20[reg20FSel(sd)[,i],],filename,quote=FALSE,col.names=FALSE,
                row.names = FALSE)
}




ds<-data.frame(name=c("tall_p1", "tall_p2_1", "tall_p2_2", "tall_p3_1", "tall_p3_2", "jurk_sandar_1", "jurk_sandar", "jurk", "rpmi_1", "rpmi_2", "cem_1", "cem_2", "ecfc-tsa", "meka", "cd133", "cd34", "cd34_new", "eryt", "eryt_f", "eryt_a", "k562_1", "k526_2"),type=c(1,1,1,1,1,1,1,1,1,1,1,1,2,3,3,3,3,4,4,4,4,4))


r="ecfc"
sd=1
n<-prc20$normData[reg20FSel(sd)[,r],swapFunD(colnames(prc20$normData,2))==r]
d<-prc20$normData[reg20FSel(sd)[,r],swapFunD(colnames(prc20$normData,2))!=r]
mean2<-function(x) sum(x)/(length(x)-1)
x<-apply(cbind(n,0),1,mean2)/apply(cbind(d,0),1,mean2)
hist(x)

plot(seq(0,5,0.2),sapply(seq(0,5,0.2),function(sd){
#a<-which(reg20FSel(sd)[,"Leukemia"])
a<-which(normalize(prc20$eigenVectors[,1])<(-sd))
b<-order(prc20$normData[,"eryt"],decreasing=TRUE)[1:length(a)]
length(intersect(b,a))/length(a)}),xlab = "", ylab = "")


plotPCMat2D(pca2Matr(prc20),c(1,4),categories,swapFun,swapFunB,swapFunC)

prc20$eigenVectors[,1]

normalize<-function(x)(x-mean(x))/(sqrt(var(x)))

eryt=bed20[normalize(prc20$eigenVectors[,1])>2,]
jurk=bed20[normalize(prc20$eigenVectors[,1])<(-2),]
stem=bed20[normalize(prc20$eigenVectors[,4])>2,]
ecfc=bed20[normalize(prc20$eigenVectors[,4])<(-2),]

printBB(eryt,"~/Dropbox/eryt_sd_2.bb")
printBB(jurk,"~/Dropbox/jurk_sd_2.bb")
printBB(stem,"~/Dropbox/stem_sd_2.bb")
printBB(ecfc,"~/Dropbox/ecfc_sd_2.bb")

printBB<-function(x,ofile){
tfile<-tempfile()
write.table(x,tfile,quote=FALSE,col.names=FALSE, row.names=FALSE)
system(paste0("~/bin/bedToBigBed ",tfile," ~/genomes/human.hg19.genome ",ofile))
}


# investigate using base R

which(over20[,"jurk"]==1)

normData<-qn(heights20)

prc<-prcomp(t(normData))$rotation

qp<-function(x,a,b){plot(x[,a],x[,b])}

qp(t(normData)%*%prc,1,2)


comp<-function(a,b){length(intersect(a,b))/length(a)}
x<-seq(0,5,0.2)
plot(x,sapply(x,function(sd){
a<-which(normalize(prc[,1])>sd)
b<-which(apply(cbind(over20[,categories[swapFunD(categories)=="Leukemia"]],0) ,1,sum)>0)
c<-order(normData[,"jurk"],decreasing = TRUE)[1:length(a)]
comp(b,c)
})
     ,xlab = "", ylab = "")


b<-which(apply(cbind(over20[,categories[swapFunD(categories)=="Leukemia"]],0) ,1,sum)>0)
c<-order(normData[,"jurk"],decreasing = TRUE)[1:length(b)]
comp(b,c)


n="cd34"
b<-which(over20[,n]>0)
c<-order(normData[,n],decreasing = TRUE)[1:length(b)]
comp(c,b)

mean(heights20[which(over20[,n]>0),n])/mean(heights20[which(over20[,n]==0),n])

sort(heights20[which(over20[,n]>0),n],decreasing=TRUE)

sort(heights20[order(normData[,n],decreasing = TRUE)[1:length(b)],n],decreasing=TRUE)

#a<-which(normalize((normData*prc[,1])[,"eryt"])>1)
a<-which(normalize(prc[,1])>2)
b<-which(apply(cbind(over20[,categories[swapFunD(categories)=="Erythroid"]],0) ,1,sum)>0)
#b<-which(over20[,"eryt"]==1)
length(intersect(a,b))/length(a)

colmean<-function(type)
    apply(cbind(heights20[,categories[swapFunD(categories)==type]],0) ,1,mean2)

selcol<-function(type,fun="==")
which(apply(cbind(over20[,categories[do.call(fun,list(type,swapFunD(categories)))]],0) ,1,sum)>0)

order(colmean("Leukemia")[selcol("Leukemia")],decreasing=TRUE)

a<-order(colmean("Leukemia"),decreasing=TRUE)#[selcol("Leukemia",fun="==")]
b<-order(colmean("Leukemia"),decreasing=TRUE)#[selcol("Leukemia",fun="!=")]

comp<-function(x){paste0(x$chr,":",x$start,"-",x$end)}

data.frame(
    locs=comp(bed20[a,][1:100,]),
    Leuk=round(colmean("Leukemia")[a][1:100],0),
    Eryt=round(colmean("Erythroid")[a][1:100],0),
    HSC=round(colmean("HSC")[a][1:100],0),
    ECFC=round(colmean("ECFC")[a][1:100],0),
    PC1=floor(normalize(prc20$eigenVectors[a,1][1:100])),
    PC4=floor(normalize(prc20$eigenVectors[a,4][1:100]))
    )

    Eryt=round(colmean("Erythroid")[b][1:100],0),
    HSC=round(colmean("HSC")[b][1:100],0),
    ECFC=round(colmean("ECFC")[b][1:100],0)
    )

normData
