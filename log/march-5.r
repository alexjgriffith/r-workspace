
install.packages("devtools")
source("http://bioconductor.org/biocLite.R")
#pkgs <- rownames(installed.packages())
#biocLite()

biocLite("Biobase")

biocLite("Biostrings")


install.packages("httr")

install_github("alexjgriffith/CCCA")


library(abind)

load("peakDensity.RData")



l<-do.call(abind,append(
    lapply(which(swapFunB(env$categories)=="Leukemia"),function(x)pd20[[x]])
   ,list(along=3)))

acu<-l[,,1]
for(i in 2:12){
    acu<-acu+l[,,i]
}



dist<-function(reg){
    d<-which(swapFunD(env$categories)==reg)
l<-do.call(abind,append(
    lapply(d,function(x)pd20[[x]])
   ,list(along=3)))
    eryt<-l[,,1]
    if(length(d)>1)
        for(i in 2:length(d)){
            eryt<-eryt+l[,,i]
        }
    eryt
}

contexts<-unique(swapFunD(env$categories))

regs<-addNames(lapply(contexts,dist),
               contexts,
               list=TRUE)


d<-lapply(regs,function(x) x[env$reg[,"Leukemia"],])

o<-order(env$prc$eigenVectors[env$reg[,"Leukemia"],1],decreasing = TRUE)[1:1000]

lapply(contexts,
       function(x) {
           png(paste0("~/Dropbox/Data/dens/pca_20_2_Leukemia_",x,".png"))
           heatmap(d[[x]][o,],Rowv=NA,Colv=NA,margin=c(0,0))
           dev.off()
       })


       


d<-lapply(regs,function(x) x[env$reg[,"Erythroid"],])
o<-order(env$prc$eigenVectors[env$reg[,"Erythroid"],1],decreasing = TRUE)[1:1000]
lapply(contexts,
       function(x) {
           png(paste0("~/Dropbox/Data/dens/pca_20_2_Erythroid_",x,".png"))
           heatmap(d[[x]][o,],Rowv=NA,Colv=NA,margin=c(0,0))
           dev.off()
       })
d<-lapply(regs,function(x) x[env$reg[,"ECFC"],])
o<-order(env$prc$eigenVectors[env$reg[,"ECFC"],1],decreasing = TRUE)[1:1000]
lapply(contexts,
       function(x) {
           png(paste0("~/Dropbox/Data/dens/pca_20_2_ECFC_",x,".png"))
           heatmap(d[[x]][o,],Rowv=NA,Colv=NA,margin=c(0,0))
           dev.off()
       })     
d<-lapply(regs,function(x) x[env$reg[,"HSC"],])
o<-order(env$prc$eigenVectors[env$reg[,"HSC"],1],decreasing = TRUE)[1:1000]
lapply(contexts,
       function(x) {
           png(paste0("~/Dropbox/Data/dens/pca_20_2_HSC_",x,".png"))
           heatmap(d[[x]][o,],Rowv=NA,Colv=NA,margin=c(0,0))
           dev.off()
       })

heatmap(d[["ECFC"]][o,],Rowv=NA,Colv=NA,margin=c(0,0))
       
setOne<-function(x,v)
    (x-min(x))/(max(x)-min(x))





layer<-function(reg){
    d<-lapply(regs,function(x) x[env$reg[,reg],])
    o<-order(env$prc$eigenVectors[env$reg[,reg],1],decreasing = TRUE)
nd<-addNames(lapply(contexts,function(x) apply(d[[x]][o,],2,function(x) if(sum(x)>0)return(setOne(x)) else return(x))),
             contexts,
             list=TRUE
             )
#plot(-10000,xlim=c(0,700),ylim=c(0,max(sapply(nd,function(x)max(apply(x,2,mean))))))
#lines(apply(nd[["Erythroid"]],2,mean),col="red")
#lines(apply(nd[["Leukemia"]],2,mean),col="blue")
#lines(apply(nd[["ECFC"]],2,mean),col="green")
#lines(apply(nd[["HSC"]],2,mean),col="orange")
    infun<-function(x)
        cbind(apply(x,2,mean),
              apply(x,2,function(x) sort(x)[floor(length(x)*.25)]),
              apply(x,2,function(x) sort(x)[floor(length(x)*.75)])
              #apply(x,2,function(x) mean(x)-sd(x)/4),
              #apply(x,2,function(x) mean(x)+sd(x)/4)
              )
#cbind(as.data.frame(
    #as.data.frame(cbind(x=seq(700),addColnames(do.call(cbind,lapply(nd,infun)),paste0(sapply(names(nd),rep,3),".",rep(c("x","low","high"),4)))))
    cbind(as.data.frame(cbind(x=seq(-350,length=700),
                        addColnames(do.call(rbind,lapply(nd,infun)),c("y","low","high"))
                        )),t=c(sapply(swapFunC(contexts),function(x) rep(x,700))))
                    #),cbind(x=seq(700)))
}







library("ggplot2")


df<-layer("Erythroid")
h <- ggplot(df, aes(x,y,colour=t))
h + geom_ribbon(aes(ymin = low, ymax = high,fill=t),alpha=0.1) +
geom_line()+theme_bw()+theme(legend.position="none")+scale_color_manual(values=sort(swapFunC(contexts)))+scale_fill_manual(values=sort(swapFunC(contexts)))
    
        
