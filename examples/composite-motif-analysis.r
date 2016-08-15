source("~/r-workspace/CCCA.r")
source("~/r-workspace/project.r")
source("~/r-workspace/preferedDistance.r")

source("~/Masters/CCCA/inst/scripts/eboxFrequency.r")

source("~/R/x86_64-unknown-linux-gnu-library/3.2/CCCA/scripts/eboxFrequency.r")

library(BSgenome.Hsapiens.UCSC.hg19)

# prefered distance pca
eboxs<-c(genEboxCombs(),"CANNTG")

gata<-"GATAA"

normalize<-CCCA::normalize

env20<-getPRC20(2)

env20<-addFasta(env20)

grepMotifs(motif,env20$fasta)

mList<-c(eboxs,gata)

cList<-lapply(mList,compliment)

locationsM<-lapply(mList,grep,env20$fasta)
locationsC<-lapply(cList,grep,env20$fasta)

save(locationsM,locationsC,env20,mList,cList,file="~/Dropbox/env20motifs.RData")

load("~/Dropbox/env20motifs.RData")

#make sure CCCA is updated
gnnebox<-findMotifDist(env20$fasta,"GATAA","CACCTG",env20$reg[,"Leukemia"],13,14,15)

plotFastMot(env20,mList,cList,locationsM,locationsC,1,12,reg=c("ECFC","Erythroid","HSC","Leukemia"),xlim=c(-32,32))
