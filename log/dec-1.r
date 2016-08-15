# pvalue = 20
# 1+ = Erythroid
# 1- = Leukemia
# 4+ = HSC + Megakaryocyte
# 4- = ECFC

source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/dec-variables.r")





multiPlot(list(plotPCMat2D(pca2Matr(prc20),c(pcsO[1],pcsO[4]),categories,swapFun,swapFunB,swapFunC),
               plotPCMat2D(pca2Matr(prc20),c(pcsO[1],pcsO[2]),categories,swapFun,swapFunB,swapFunC),
               plotPCMat2D(pca2Matr(prc5),c(pcsO[1],pcsO[3]),categories,swapFun,swapFunB,swapFunC),
               plotPCMat2D(pca2Matr(prc5),c(pcsO[3],pcsO[7]),categories,swapFun,swapFunB,swapFunC)
               ),2,2)


a<-qhw(fasta20,reg20,"PC1+1","PC1-1",20)
a<-qhw(fasta20,reg20,"PC1-1","PC1+1",20)
a<-qhw(fasta20,reg20,"PC2+1","PC2-1",20)
a<-qhw(fasta20,reg20,"PC2-1","PC2+1",20)
a<-qhw(fasta20,reg20,"PC4+1","PC4-1",20)
a<-qhw(fasta20,reg20,"PC4-1","PC4+1",20)
a<-qhw(fasta5,reg5,"PC1+1","PC1-1",5)
a<-qhw(fasta5,reg5,"PC1-1","PC1+1",5)
a<-qhw(fasta5,reg5,"PC3+1","PC3-1",5)
a<-qhw(fasta5,reg5,"PC3-1","PC3+1",5)
rHSC<-addColnames(cbind(reg5[,"PC3+1"]&reg5[,"PC7+1"],
                  reg5[,"PC3-1"]&reg5[,"PC7-1"])|(reg5[,"PC3-1"]&reg5[,"PC7+1"])|(reg5[,"PC3+1"]&reg5[,"PC7-1"]),c("HSC","NotHSC"))
a<-qhw(fasta5,rHSC,"HSC","NotHSC",5)
rOther<-addColnames(cbind(reg5[,"PC3+1"]&reg5[,"PC5+1"],
                  reg5[,"PC3-1"]&reg5[,"PC5-1"])|(reg5[,"PC3-1"]&reg5[,"PC5+1"])|(reg5[,"PC3+1"]&reg5[,"PC5-1"]),c("Other","NotOther"))
a<-qhw(fasta5,rOther,"Other","NotOther",5)
a<-qhw(fasta20,reg20,"PC1+1","PC1-1",20,6)
a<-qhw(fasta20,reg20,"PC1-1","PC1+1",20,6)
a<-qhw(fasta20,reg20,"PC2+1","PC2-1",20,6)
a<-qhw(fasta20,reg20,"PC2-1","PC2+1",20,6)
a<-qhw(fasta20,reg20,"PC4+1","PC4-1",20,6)
a<-qhw(fasta20,reg20,"PC4-1","PC4+1",20,6)
a<-qhw(fasta5,reg5,"PC1+1","PC1-1",5,6)
a<-qhw(fasta5,reg5,"PC1-1","PC1+1",5,6)
a<-qhw(fasta5,reg5,"PC3+1","PC3-1",5,6)
a<-qhw(fasta5,reg5,"PC3-1","PC3+1",5,6)
a<-qhw(fasta5,HSC5,"HSC","NotHSC",5,6)
a<-qhw(fasta5,Other5,"Other","NotOther",5,6)


combs20<-stampDF(cbind(strsplit("PC1-1_PC1+1 PC1+1_PC1-1 PC2+1_PC2-1 PC2-1_PC2+1 PC4+1_PC4-1 PC4-1_PC4+1"," ")[[1]],c(rep(6,6),rep(8,6)),c(rep("TRANSFAC_Fams",12),rep("JASPAR_Fams",12))),pvalue=20)

combs5<-stampDF(cbind(strsplit("PC1-1_PC1+1 PC1+1_PC1-1 PC3+1_PC3-1 PC3-1_PC3+1 Other_NotOther HSC_NotHSC"," ")[[1]],c(rep(6,6),rep(8,6)),c(rep("TRANSFAC_Fams",12),rep("JASPAR_Fams",12))),pvalue=5)

allCombs<-rbind(combs5,combs20)

db<-multiStamp(allCombs)

write.table(db,"~/thesis-november/motifAnnotations.table",quote=FALSE,row.names=FALSE)


