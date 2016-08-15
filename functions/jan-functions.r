reg20FSel<-function(sd)applySel(qPeakP(prc20$eigenVectors,gop(opt20all,sd)),convPvar(swapS20)(sd) )

reg5FSel<-function(sd)applySel(genReg5pCombined(qPeakP(prc5$eigenVectors,gop(opt5all,sd)),sd) ,convPvar(swapS5)(sd))

gop<-function(x,sd)gsub("%var%",sd,x)

orderBed<-function(ret)
    ret[order(as.character(ret[,1]),ret[,3]),]


normalize<-function(x)(x-mean(x))/(sqrt(var(x)))

printBB<-function(x,ofile){
tfile<-tempfile()
write.table(x,tfile,quote=FALSE,col.names=FALSE, row.names=FALSE)
system(paste0("~/bin/bedToBigBed ",tfile," ~/genomes/human.hg19.genome ",ofile))
}

makeLogic<-function(loc,size){
    x=rep(FALSE,size)
    x[loc]<-TRUE
    x
}
mwr<-function(regA,regB){
    cbind(regA,
          All=rep(TRUE,dim(regA)[1]),
          None=makeLogic(setdiff(seq(dim(regB)[1]),unlist(apply(regB,2,which))),dim(regB)[1]))
}

# redefined from november-functions
qhwJ<-function(fasta,r,p1,p2,pvalue,len=8,sd=1){
    #paste("~/thesis-january/",p1,"_",p2,"_pvalue=",pvalue,"_len=",len,"_sd=",sd,".motif,sep="")
    homerWrapper(fasta,r[,p1],r[,p2],"~/Masters/mulcal/inst/lib/homer-4.7/bin/homer2",paste("~/thesis-january/",p1,"_",p2,"_pvalue=",pvalue,"_len=",len,"_sd=",sd,".motif",sep=""),opts=paste("-S 25 -len ",len,sep=""))
}


#' has been add to CCCA/R/motifComparision
#'
#' @export
motifDistances<-function(loc,fasta,motifsPre=NULL,width=150){
    allLoc<-sort(unique(unlist(loc)))
    if(is.null(motifsPre))
        motifsPre=names(loc)
    motifs<-sapply(motifsPre,IUPACtoBase)
    compl<-sapply(motifs,compliment);
    shift<-function(x,width) min(abs(x-width))
    findEbox<-function(ebox,width){
        mlo<-gregexpr(ebox,fasta[allLoc],ignore.case = TRUE)
        sapply(mlo ,shift,width)
    }
    findMin<-function(x) motifs[which.min(x)]
    regionsA <-do.call(cbind,lapply(motifs,findEbox,width))
    regionsB <-do.call(cbind,lapply(compl,findEbox,width))
    combined<-abind(regionsA,regionsB,along=3)
    low3D<-Vectorize(function(i,j)min(combined[i,j,]))
    regions<-outer(seq(dim(combined)[1]),seq(dim(combined)[2]),low3D)
    distances<-apply(regions,1,findMin)
    lapply(motifs,function(motif,close){allLoc[motif==close]}, distances)
}

colnameUR<-function(data,fun,...) addColnames(data,fun(colnames(data),...))


# using only 1 and 4
swapS20<-"PC1-%pcv% jurk PC1+%pcv% eryt PC2+%pcv% NA PC2-%pcv% NA PC4+%pcv% cd34 PC4-%pcv% ecfc Leukemia jurk ECFC ecfc HSC cd34 Erythroid eryt "

swapS5<-"PC1-%pcv% jurk PC1+%pcv% eryt PC3+%pcv% ecfc PC3-%pcv% NA HSC NA NotHSC NA Other cd34 NotOther NA Leukemia jurk ECFC ecfc HSC cd34 Erythroid eryt "

convPvar<-function(swap)function(x)simpleSwapFun(gsub("%pcv%",x,swap))
conv20s3<-convPvar(swapS20)(3)
conv5s3<-convPvar(swapS5)(3)
conv20s1<-convPvar(swapS20)(1)
conv5s1<-convPvar(swapS5)(1)
conv20sp5<-convPvar(swapS20)(0.5)
conv5sp5<-convPvar(swapS5)(0.5)


#conv5pr<-simpleSwapFun("PC1-3 PC1-3 PC1+3 PC1+3 PC3+3 ecfc PC3-3 NA HSC NA NotHSC NA Other cd34 NotOther NA")

#selection<-function(x,conv){    
#    x[,conv(colnames(x))!="NA"]    
#}
selection<-function(x,conv)
    applySwapFun(x,Compose(conv,swapFunD))[,unique(swapFunD(categories))]

applySel<-function(x,conv){
    newColNames<-swapFunB(conv(colnames(x)))
    colnames(x)<-newColNames
    x[,Filter(function(x) x!="NA",newColNames)]
    
}
    


mfn16Make<-function(pvalue=5,heights="udm",control="combined"){
    paste("~/Dropbox/UTX-Alex/br-data/",heights,"/",heights,"_pvalue=",pvalue,"_",control,"_mock.txt",sep="")
}


genReg5pCombined<-function(reg5,sd){
    HSC5<-addColnames(cbind(reg5[,paste0("PC3+",sd)]&reg5[,paste0("PC7+",sd)],
                            reg5[,paste0("PC3-",sd)]&reg5[,paste0("PC7-",sd)])|(reg5[,paste0("PC3-",sd)]&reg5[,paste0("PC7+",sd)])|(reg5[,paste0("PC3+",sd)]&reg5[,paste0("PC7-",sd)]),c("HSC","NotHSC"))
    
    Other5<-addColnames(cbind(reg5[,paste0("PC3+",sd)]&reg5[,paste0("PC5+",sd)],
                              reg5[,paste0("PC3-",sd)]&reg5[,paste0("PC5-",sd)])|(reg5[,paste0("PC3-",sd)]&reg5[,paste0("PC5+",sd)])|(reg5[,paste0("PC3+",sd)]&reg5[,paste0("PC5-",sd)]),c("Other","NotOther"))
    cbind(reg5[,1:4],HSC5,Other5)
}


genReg5pCombined3sd<-function(reg5){
    HSC5<-addColnames(cbind(reg5[,"PC3+3"]&reg5[,"PC7+3"],
                            reg5[,"PC3-3"]&reg5[,"PC7-3"])|(reg5[,"PC3-3"]&reg5[,"PC7+3"])|(reg5[,"PC3+3"]&reg5[,"PC7-3"]),c("HSC","NotHSC"))
    
    Other5<-addColnames(cbind(reg5[,"PC3+3"]&reg5[,"PC5+3"],
                              reg5[,"PC3-3"]&reg5[,"PC5-3"])|(reg5[,"PC3-3"]&reg5[,"PC5+3"])|(reg5[,"PC3+3"]&reg5[,"PC5-3"]),c("Other","NotOther"))
    cbind(reg5[,1:4],HSC5,Other5)
}



genReg5pCombined1sd<-function(reg5){
    HSC5<-addColnames(cbind(reg5[,"PC3+1"]&reg5[,"PC7+1"],
                            reg5[,"PC3-1"]&reg5[,"PC7-1"])|(reg5[,"PC3-1"]&reg5[,"PC7+1"])|(reg5[,"PC3+1"]&reg5[,"PC7-1"]),c("HSC","NotHSC"))
    
    Other5<-addColnames(cbind(reg5[,"PC3+1"]&reg5[,"PC5+1"],
                              reg5[,"PC3-1"]&reg5[,"PC5-1"])|(reg5[,"PC3-1"]&reg5[,"PC5+1"])|(reg5[,"PC3+1"]&reg5[,"PC5-1"]),c("Other","NotOther"))
    cbind(reg5[,1:4],HSC5,Other5)
}



geneMatrix<-function(over,reg,regions,geneList,id="name2"){
    geneCount<-function(y,allGenes){
        d<-setdiff(allGenes,y)
        t<-rbind(data.frame(names=y,count=rep(1,length(y))),data.frame(names=d,count=rep(0,length(d))))
        t[order(as.character(t[,1])),2]
    }
    
    genes<-lapply(seq(dim(reg)[2]),function(x) greatGeneAssoc(over[reg[,x],c(1,2,3)],regions,geneList)[,id])
    allGenes<-sort(unique(as.character(unlist(genes))))
    geneCountWrapper<-Compose(as.character,unique,function(x)geneCount(x,allGenes))
    geneMatrix<-addNames(do.call(cbind,lapply(genes,geneCountWrapper)),colnames(reg),allGenes)
    return(geneMatrix)
}


writegenes<-function(name,frame,header,conv){
    frameName<-function(name, frame){
        as.character(rownames(frame)[which(frame[,name]==1)])
    }    
    fileName<-paste(header,conv(name),".table",sep="")
    write.table(frameName(name,frame),fileName,quote=FALSE,col.names=FALSE,row.names = FALSE)
}

writegenesNoConv<-function(name,frame,header,conv){
    frameName<-function(name, frame){
        as.character(rownames(frame)[which(frame[,name]==1)])
    }    
    fileName<-paste(header,conv(name),".table",sep="")
    write.table(frameName(name,frame),fileName,quote=FALSE,col.names=FALSE,row.names = FALSE)
}
