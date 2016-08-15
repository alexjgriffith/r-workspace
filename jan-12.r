#' Jan 12 2016
#' 
#' Thesis work Using CCCA. Specificly working on assesing the
#' biological results of analysis in the hematopoietic domain
#'
#' \enumerate{Motes:
#'   \item four catagories of interest <tall eryt ecfc stem>
#' }
#' \enumerate{Goals:
#'   \item  Check Ebox Files <to-do>
#'   \item  Generate "Combined Mock" bed files <to-do>
#'           a. Bed files with unique locations for each contributor
#'           b. Bed files seperated using pca sep 3
#'   \item  Find Motifs and Anotate <to-do>
#'           a. First make a unified background and save it as
#'              all in the bed file location and a file which
#'              does not apear in any pca direction
#'           b. Use homer for each partition at 5 and 20 mac as
#'              well as the unique bed files
#'   \item  Find genes for all 16 combinations <to-do>
#' }
#'
#' @author Alexander Griffith 
NULL

# load source files
source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/jan-functions.r")
source("~/Masters/CCCA/inst/scipts/eboxFrequency.r")
source("~/Masters/mulcal/newR/rGREAT.r")
source("~/Dropbox/R/makeLatexTable.R")

# generate variables
source("~/r-workspace/jan-variables.r")
source("~/r-workspace/fasta-data.r")

# the motifDistances function is quite slow since it relies on gregexpr
# interfaced with dnastring type
eboxD5<-motifDistances(eboxLoc5,fasta5)
eboxD20<-motifDistances(eboxLoc20,fasta20)

# Count how many locations are in the eboxLoc list for both a macs cut off
# of 5 and 20
allDataTotal<-data.frame(all5=do.call(rbind,lapply(eboxLoc5,length)),
                    near5=do.call(rbind,lapply(eboxD5,length)),
                    all20=do.call(rbind,lapply(eboxLoc20,length)),
                    near20=do.call(rbind,lapply(eboxD20,length)))


allData<-addRownames(do.call(cbind,mapply(function(size,type)round(allDataTotal[type]/size,4),  sapply(list(fasta5,fasta5,fasta20,fasta20),length),colnames(allDataTotal))),rownames(allDataTotal))



motifFrequency<-function(regI,mlocI,sel,subs=NULL,motifs=NULL,cols=NULL){
    getNC<-function(x)
        if(is.null(colnames(x))){names(x)}else{colnames(x) }
    if(is.null(motifs))
        motifs=names(mlocI)
    if(is.null(cols))
        cols=getNC(regI)
    motifAllFrequencyAlt<-motifFreqGen(function(){function(x,y){
        length(intersect(x,y))
    }})
    if (! is.list(regI))
        temp<-lapply(seq(dim(regI)[2]),function(x) which(regI[,x]))
    else
        temp<-regI
    if(is.null(subs))
        subs=union(unlist(temp),unlist(mlocI))
    reg<-lapply(temp,intersect,subs)
    mloc<-lapply(mlocI,intersect,subs)
    all<-sapply(mloc,length)   
    count<-motifAllFrequencyAlt(motifs,cols,mloc,reg)
    countM<-rbind(sapply(reg,length),0)
    colnames(countM)<-cols;
    list(matrix=applySel(count,sel)[,unique(swapFunD(categories))],
         colSize=applySel(countM,sel)[1,unique(swapFunD(categories))],
         rowSize=all)   
}


motF<-list(    
    bed20D=motifFrequency(conPeakLoc20,eboxD20,conv20,which(conUniqueLocs20)),
    bed20A=motifFrequency(conPeakLoc20,eboxLoc20,conv20,which(conUniqueLocs20)),
    bed5D=motifFrequency(conPeakLoc5,eboxD5,conv5,which(conUniqueLocs5)),    
    bed5A=motifFrequency(conPeakLoc5,eboxLoc5,conv5,which(conUniqueLocs5)),
    pca20D=motifFrequency(reg20,eboxD20,conv20),    
    pca20A=motifFrequency(reg20,eboxLoc20,conv20),
    pca5D=motifFrequency(reg5p,eboxD5,conv5),
    pca5A=motifFrequency(reg5p,eboxLoc5,conv5))

save(motF,file="~/thesis-january/ebox-overlap.RData")
#load("~/thesis-january/ebox-overlap.RData")


### make a snapshot of bed data used
sreg20<-apply(selection(reg20,conv20),2,which)
sreg5p<-apply(selection(reg5p,conv5),2,which)
# conPeakLoc5
# conPeakLoc20

brss<-"Leukemia tall Erythroid eryt ECFC ecfc HSC stem"
brOut<-simpleSwapFun(brss)
brIn<-simpleSwapFun(paste(pairSwitch(strsplit(brss," ")[[1]]),collapse =" "))

bedData<-lapply(list(bed20=conPeakLoc20,bed5=conPeakLoc5,pca20=sreg20,pca5=sreg5p), function(x)lapply(x,function(x) bed20[x,]))

#save(bedData,file="~/thesis-january/ind-bed.RData")

saveBedFiles(bedData)
bedFileString<-function(s,v,p,c="combined")
    paste("bed_",s,"_",brOut(v),"_pvalue=",p,"_",c,"_mock.bed",sep="")

saveBedFiles<-function(bedData){
    writebed<-function(data,file){
        write.table(data,file,quote=FALSE,row.names=FALSE,sep="\t")
    }
    for (s in c("pca","bed"))    
        for(p in c(20,5))
            for(v in unique(swapFunD(categories))){
                data<-bedData[[paste(s,p,sep="")]][[v]][1:3,]
                file<-paste("~/Dropbox/UTX-Alex/br-data/bed/",
                            bedFileString(s,v,p),sep="")
                writebed(data,file);
            }
}


### find an all and a none region for pca5 and pca20 subset that we care about



with(list(noneBed=bed5[setdiff(seq(length(fasta5)),sort(unique(unlist(apply(selection(reg5p,conv5),2,which))))),],
          fileName=paste("~/Dropbox/UTX-Alex/br-data/bed/",bedFileString("pca","none",5),sep="")),{
              write.table(noneBed,fileName,quote=FALSE,row.names=FALSE)
})


oldFunc<-function(x,y) x+y

newFunc<-Curry(lapply,FUN=sum)

newFunc(list(1,2,3))

#.b<-function(...)paste(...,col="",sep="")
#.b("~/Dropbox/UTX-Alex/br-data/bed/",c(bedFileString("pca","none",5),bedFileString("pca","none",5)))
# this is paste0!!



# test if jit compilaion allows infinite recursion

testFun<-cmpfun(function(x){
    testFun<-function(x) {if(x<1){return ("complete")}else{testFun(x-1)}}
    testFun(x)
})

tf<-function(y) lapply(seq(y),function(x) x+1)
tfc<-cmpfun(tf)

system.time(tf(10000000))

system.time(tfc(10000000))
