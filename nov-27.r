source("~/r-workspace/nov-functions.r")

readMatrixData<-function(filename,header=TRUE){
    dA<-read.table(filename,header=header)
    list(bed=dA[,1:3],data=dA[,4:dim(dA)[2]])
}

opts<-apply(truthTable(4),1,function(x) which(as.logical(x)))

peakData<-readMatrixData(mfn4(pvalues[1]))$data
pops<-opts[which(sapply(opts ,function(x) length(x))>0)]


comb<-unlist(lapply(pops,function(x) paste(categoriesB[x],collapse = "_")))
ret<-unlist(lapply(pops,function(x) length(which(apply(cbind(0,peakData[,x]),1,sum)==length(x))) ))


## overlaps for venn diagram
allIn<-function(x,y){
    which(andM(t(apply(as.matrix(y),1,"==",1))[,which(as.logical(x))]))    
}

df<-data.frame(name=comb,count=ret,depth=unlist(lapply(pops,length)),nameColumns(truthTable(4)[2:16,],categoriesB))

tb<-sapply(seq(1,15), function(i){
    toS<-allIn(df[i,4:7], df[,4:7])
    sum(df$count[toS]*(-1)^(df$depth[toS]-df$depth[i])) 
}
)
df$count2=tb




paste(overlap(peakData),sep=" & ")

heightData<-read.table(mfn22(pvalues[1],heights=TRUE),header=TRUE)
colnames(heightData)<-categories


data<-pca(heightData)

plotPCs(data$eigenVectors,c(1,3),data$normData,categories)

dta1<-cor(heightData[,categories])
dta<-round(dta1,2)
dta1<-applySwapFun(dta1,swapFun)
rownames(dta)<-colnames(dta1)
colnames(dta)<-colnames(dta1)
cat(makeLatexTable(dta))

makeLatexTable<-function(x){
   
    sidecall<-rownames(x)
    dimen<-dim(x)
    #x<-as.character(unlist(x))
    top<-paste("\\begin{table} \n\\begin{center}\n\\begin{tabular}{|c|",do.call(paste, as.list(c(rep("c",dimen[2]),sep=""))),"|}\n \\hline \n",sep="")
    header<-do.call(paste, as.list(c("\t& ",paste(colnames(x)[1:(dimen[2]-1)]," & ",sep=""),colnames(x)[dimen[2]], "  \\\\\n\\hline",sep="")))
    tail<-"\\end{tabular}\n\\end{center}\n\\end{table}\n"
    col<-paste(top,header,sep="")
    for (i in seq(dimen[1])){
        col<-paste(col,rownames(x)[i]," & ")
        for (j in seq(dimen[2])){
            if(j!=dimen[2])
                col<-paste(col," ",x[i,j]," & ",sep="")
            else
                col<-paste(col,x[i,j]," \\\\ \n",sep="")
        }}
    col<-paste(col," \\hline \n",tail,sep="")
    col
}

