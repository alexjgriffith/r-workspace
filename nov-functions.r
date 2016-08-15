library(CCCA)
library(functional)
library(parallel)
library(ggplot2)
library(grid)
library(abind)

getCor<-function(i){
    file<-paste("~/thesis-november/22x22-pvalue=",i,"_heights.matrix",sep="")
    test<-read.table(file,header=T)
    colnames(test)<-categories
    cor(test)
}


readTagMatrix<-function(i,cats=NULL,fun=function(x) x){
    getOve(i,cats=cats,fun,loc="",tail="")
}

getOve<-function(i,cats=NULL,fun=overlap,loc="~/thesis-november/4x4-pvalue=",categories=NULL,tail=".matrix"){
    file<-paste(loc,i,tail,sep="")
    test<-read.table(file,header=T)
    test<-test[4:length(test)]
    if(! is.null(categories))
        test<-test[,categories]
    depth<-dim(test)[1]
    if(is.null(cats))
        cats<-colnames(test)
    out<-matrix(numeric(depth*length(cats)),ncol=length(cats))
    colnames(out)<-cats
    for(i in cats){
        if(i %in% colnames(test))
            out[,i]=test[,i]
    }
    fun(out)
}


overlap<-function (f) 
{
    permutations<-function(n)
        do.call(rbind,lapply(seq(n),function(x) cbind(seq(n),x)))

    overlapNorm <- function(x, f) {
        a <- f[,x[1]]
        b <- f[,x[2]]
        nu<-sum(a & b)
        de<-sum(a | b)
        if(de==0)
            return (1)
        else
            return (nu/de)
    }    
    width <- dim(f)[2]
    m <- matrix(apply( apply(permutations(seq(width)),1, as.numeric), 2, overlapNorm, f ), ncol = width)
    colnames(m) <- colnames(f)
    rownames(m) <- colnames(f)
    m
}


getPvalue<-function(n=30){    
    pvalues<-seq(from=0,by=2.5,length=n)
    pvalues
}

createZip<-function(x)
    Map(function(i,j)cbind(x[i],x[j]),seq(1,length(x),2),seq(2,length(x),2))
getSwapCats<-function(cats,names){    
    swCat<-cats
    names(swCat)<-names
    function(x){as.character(unlist(lapply(x,function(x) swCat[x])))}
}
stringToSwap<-function(x)do.call(getSwapCats,splitZip(createZip(strsplit(x," ")[[1]])))
splitZip<-function(inList)
    list(sapply(inList,"[",1),sapply(inList,"[",2))
getSwapCats<-function(cats,names){    
    swCat<-cats
    names(swCat)<-names
    function(x){as.character(unlist(lapply(x,function(x) swCat[x])))}
}
simpleSwapFun<-function(char){
stringToSwap(paste (rev(strsplit(char, " ")[[1]]),collapse=" "))
}

swapFunB<-stringToSwap(paste (rev(strsplit("rpmi Leukemia tall_p3 Leukemia tall_p2 Leukemia tall_p1 Leukemia tall_p2_1 Leukemia tall_p2_2 Leukemia tall_p3_1 Leukemia tall_p3_2 Leukemia jurk_sandar_1 Leukemia jurk_sandar Leukemia jurk Leukemia jurk_1 Leukemia jurk_2 Leukemia rpmi_1 Leukemia rpmi_2 Leukemia cem_1 Leukemia cem_2 Leukemia cem Leukemia ecfc-tsa ECFC ecfc_old ECFC ecfc.tsa ECFC ecfc ECFC meka MEKA cd133 HSC cd34 HSC cd34_new HSC eryt Erythroid eryt_1 Erythroid eryt_2 Erythroid eryt_f Erythroid eryt_a Erythroid k562 Erythroid k562_1 Erythroid k562_2 Erythroid", " ")[[1]]),collapse=" "))

swapFun<-stringToSwap(paste (rev(strsplit("rpmi RPMI tall_p3 Prima5 tall_p2 Prima2 tall_p1 Prima5 tall_p2_1 Prima2 tall_p2_2 Prima2 tall_p3_1 Prima5 tall_p3_2 Prima5 jurk_sandar_1 Jurkat jurk_sandar Jurkat jurk Jurkat jurk_1 Jurkat jurk_2 Jurkat rpmi_1 RPMI rpmi_2 RPMI cem_1 CEM cem_2 CEM cem CEM ecfc-tsa ECFC ecfc_old ECFC ecfc.tsa ECFC ecfc ECFC meka MEKA cd133 HSC cd34 HSC cd34_new HSC eryt Erythroid eryt_1 Erythroid eryt_2 Erythroid eryt_f Erythroid eryt_a Erythroid k562 K562 k562_1 K562 k562_2 K562", " ")[[1]]),collapse=" "))

swapFunC<-stringToSwap(paste (rev(strsplit("Leukemia blue Erythroid red HSC orange ECFC green MEKA brown", " ")[[1]]),collapse=" "))

mockToNorm<-simpleSwapFun("cd133_mock cd133 cd34_mock cd34 cd34_new_mock cd34_new cem_mock cem_1 ecfc-tsa_mock ecfc-tsa eryt_a_mock eryt_a eryt_f_mock eryt_f jurk_mock jurk jurk_sandar_mock jurk_sandar k562_mock_1 k562_1 k562_mock_2 k562_2 meka_mock meka rpmi_mock_1 rpmi_1 tall_p1_mock tall_p1 tall_p2_mock tall_p2_1 tall_p3_mock tall_p2_1")

## combines Meka and HSC
swapFunD<-stringToSwap(paste (rev(strsplit("rpmi Leukemia tall_p3 Leukemia tall_p2 Leukemia tall_p1 Leukemia tall_p2_1 Leukemia tall_p2_2 Leukemia tall_p3_1 Leukemia tall_p3_2 Leukemia jurk_sandar_1 Leukemia jurk_sandar Leukemia jurk Leukemia rpmi_1 Leukemia rpmi_2 Leukemia cem_1 Leukemia cem_2 Leukemia cem Leukemia ecfc-tsa ECFC ecfc_old ECFC ecfc.tsa ECFC ecfc ECFC meka HSC cd133 HSC cd34 HSC cd34_new HSC eryt Erythroid eryt_f Erythroid eryt_a Erythroid k562 Erythroid k562_1 Erythroid k562_2 Erythroid", " ")[[1]]),collapse=" "))


#sdf<-shiftdf("ecfc-tsa",0,-0.05)
#sdf<-shiftdf("cd34",-0.12,0,sdf)
#df<-data.frame(x=matr[,"PC1",1],y=matr[,"PC3",1],categories=swapFun(categories),Conditions=cols)
#postext<-shiftCols(df$x,df$y,categories,sdf)
#p<-ggplot(df,aes(x=x,y=y,col=Conditions,label=categories))+ylab("PC3")+xlab("PC1")+geom_point(size=3,shape=1)+theme_bw()+geom_text(x=postext$x,y=postext$y,show_guide=F,size=4)+scale_color_manual(values=c("green","red" , "orange","blue","brown"))
shiftdf<-function(category,x,y ,df=NULL){
    dft<-data.frame(categories=category, x=x, y=y)
    if(is.null(df))
        return(dft)
    else
        return(rbind(df,dft))
}
shiftCols<-function(x,y,categories,shiftdf){
    nx<-diff(range(x))/2
    ny<-diff(range(y))/2
    cats<-as.character(unlist(categories))
    cat<-as.character(unlist(shiftdf$categories))
    for(i in seq(dim(shiftdf)[1])){
        if(cat[i] %in% cats){
            print(shiftdf$y*ny)
           x[which(cat[i]==cats)]= x[which(cat[i]==cats)]+shiftdf$x[i]*nx
           y[which(cat[i]==cats)]= y[which(cat[i]==cats)]+shiftdf$y[i]*ny
       }
    }
    data.frame(x=x,y=y,categories=cats)
}

applySwapFun<-function(data,swapFun){
    colnames(data)<-swapFun(colnames(data))
    data
}


readCategories<-Compose(read.table,unlist,as.character)

contFile<-"~/Dropbox/UTX-Alex/jan/contcats"
conts<-readCategories(contFile)
conts<-conts[conts!="ecfc_old_mock"]

seql<-Compose(length,seq)
catFile<-"~/Dropbox/UTX-Alex/jan/catagories"
categories<-readCategories(catFile)
categoriesB<-c("cd34_new","eryt","jurk","cem_1")
#AFS_22x22_0<-makeAFS(makePeakFiles(categories),categories)
makePeakFiles<-function(categories,mock="combined")paste("~/Dropbox/Data/august_peaks/",categories,"~",mock,"_mock_peaks.xls",sep="")

wrapAFS<-function(makePeakFiles,categories,mock="combined"){
    function(pvalue=0)
        makeAFS(makePeakFiles(categories,mock),categories,pValue=pvalue)
}
renameColumn<-function(x,colName,newColName)
    {
        UseMethod("renameColumn",x)
    }
renameColumn.data.default<-function(x,colName,newColName){
    warning("renameColumn only works for data.frame currently.")
    x
}
renameColumn.data.frame<-function(x,colName=NULL,newColName=NULL){
    if(! is.null(colName) | ! is.null(newColName))
    whichName<-(colnames(x)==colName)
    if(! any(whichName))
        warning(colName," was not found in colnames(x)!")
    else{
    newCols<-replace(colnames(x),whichName,newColName)
    colnames(x)<-newCols}
    x
}

nameColumns<-function(x,names){
    colnames(x)<-names
    x
}
    

#rawDataFiles<-function(categories)paste("/mnt/brand01-00/mbrand_analysis/data_sets/",categories,"/",categories,"_unique_nodupes.bed",sep="")
#makePeakFiles<-function(categories)paste("~/.peaktemp/",categories,"~combined_mock_peaks.xls",sep="")
#makePeakFiles<-function(categories)paste("~/Dropbox/Data/august_peaks/",categories,"~combined_mock_peaks.xls",sep="")
all(sapply(makePeakFiles(categories),file.exists))


dend<-function(d,dist=function(x) as.dist(1-cor(x)),norm=qn,linkage=NULL,distanceName=NULL,colours=c("green","red","blue","orange"),n=4,swapFun=function(x) x){

    name=paste("Dendrogram ",
          ifelse(! is.null(distanceName),paste("Distance=",distanceName," ",sep=""),paste("Distance=","dissimilarity ",sep="")),
          ifelse(! is.null(linkage),paste("Linkage=",linkage,sep=""),
                 paste("Linkage=","complete ",sep=""))
          ,sep="")


    source("~/Masters/CCCA/inst/scipts/A2R.r")
    if(is.null(linkage))
        hc<-hclust(dist(norm(applySwapFun(d,swapFun))))
    else
        hc<-hclust(dist(norm(applySwapFun(d,swapFun))),linkage)
    if(n>1)
        A2Rplot(hc,k=n,col.down=colours,boxes=FALSE,main=name)
    else
        plot(hc,hang=-1,sub="",xlab="Cell Types",ylab=NULL,)
}



mergeFun<-function(ma,swapFun=swapFunB){
    newCols<-swapFun(colnames(ma))
    #print(newCols)
    unc<-unique(newCols)
    outL<-c()
    for(i in unc){
        pos<-which(i==newCols)
        #print(pos)
        #print(data.frame(c=colnames(ma)[pos],v=do.call(rbind,lapply(pos,function(x) sum(ma[,x])))))
        if(length(pos)>1)
            temp<-orM(ma[,pos])
        else{
            temp<-ma[,pos]
            #print(sum(temp))
        }
        outL<-cbind(outL,temp)}
    
    colnames(outL)<-unc
    outL
}


plotPCMat2D<-function(matr,pcs,categories,swapFun,swapFunB,swapFunC,...){
        df<-data.frame(x=matr[,pcs[1]],y=matr[,pcs[2]],categories=swapFun(categories),Conditions=swapFunB(categories))
        plotPCMatAux(df,pcs,categories,swapFunC(unique(sort(swapFunB(categories)))),...)
    }


plotPCMatAux<-function(df,pcs,categories,colours,sdf=NULL,text=NULL,legend=NULL){
    if(is.null(sdf))
       postext<-df
   else
       postext<-shiftCols(df$x,df$y,categories,sdf)
    p<-ggplot(df,aes(x=x,y=y,col=Conditions,label=categories))+ylab(pcs[2])+xlab(pcs[1])+geom_point(size=6,shape=20)+theme_bw()
    if(! is.null(text))
        p<-p+geom_text(x=postext$x,y=postext$y,show_guide=F,size=5)
    if(! is.null(colours))
        p<-p+scale_color_manual(values=colours)
    if(is.null(legend))
        p<-p+theme(legend.position="none")+scale_x_continuous(breaks=NULL)+scale_y_continuous(breaks=NULL)
    p
}

plotPCMat<-function(matr,pcs,categories,swapFun,swapFunB,n=1,sdf=NULL,colours=c("green","red" , "orange","blue","brown"),text=NULL,legend=NULL){    
    df<-data.frame(x=matr[,pcs[1],n],y=matr[,pcs[2],n],categories=swapFun(categories),Conditions=swapFunB(categories))
    if(is.null(sdf))
       postext<-df
   else
       postext<-shiftCols(df$x,df$y,categories,sdf)
    p<-ggplot(df,aes(x=x,y=y,col=Conditions,label=categories))+ylab(pcs[2])+xlab(pcs[1])+geom_point(size=4,shape=20)+theme_bw()
    if(! is.null(text))
        p<-p+geom_text(x=postext$x,y=postext$y,show_guide=F,size=4)
    if(! is.null(colours))
        p<-p+scale_color_manual(values=colours)
    if(is.null(legend))
        p<-p+theme(legend.position="none")+scale_x_continuous(breaks=NULL)+scale_y_continuous(breaks=NULL)
    p
}

pass<-function(x) x

mfnMake<-function(fun)
    function(pvalue=0,heights=FALSE){    
        if(heights)
            heights="_heights"
        else
            heights=""
        fun(pvalue,heights)
}

mfn4<-mfnMake(function(pvalue,heights)
    paste("~/thesis-november/4x4-pvalue=",pvalue,heights,".matrix",sep=""))

mfn22<-mfnMake(function(pvalue,heights)
    paste("~/thesis-november/22x22-pvalue=",pvalue,heights,".matrix",sep=""))

mfn6<-mfnMake(function(pvalue,heights)
    paste("~/thesis-november/6x6-pvalue=",pvalue,heights,".matrix",sep=""))

mfnC<-mfnMake(function(pvalue,heights)
    paste("~/thesis-november/cont-pvalue=",pvalue,heights,".matrix",sep=""))


mergeFun<-function(ma,swapFun){
    newCols<-swapFun(colnames(ma))
    #print(newCols)
    unc<-unique(newCols)
    outL<-c()
    for(i in unc){
        pos<-which(i==newCols)
        #print(pos)
        #print(data.frame(c=colnames(ma)[pos],v=do.call(rbind,lapply(pos,function(x) sum(ma[,x])))))
        if(length(pos)>1)
            temp<-orM(ma[,pos])
        else{
            temp<-ma[,pos]
            print(sum(temp))
        }
        outL<-cbind(outL,temp)}
    
    colnames(outL)<-unc
    outL
}

pvalues<-seq(from=0,by=2.5,length=30)

multiPlot<-function(plots,rows,cols){
pushViewport(viewport(layout=grid.layout(rows,cols)))
for(i in seq(length(plots))){
    r<-modulous((i-1),cols)+1
    c<-floor((i-1)/cols)+1
    print(plots[[i]],vp=viewport(layout.pos.row=c,layout.pos.col=r))
    if(i>=rows*cols)
        break
}
}


truthTable<-function(n){
    count<-function(s,n) {
        rep(c(rep(0,n/(2^s)),rep(1,n/(2^s))),2^(s-1))
    }
    p<-2^n
    do.call(cbind,lapply(seq(n),count,p))
   
}

addColnames<-function(matrix,colnames){
    colnames(matrix)<-colnames
    matrix
}

addRownames<-function(matrix,rownames){
    rownames(matrix)<-rownames
    matrix
}

addNames<-function(matrix,colnames,rownames=colnames,list=NULL){
    if(is.null(list)){
        colnames(matrix)<-colnames
        rownames(matrix)<-rownames
    }
    else{
        names(matrix)<-colnames
    }
    matrix        
}


stackedContribWrapper<-function(heightFile,tagFile,cats,pcs,swapFunB,swapFunC,sum=TRUE){
    source("~/Dropbox/thesis/R/contrib_new.R")
    tags<-renameColumn(read.table(tagFile,header=T),"ecfc.tsa","ecfc-tsa")
    
    tags<-tags[,4:dim(tags)[2]]    
    tags<-(function(){
        n<-1
        newMatrix<-do.call(cbind,lapply(seq(length(categories)),function(i){ if(categories[i] %in% colnames(tags)) {r<-tags[,n]; n<<-n+1;r} else rep(0,dim(tags)[1])}))
        colnames(newMatrix)<-cats
        newMatrix
    })()
    heights<-addColnames(read.table(heightFile,header=T),categories)
    prc<-pca(heights)
    t2<-mergeFun(tags) 
    lapply(pcs,function(pc)        
        stackedContrib(prc$eigenVectors[,pc[1]],"contrib2",t2,swapFun=swapFunB,colourOveride=swapFunC,sum=sum,blank=TRUE)
 )
}

readMatrixData<-function(filename,header=TRUE,categories=NULL){
    dA<-read.table(filename,header=header)
    list(bed=dA[,1:3],data=dA[,4:dim(dA)[2]])
}

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

pcsO<-paste("PC",seq(16),sep="")


pca2Matr<-function(x,n=pass){
    t(x$normData)%*%apply(x$eigenVectors,2,n)
}



removeColumn<-function(df,c){
    if(is.vector(df))
        removeColumn.vector(df,c)
    else if (is.matrix(df) & ! is.numeric(c))
        removeColumn.data.frame(df,c)
    else
        UseMethod("removeColumn",df)    
}

removeColumn.vector<-function(v,c){
    Filter(function(x) ! x %in% c,v)
}

removeColumn.data.frame<-function(df,c){
    df[,removeColumn(colnames(df),c)]
}


X2<-function(n,m,rot=TRUE){
    a<-rep(seq(m),n)
    b<-unlist(lapply(seq(n),rep,m))
    if(rot)
        c(rbind(a,b))
    else
        c(rbind(b,a))
}

getRange<-function(list,limit){
    range(sort(list)[floor(length(list)*(1-limit)+1):ceiling(length(list)*limit)])    
}
    

heatmapWrapper<-function(x,range=range(x)){
    to<-hclust(dist(x),"single")$order
    print(categories[to])
heatmap(x[rev(to),to],hclustfun=function(x) hclust(x,"single"),keep.dendro=FALSE,scale="none",zlim=range,Rowv=NA,Colv=TRUE,cexRow=1.5,cexCol=1.5,margin=c(6,6))    
}

MSD<-function(m1,m2){
    mean((dist(m1)-dist(m2))^2)
}

limitRange<-function(v)
    (v-range(v)[1])/diff(range(v))
