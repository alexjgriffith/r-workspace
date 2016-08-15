# requires utils.r
library(ggplot2)

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


pca2Matr<-function(x,n=pass){
    t(x$normData)%*%apply(x$eigenVectors,2,n)
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
