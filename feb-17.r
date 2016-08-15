source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/jan-functions.r")
source("~/r-workspace/feb-functions.r")
source("~/r-workspace/feb-variables.r")


env=getPRC20()

menv<-mergeFun(env$over[,4:dim(env$over)[2]],swapFunD)

unique<-apply(menv,1,sum)==1
cols<-swapFunC(unlist(apply(menv[unique,],1,function(x) names(which(x==1)))))
df<-data.frame(x=normalize(env$prc$eigenVectors[unique,1]),y=normalize(env$prc$eigenVectors[unique,2]),c=cols)
a<-ggplot(df,aes(x=x,y=y,col=c))+geom_point(alpha=1)+scale_color_manual(values=levels(df$c))+theme(legend.position="none")+scale_x_continuous(breaks=NULL)+scale_y_continuous(breaks=NULL)+
                theme(panel.background=element_blank(),
                      panel.grid.major=element_blank(),
                      panel.grid.minor=element_blank(),
                      panel.border=element_blank(),
                      axis.title.y=element_blank(),
                      axis.title.x=element_blank())
ggsave("~/Dropbox/Data/new-figs/pca-dist-1-2-pvalue_20.png",a)
df<-data.frame(x=normalize(env$prc$eigenVectors[unique,1]),y=normalize(env$prc$eigenVectors[unique,4]),c=cols)
b<-ggplot(df,aes(x=x,y=y,col=c))+geom_point(alpha=1)+scale_color_manual(values=levels(df$c))+theme(legend.position="none")+scale_x_continuous(breaks=NULL)+scale_y_continuous(breaks=NULL)+
                theme(panel.background=element_blank(),
                      panel.grid.major=element_blank(),
                      panel.grid.minor=element_blank(),
                      panel.border=element_blank(),
                      axis.title.y=element_blank(),
                      axis.title.x=element_blank())
ggsave("~/Dropbox/Data/new-figs/pca-dist-1-4-pvalue_20.png",b)


menv<-mergeFun(env$over[,4:dim(env$over)[2]],swapFunD)
unique<-apply(menv,1,sum)==1
cols<-unlist(apply(menv,1,function(x) if(sum(x)==1) swapFunC(names(which(x==1))) else "black"))
df<-data.frame(x=normalize(env$prc$eigenVectors[,1]),y=normalize(env$prc$eigenVectors[,2]),c=cols)
c<-ggplot(df,aes(x=x,y=y,col=c))+geom_point(alpha=1)+scale_color_manual(values=levels(df$c))+theme(legend.position="none")+scale_x_continuous(breaks=NULL)+scale_y_continuous(breaks=NULL)+
                theme(panel.background=element_blank(),
                      panel.grid.major=element_blank(),
                      panel.grid.minor=element_blank(),
                      panel.border=element_blank(),
                      axis.title.y=element_blank(),
                      axis.title.x=element_blank())
ggsave("~/Dropbox/Data/new-figs/pca-dist-1-2-pvalue_20-alt.png",c)
df<-data.frame(x=normalize(env$prc$eigenVectors[,1]),y=normalize(env$prc$eigenVectors[,4]),c=cols)
d<-ggplot(df,aes(x=x,y=y,col=c))+geom_point(alpha=1)+scale_color_manual(values=levels(df$c))+theme(legend.position="none")+scale_x_continuous(breaks=NULL)+scale_y_continuous(breaks=NULL)+
                theme(panel.background=element_blank(),
                      panel.grid.major=element_blank(),
                      panel.grid.minor=element_blank(),
                      panel.border=element_blank(),
                      axis.title.y=element_blank(),
                      axis.title.x=element_blank())
ggsave("~/Dropbox/Data/new-figs/pca-dist-1-4-pvalue_20-alt.png",d)



env=getPRC5()
menv<-mergeFun(env$over[,4:dim(env$over)[2]],swapFunD)
unique<-apply(menv,1,sum)==1
cols<-swapFunC(unlist(apply(menv[unique,],1,function(x) names(which(x==1)))))
df<-data.frame(x=normalize(env$prc$eigenVectors[unique,1]),y=normalize(env$prc$eigenVectors[unique,3]),c=cols)
e<-ggplot(df,aes(x=x,y=y,col=c))+geom_point(alpha=1)+scale_color_manual(values=levels(df$c))+theme(legend.position="none")+scale_x_continuous(breaks=NULL)+scale_y_continuous(breaks=NULL)+
                theme(panel.background=element_blank(),
                      panel.grid.major=element_blank(),
                      panel.grid.minor=element_blank(),
                      panel.border=element_blank(),
                      axis.title.y=element_blank(),
                      axis.title.x=element_blank())
ggsave("~/Dropbox/Data/new-figs/pca-dist-1-3-pvalue_5.png",e)
menv<-mergeFun(env$over[,4:dim(env$over)[2]],swapFunD)
unique<-apply(menv,1,sum)==1
cols<-unlist(apply(menv,1,function(x) if(sum(x)==1) swapFunC(names(which(x==1))) else "black"))
df<-data.frame(x=normalize(env$prc$eigenVectors[,1]),y=normalize(env$prc$eigenVectors[,3]),c=cols)
f<-ggplot(df,aes(x=x,y=y,col=c))+geom_point(alpha=1)+scale_color_manual(values=levels(df$c))+theme(legend.position="none")+scale_x_continuous(breaks=NULL)+scale_y_continuous(breaks=NULL)+
                theme(panel.background=element_blank(),
                      panel.grid.major=element_blank(),
                      panel.grid.minor=element_blank(),
                      panel.border=element_blank(),
                      axis.title.y=element_blank(),
                      axis.title.x=element_blank())
ggsave("~/Dropbox/Data/new-figs/pca-dist-1-3-pvalue_5-alt.png",f)


