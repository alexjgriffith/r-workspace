source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/jan-functions.r")
source("~/r-workspace/feb-functions.r")

env<-getPRC20();

sdf<-shiftdf("ecfc-tsa",0,-0.05)
sdf<-shiftdf("cd34",-0.15,0,sdf)
sdf<-shiftdf("cd133",-0.15,0,sdf)
sdf<-shiftdf("cd34_new",0,+.05,sdf)
sdf<-shiftdf("meka",0,-.05,sdf)
sdf<-shiftdf("eryt_f",-.2,0,sdf)
sdf<-shiftdf("eryt",0,.05,sdf)
sdf<-shiftdf("eryt_a",-.2,0,sdf)
sdf<-shiftdf("k562_1",-.15,0,sdf)
sdf<-shiftdf("k562_2",-.15,.05,sdf)
sdf<-shiftdf("jurk_sandar_1",-0.07,0.05,sdf)
sdf<-shiftdf("rpmi_1",.15,-.01,sdf)
sdf<-shiftdf("rpmi_2",.17,0.01,sdf)
sdf<-shiftdf("tall_p1",.01,+.05,sdf)
sdf<-shiftdf("tall_p2_1",.1,0,sdf)
sdf<-shiftdf("tall_p2_2",.03,-.07,sdf)
sdf<-shiftdf("jurk",0.15,0,sdf)
sdf<-shiftdf("tall_p3_1",-0.1,+0.05,sdf)
sdf<-shiftdf("tall_p3_2",0.15,+0.05,sdf)
sdf<-shiftdf("cem_1",-0.11,.05,sdf)
sdf<-shiftdf("cem_2",.15,.05,sdf)
sdf<-shiftdf("jurk_sandar",0,-.05,sdf)
cols<-swapFunB(categories)
#df<-data.frame(x=pca2Matr(env$prc,normalize)[,"PC1"],y=pca2Matr(env$prc,normalize)[,"PC4"],categories=swapFunM(env$categories),Conditions=cols)
#postext<-shiftCols(df$x,df$y,categories,sdf)
#a<-ggplot(df,aes(x=x,y=y,col=Conditions,label=categories))+ylab("PC3")+xlab("PC1")+geom_point(size=6,shape=20)+theme_bw()+geom_text(x=postext$x,y=postext$y,show.legend = F,size=4)+scale_color_manual(values=c("green","red" , "orange","blue","brown"))+theme(legend.position="none")

a<-plotPCMat2D(pca2Matr( env$prc) ,c("PC1","PC4"), categories=categories, sdf=sdf,swapFun=swapFun,swapFunB=swapFunB,swapFunC=swapFunC,text=TRUE)
                                        #scale_color_manual(values=c("green","red" , "orange","blue","brown"))

ggsave("~/Dropbox/Data/new-figs/pca_1-4_pvalue_20_sd_1_labled-a.png",a)
a<-ggplot(df,aes(x=x,y=y,col=Conditions,label=categories))+ylab("PC3")+xlab("PC1")+geom_point(size=6,shape=20)+theme_bw()+geom_text(x=postext$x,y=postext$y,show.legend = F,size=4)+scale_color_manual(values=c("green","red" , "orange","blue","orange"))+theme(legend.position="none")
                                        #scale_color_manual(values=c("green","red" , "orange","blue","brown"))
ggsave("~/Dropbox/Data/new-figs/pca_1-4_pvalue_20_sd_1_labled-b.png",a)



#############

env<-getPRC20(3);

sdf<-shiftdf("ecfc-tsa",0,-0.05)
sdf<-shiftdf("cd34",-0.15,0,sdf)
sdf<-shiftdf("cd133",-0.15,0,sdf)
sdf<-shiftdf("cd34_new",0,+.05,sdf)
sdf<-shiftdf("meka",0,-.05,sdf)
sdf<-shiftdf("eryt_f",0,-.05,sdf)
sdf<-shiftdf("eryt",-.2,0.05,sdf)
sdf<-shiftdf("eryt_a",-.2,-.03,sdf)
sdf<-shiftdf("k562_1",-.1,-.01,sdf)
sdf<-shiftdf("k562_2",-.1,.05,sdf)
sdf<-shiftdf("jurk_sandar_1",0.1,0,sdf)
sdf<-shiftdf("rpmi_1",-.1,-.05,sdf)
sdf<-shiftdf("rpmi_2",0.1,-0.05,sdf)
sdf<-shiftdf("tall_p1",-.03,-.05,sdf)
sdf<-shiftdf("tall_p2_1",0.03,-0.05,sdf)
sdf<-shiftdf("tall_p2_2",-.08,-.05,sdf)
sdf<-shiftdf("jurk",0.1,0.05,sdf)
sdf<-shiftdf("tall_p3_1",-0.1,0,sdf)
sdf<-shiftdf("tall_p3_2",0.13,0,sdf)
sdf<-shiftdf("cem_1",+0.07,.05,sdf)
sdf<-shiftdf("cem_2",-.07,0.04,sdf)
sdf<-shiftdf("jurk_sandar",0,+.05,sdf)


#a<-plotPCMat2D(pca2Matr( env$prc) ,c("PC1","PC2"), categories=categories, sdf=sdf,swapFun=swapFun,swapFunB=swapFunM,swapFunC=swapFunC,text=TRUE)

cols<-swapFunB(categories)
df<-data.frame(x=pca2Matr(env$prc,normalize)[,"PC1"],y=pca2Matr(env$prc,normalize)[,"PC2"],categories=swapFunM(env$categories),Conditions=cols)
postext<-shiftCols(df$x,df$y,categories,sdf)
a<-ggplot(df,aes(x=x,y=y,col=Conditions,label=categories))+ylab("PC3")+xlab("PC1")+geom_point(size=4,shape=20)+geom_text(x=postext$x,y=postext$y,show.legend = F,size=4)+scale_color_manual(values=c("green","red" , "orange","blue","brown"))+theme(legend.position="none")+theme_bw()+theme(legend.position="none")+scale_x_continuous(breaks=NULL)+scale_y_continuous(breaks=NULL)

                                        #scale_color_manual(values=c("green","red" , "orange","blue","brown"))

ggsave("~/Dropbox/Data/new-figs/pca_1-2_pvalue_20_labled-a.png",a)

#a<-ggplot(df,aes(x=x,y=y,col=Conditions,label=categories))+ylab("PC3")+xlab("PC1")+geom_point(size=6,shape=20)+theme_bw()+geom_text(x=postext$x,y=postext$y,show.legend = F,size=4)+scale_color_manual(values=c("green","red" , "orange","blue","orange"))+theme(legend.position="none")
                                        #scale_color_manual(values=c("green","red" , "orange","blue","brown"))

ggsave("~/Dropbox/Data/new-figs/pca_1-2_pvalue_20_labled-b.png",a)


###########

env<-getPRC5();

sdf<-shiftdf("ecfc-tsa",0,-0.05)
sdf<-shiftdf("cd34",-0.15,0,sdf)
sdf<-shiftdf("cd133",-0.15,0,sdf)
sdf<-shiftdf("cd34_new",0,+.05,sdf)
sdf<-shiftdf("meka",0,-.05,sdf)
sdf<-shiftdf("eryt_f",0,0,sdf)
sdf<-shiftdf("eryt",0,0,sdf)
sdf<-shiftdf("eryt_a",0,0,sdf)
sdf<-shiftdf("k562_1",0,0,sdf)
sdf<-shiftdf("k562_2",0,0,sdf)
sdf<-shiftdf("jurk_sandar_1",0,0,sdf)
sdf<-shiftdf("rpmi_1",-.05,.5,sdf)
sdf<-shiftdf("rpmi_2",.5,.5,sdf)
sdf<-shiftdf("tall_p1",0,0,sdf)
sdf<-shiftdf("tall_p2_1",0,.0,sdf)
sdf<-shiftdf("tall_p2_2",0,0,sdf)
sdf<-shiftdf("jurk",0,0,sdf)
sdf<-shiftdf("tall_p3_1",0,0,sdf)
sdf<-shiftdf("tall_p3_2",0,0,sdf)
sdf<-shiftdf("cem_1",+0,0,sdf)
sdf<-shiftdf("cem_2",0,0,sdf)
sdf<-shiftdf("jurk_sandar",0,0,sdf)
cols<-swapFunB(categories)
df<-data.frame(x=pca2Matr(env$prc,normalize)[,"PC1"],y=pca2Matr(env$prc,normalize)[,"PC3"],categories=swapFunM(env$categories),Conditions=cols)
postext<-shiftCols(df$x,df$y,categories,sdf)
ggplot(df,aes(x=x,y=y,col=Conditions,label=categories))+ylab("PC3")+xlab("PC1")+geom_point(size=6,shape=20)+theme_bw()+geom_text(x=postext$x,y=postext$y,show.legend = F,size=4)+scale_color_manual(values=c("green","red" , "orange","blue","brown"))+theme(legend.position="none")
                                        #scale_color_manual(values=c("green","red" , "orange","blue","brown"))


ggsave("~/Dropbox/Data/new-figs/pca_1-2_pvalue_20_labled-a.png",a)
a<-ggplot(df,aes(x=x,y=y,col=Conditions,label=categories))+ylab("PC3")+xlab("PC1")+geom_point(size=6,shape=20)+theme_bw()+geom_text(x=postext$x,y=postext$y,show.legend = F,size=4)+scale_color_manual(values=c("green","red" , "orange","blue","orange"))+theme(legend.position="none")
                                        #scale_color_manual(values=c("green","red" , "orange","blue","brown"))
ggsave("~/Dropbox/Data/new-figs/pca_1-2_pvalue_20_labled-b.png",a)




sdf<-shiftdf("eryt_f",0,0,sdf)
sdf<-shiftdf("eryt",0,0,sdf)
sdf<-shiftdf("eryt_a",0,0,sdf)
sdf<-shiftdf("k562_1",0,0,sdf)
sdf<-shiftdf("k562_2",0,0,sdf)
sdf<-shiftdf("jurk_sandar_1",0,0,sdf)
sdf<-shiftdf("rpmi_1",0,0,sdf)
sdf<-shiftdf("rpmi_2",0,0,sdf)
sdf<-shiftdf("tall_p1",0,0,sdf)
sdf<-shiftdf("tall_p2_1",0,.0,sdf)
sdf<-shiftdf("tall_p2_2",0,0,sdf)
sdf<-shiftdf("jurk",0,0,sdf)
sdf<-shiftdf("tall_p3_1",0,0,sdf)
sdf<-shiftdf("tall_p3_2",0,0,sdf)
sdf<-shiftdf("cem_1",+0,0,sdf)
sdf<-shiftdf("cem_2",0,0,sdf)
sdf<-shiftdf("jurk_sandar",0,0,sdf)





#############


env<-getPRC20(2);

sdf<-shiftdf("ecfc-tsa",0,-0.05)
sdf<-shiftdf("cd34",-0.15,0,sdf)
sdf<-shiftdf("cd133",-0.15,0,sdf)
sdf<-shiftdf("cd34_new",0,+.05,sdf)
sdf<-shiftdf("meka",0,.05,sdf)
sdf<-shiftdf("eryt_f",0,.05,sdf)
sdf<-shiftdf("eryt",-.2,-0.05,sdf)
sdf<-shiftdf("eryt_a",-.2,.03,sdf)
sdf<-shiftdf("k562_1",-.1,.01,sdf)
sdf<-shiftdf("k562_2",-.1,-.05,sdf)
sdf<-shiftdf("jurk_sandar_1",0.1,0,sdf)
sdf<-shiftdf("rpmi_1",-.1,.05,sdf)
sdf<-shiftdf("rpmi_2",0.1,0.05,sdf)
sdf<-shiftdf("tall_p1",-.03,.05,sdf)
sdf<-shiftdf("tall_p2_1",0.03,0.05,sdf)
sdf<-shiftdf("tall_p2_2",-.08,.05,sdf)
sdf<-shiftdf("jurk",0.1,-0.05,sdf)
sdf<-shiftdf("tall_p3_1",-0.1,0,sdf)
sdf<-shiftdf("tall_p3_2",0.13,0,sdf)
sdf<-shiftdf("cem_1",+0.07,-.05,sdf)
sdf<-shiftdf("cem_2",-.07,-0.04,sdf)
sdf<-shiftdf("jurk_sandar",0,-.05,sdf)


#a<-plotPCMat2D(pca2Matr( env$prc) ,c("PC1","PC2"), categories=categories, sdf=sdf,swapFun=swapFun,swapFunB=swapFunM,swapFunC=swapFunC,text=TRUE)

cols<-swapFunB(categories)
df<-data.frame(x=pca2Matr(env$prc,normalize)[,"PC1"],y=-pca2Matr(env$prc,normalize)[,"PC2"],categories=swapFunM(env$categories),Conditions=cols)

df<-data.frame(x=pca2Matr(env$prc,normalize)[,"PC1"],y=-pca2Matr(env$prc,normalize)[,"PC2"],categories=swapFun(env$categories),Conditions=cols)

postext<-shiftCols(df$x,df$y,categories,sdf)
a<-ggplot(df,aes(x=x,y=y,col=Conditions,label=categories))+ylab("-PC2")+xlab("PC1")+geom_point(size=4,shape=20)+geom_text(x=postext$x,y=postext$y,show.legend = F,size=4)+scale_color_manual(values=c("green","red" , "orange","blue","brown"))+theme(legend.position="none")+theme_bw()+theme(legend.position="none")+scale_x_continuous(breaks=NULL)+scale_y_continuous(breaks=NULL)
a
                                        #scale_color_manual(values=c("green","red" , "orange","blue","brown"))

ggsave("~/Dropbox/Data/new-figs/reverse_pca_1-2_pvalue_20_sd_2_labled-b.png",a)
