pvalues<-seq(from=5,to=75,by=2.5)
over<-readAFS(pvalues[9],"single")
bed<-orderBed(over)
heights<-readUDM(pvalues[9],"treatment","single")
prc=pca(heights)
