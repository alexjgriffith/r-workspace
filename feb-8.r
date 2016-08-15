source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/jan-functions.r")
source("~/r-workspace/feb-functions.r")



cs<-makeForkCluster(22,reince=1)

with(list(),{
env20<-getPRC20(2)
tr20<-with(env20,{
    data<-orderBed(bed)
    tr<-peakDensity(data,rawDataFiles(categories),20,350,TRUE,clust=cs)
    tr
})
env5<-getPRC5(2)
tr5<-with(env5,{
    data<-orderBed(bed)
    tr<-peakDensity(data,rawDataFiles(categories),20,350,TRUE,clust=cs)
    tr
})
save(tr5,tr20,file="~/Dropbox/heights.RData")
 })
     
stopCluster(cs)

