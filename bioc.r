#source("http://bioconductor.org/biocLite.R")
#pkgs <- rownames(installed.packages())
#biocLite()

#Aux Functions
collapse<-function(...){paste(...,sep="")}
lcollapse<-function(x){br<-"";for(i in x){br<-paste(br,i,sep="")};br}
collect<-function(x,fn,...){lcollapse(sapply(x,fn,...))}
modulous<-function(x,m)
    {t1<-floor(x/m)
     (x-t1*m)}

#Logo generators
require("seqLogo")

getPWM<-function(motif,data){
    sets<-data$poisit[which(data$poisit==motif),]
    sets<-sets[order(sets[2]),3:length(sets[,1])]
    apply(sets,1,as.numeric)} 

printLogos<-function(data,shortName){
    for(i in data$name){
        x<-getPWM(i,data)
        logos<-makePWM(as.data.frame(x,row.names=c("A","B","C","D")))
        #png(paste("logo_",shortName,"_",gsub(">","",i),".png",sep=""),width=240,height=120)
        seqLogo(logos,xaxis=FALSE,yaxis=FALSE,ic.scale=FALSE)
        #dev.off()
        print(i)}}

fileLocation<-"/mnt/brand01-00/mbrand_analysis/projects/august/data/motifpeaks/"
sets<-c("normal-abnormalStem_not_normal-abnormalStem.pwm","abnormal-normalStem_not_abnormal-normalStem.pwm")
shortNames<-c("unormal","uabnormal")
data<-loadData(collapse(fileLocation,sets[1]))


data<-loadData("test.data")
motifs<-c()
for(j in seq(length(sets))){
    data<-loadData(paste(fileLocation,sets[j],sep=""))
    #printLogos(data,shortNames[j])
    motifs<-cbind(motifs,c(shortNames[j],data$name))}
motifs<-as.data.frame(t(as.data.frame(t(motifs[2:dim(motifs)[1],]), row.names=motifs[1,])))

br<-c()
for (i in seq(dim(a)[1]))
    {
        print(i)
        for (j in shortNames)
            print(j)
            print(motifs[j][i,1])
            br<-c(br,motifs[j][i,1], htmlTable(matrix(runif(3,0,1),3,1)))
    }

        

#Html ggenerator 
buildAnotations<-function(...){
    value<-c(...)
    l=length(value)
    br=""
    if(modulous(l,2)==0)
        collect(seq(from=1,to=l,by=2),function(x){collapse(" ",value[x],'="',value[x+1],'"')})}

        
htmlTags<-function(tag,value=FALSE,anotations=FALSE){
    closeT<-"/>"
    an<-""
    if(FALSE != anotations[1])
        an<-lcollapse(anotations)
    if(FALSE != value)
        closeT<-collapse(an,">",value,"</",tag,">")
    else
        closeT<-collapse(an,closeT)
    collapse("<",tag,closeT)}


htmlDoc<-function(...,doctype="<!DOCTYPE html>",tag="html"){
    html<-htmlTags(tag,collapse(...))
    collapse(doctype,html)}

htmlTable<-function(x,...){
    br<-""
    shape<-dim(x)
    for(i in seq(shape[1]))
        {
            tr<-""
        for(j in seq(shape[2]))
            tr<-collapse(tr,htmlTags("td",x[i,j]))
            br<-collapse(br,htmlTags("tr",tr))
         }   
    htmlTags("table",br,...)}

htmlImage<-function(location,...){
    htmlTags("img",anotations=c(buildAnotations("src",location)),...)}

website<-htmlDoc(htmlTags("body",collapse(htmlTags("H1","Welcome"),
         htmlTags("p","This is a quick sample pagargaph"),
         htmlTable(matrix(c(1,2,3,4),2,2)))))

write(website,"test.html")






nb<-buildAnotations("cellpadding", "0")
cat(htmlTable(matrix(c(htmlTable(matrix(c(1,2,3),3,1),nb),htmlTable(matrix(c(1,2,3),3,1),nb),htmlTable(matrix(c(1,2,3),3,1),nb),htmlTable(matrix(c(1,2,3),3,1),nb)),2,2)))


matrix(list(list(1,2,3),list(1,2),list(1),list(1,2,3,4)),ncol=2)[1,1]
