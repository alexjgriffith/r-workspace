source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/jan-functions.r")
source("~/r-workspace/feb-functions.r")
source("~/r-workspace/feb-variables.r")

source("~/Masters/CCCA/inst/scipts/eboxFrequency.r")

db=read.table("~/thesis-feb/motifAnnotations.table",header=TRUE)
## work on prefered distance stuff

env<-

eboxs<-genEboxCombs()

motifs<-unique(c(eboxs,as.character(db$motif[makeLogic(grep("_ALL",db$dataset),length(db$motif)) & 5==db$MACS & db$size==6])))

e2<-getMotifInfo(addFasta(getPRC5(1)),motifs)

farc<-addNames(outer(seql(e2$eboxLoc),seql(e2$eboxLoc),Vectorize(
    function(i,j){
        length(intersect(intersect(e2$eboxLoc[[i]],which(e2$reg[,"Leukemia"])),e2$eboxLoc[[j]]))
    }
)),motifs,motifs)


farc[1,sapply(seql(farc[,1]),
       function(l)
           order( (farc[,l]/sapply(e2$eboxLoc,length))[ seql(farc[,l])!=l] ))]


farc[1,sapply(seql(farc[,1]),
              function(l)order( farc[seql(farc[,l])!=l,l]))]

    takeFirst<-function(x){
        ret=rep(TRUE,length(x))
        first=c()
        for (i in seql(x)){
            if(x[i] %in% first)
                ret[i]<-FALSE
            else
                first<-c(first,x[i])
            
        }
        ret
    }

matrixToTable<-function(matrix){
    cnames<-colnames(matrix)
    rnames<-rownames(matrix)
    rap<-unlist(lapply(rnames,function(x) rep(x,length(cnames))))
    ret<-data.frame(a=rap, b=cnames ,dist=c(matrix))
    ret<-ret[ret$a!=ret$b,]
    sel<-takeFirst(mapply(function(a,b) {paste(sort(c(a,b)),collapse ="")},as.character(ret$a),as.character(ret$b)))
    ret[sel,]
}

ret<-matrixToTable(farc)


df<-data.frame()
for( m1 in motifs){
    for(m2 in motifs){
        df<-rbind(df,data.frame(m1=m1,m2=m2,dist=length(intersect(e2$eboxLoc[[m1]],e2$eboxLoc[[m2]]))))
    }
}

ret<-df[df$m1!=df$m2,]
sel<-takeFirst(mapply(function(a,b) {paste(sort(c(a,b)),collapse ="")},as.character(df$m1),as.character(df$m2)))
ret<-ret[sel,]


oret<-ret[rev(order(ret$dist)),]

i=5
print(oret[i,])
s<-e2$fasta[intersect(intersect(e2$eboxLoc[[oret[i,]$m1]],e2$eboxLoc[[oret[i,]$m2]]),which(e2$reg[,"Leukemia"]))]

mot<-as.character(oret[i,]$m1)
comp<-consenusIUPAC(compliment(IUPACtoBase(mot)))
abs(150-cbind(unlist(gregexpr(IUPACtoBase(as.character(oret[i,]$m1)) ,s,ignore.case = TRUE)),
unlist(gregexpr(IUPACtoBase(comp) ,s,ignore.case = TRUE))))

tocheck<-rev(order(sapply(seq(2,1000),function(i){
s<-e2$fasta[intersect(intersect(e2$eboxLoc[[oret[i,]$m1]],e2$eboxLoc[[oret[i,]$m2]]),which(e2$reg[,"Leukemia"]))]
b<-unlist(gregexpr(IUPACtoBase(mot) ,s,ignore.case = TRUE))
a<-getHeights(150-b[b!=-1])
max(abs(a))
   })))


i=170
s<-e2$fasta[intersect(intersect(e2$eboxLoc[[oret[i,]$m1]],e2$eboxLoc[[oret[i,]$m2]]),which(e2$reg[,"Leukemia"]))]
b<-unlist(gregexpr(IUPACtoBase(mot) ,s,ignore.case = TRUE))
a<-getHeights(150-b[b!=-1])
plot(a)

hist(sapply(gregexpr(IUPACtoBase(as.character(oret[i,]$m1)) ,s,ignore.case = TRUE), function(x) (150-x)[which.min(abs(150-x))]))
