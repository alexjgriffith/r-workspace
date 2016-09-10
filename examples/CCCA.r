library(CCCA)
source("~/r-workspace/project/project.r")
source("~/r-workspace/project/project-variables.r")
source("~/r-workspace/project/ccca.r")

library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

ccca<-getPRC20(2)



ccca<-addFasta(ccca,BSgenome.Hsapiens.UCSC.hg19)

save(ccca,file="~/Dropbox/Data/ccca_20.RData")

system("du -h ~/Dropbox/Data/ccca_20.RData")

system("du -h fname.xz")

system("scp ~/Dropbox/Data/ccca_20.RData m4l:~/")

con <- pipe("xz -T8 -6 -e > fname.xz", "wb")

save(ccca,file=con)


addFasta<-function (env, genome = BSgenome.Hsapiens.UCSC.hg19, width = 300) 
{
    if (is.null(env$bed$chr) | is.null(env$bed$start)) 
        stop("addFasta env list must contain bed$chr bed$start and bed$end")
    if (!require(Biostrings)) 
        stop("Must install the Biostrings package from Bioconductor.\nsource(\"https://bioconductor.org/biocLite.R\"); biocLite(\"Biostrings\")")
    env$fasta <- getSeq(genome, env$bed$chr, start = (env$bed$start + 
        env$bed$end)/2 - floor(width/2), width = width)
    env
}

library(mT1)
load("ccca_20.RData")
motifs<-unique(c("CANNTG","GATAA",mT1_jaspar$jsublM))
objMT1<-mT1(ccca$fasta,motifs,verbose=TRUE)
save(objMT1)
