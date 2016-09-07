## Install the necissary packages
library(parallel)

cs<-makeForkCluster(4)
packs<-list(
    "devtools",
    "functional",
            "abind",
            "ggplot2",
            "grid",
            "httr" )

parLapply(cs,packs,function(x) 
install.packages(x,type="source",repos='http://cran.us.r-project.org'))

   
source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
biocLite("BSgenome.Hsapiens.UCSC.hg19")


library(devtools)
install_github("alexjgriffith/CCCA")
install_github("alexjgriffith/mT1")
install_github("alexjgriffith/within")

