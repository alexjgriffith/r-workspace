## work for zoe:

library("BSgenome.Scerevisiae.UCSC.sacCer3")
library("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
library("org.Sc.sgd.db")
library("Biostrings")

x<-org.Sc.sgd.db

columns(org.Sc.sgd.db)

geneL<-select(x,keys(x,"GENENAME"),c("CHR","CHRLOC","CHRLOCEND"),"GENENAME")


select(x, "ARG80", "ENTREZID", "GENENAME")[["ENTREZID"]]


txdb<-TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
    
txid <- select(txdb, eid, "TXNAME", "GENEID")[["TXNAME"]]


geneList<-select(txdb,keys(txdb,"GENEID"),c("EXONRANK","EXONCHROM","EXONSTART","EXONEND","EXONSTRAND","TXSTART","TXEND"),"GENEID")

geneList<-select(txdb,keys(txdb,"CDSID"),c("GENEID","EXONRANK","EXONCHROM","EXONSTART","EXONEND","EXONSTRAND"),"CDSID")


columns(txdb)

keytypes(txdb)

fasta<-getSeq(BSgenome.Scerevisiae.UCSC.sacCer3,geneList$EXONCHROM,start=geneList$EXONSTART,end=geneList$EXONEND,strand=geneList$EXONSTRAND)

cat(paste0(sapply(fasta,function(s) as.character(s)),collapse="\n"),file="~/Dropbox/clean.csv")

gregexpr("CAGCTA",fasta,ignore.case = TRUE)

## Get copy of CFP plasmid
