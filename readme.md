
# R Workspace

## Development work done for CCCA

## Alexander Griffith MASc, UOttawa

___
### Version 0.1
This repo contains working examples of the application of [CCCA](https://github.com/alexjgriffith/CCCA). Not all of this is stable, and a lot of the functionality presented here relies on functions from the [project/ccca.r]() folder. Explicit examples can be found in the example folder and day to day use of the package can be found in the log folder. Note that much of the source is undocumented. I have included this documentation to cover some of the core example applications of CCCA.

___
### Libraries and PreReqs
The base CCCA pacakge requires few dependencies, however for usefull downstream analysis several bioconductor packages are recomended, including Biostrings, BSgenome.Hsapiens.UCSC.hg19, RDAVIDWebService, and org.Hs.eg.db.
```R
library(CCCA)
library(xlsx) # Bindings for the JAVA xlsx package
## from within the r-workspace repo
source(project/ccca.r) # extra wrappers around CCCA
source(project/project.r) # project specific details
## We are going to generate our source variables from  scratch
## source(project/project-variables.r)
```
___
### Making the AFS

* Generate the AFS from peakfiles, uses cats to name the columns
* This sample restricts the combined peaks to a p-vale of 1E-20
* NOTE: The peaks have to be called in advance, I used MACS 2.01

```R
## categories
cats<-c("jurk","eryt","cd34","cem_1")
## list of peak locations (in MACS' xls output)
peakfiles<-paste0("samplePeaks/",cats,"~combined.xls")
## Calling the AFS function with locations to the peak files
## and categories for naming
AFS<-makeAFS(peakfiles,cats,pValue=20)
```
___
### Making the UDM (Unified Density Matrix)
The UDM is a matrix of peak pile up counts. The core is implemented in c and can be run in parallel to speed up analysis. Possible future improvements include reading from BAM files rather than BED files, and using memory map rather tha fread.

```R
## locations of raw reads that have been mapped to the genome
rawDatafiles<-paste0("sampleReads/",cats,".bed")
## generating the UDM using 4 cores
UDM<-peakDensity(AFS,rawDatafiles,n=4)
## pca is a wrapper around `prcomp` that applies qn and returns
## both the PCs and the normalized data
prc<-pca(UDM)
## Create an env named env with all of the compenents we have
## generated so far
env<-env(
	over=AFSm
	bed=AFS[,c("chr","start","end")]
	heights=UDM
	prc=prc)
	
```
___
### Ploting the Weighted Locations of the Data Sets
Coming soon ...
___
### Isolating Peaks
```R
## Define the regions that are most important for each condition	
env$reg<-with(env,{
	Leukemia<-CCCA:::normalize(prc$eigenVEctors[,1])<1
	Erythroid<-CCCA:::normalize(prc$eigenVEctors[,1])>1
	HSC<-CCCA:::normalize(prc$eigenVEctors[,2])<1
	cbind(Leukemia=Leukemia & !(Eythroid | HSC),
		Erythroid=Erythroid & !(Leukemia | HSC),
		HSC=HSC & !(Erythroud| Leukemia)
})
```

___
### Motif De Novo Search
Coming soon ....

___
### Ebox Locations
Coming soon ....

___
### Preferred Distances
Coming soon ....
___
### Gene Association
For gene assosiation I implemented the default [GREAT](http://bejerano.stanford.edu/great/public/html/) algorithm in R. I use the `org.Hs.eg.db` library to gernete the entrez ids (EZ IDs) from the official gene symbols. This example uses `hg19.RefSeqGenes.csv` for its gene list. You will have to adapt the process based on the gene list you want to use.

* Note that `geneRegion` is currently reliable but highly unoptimized. I am currently updating a version with improved speed.
* `geneRegion` only has to be called once for any geneList provided, I recomend saving the resulting list using R's build in `save` function and loading in the results when needed.

```R
library(org.Hs.eg.db) # used to transform GENE NAMES to EZ IDS
## Gene file, change if not for hg19
geneFile<-system.file("exdata","hg19.RefSeqGenes.csv",package="CCCA")
geneList<-read.delim(geneFile)
## Take a subset of the geneList <chr><tss><strand>
geneDF<-local({
	t<-list()
	t$chr<-as.vector(geneList$chrom)
	t$strand<-geneList$strand
	## switch the TSS to the txEnd for genes along the -ve strand
	t$tss<-geneList$txStart
	t$tss[t$strand=="-"]<-geneList$txEnd[t$strand=="-"]
	as.data.frame(t)
})
## Define the genomic regions, accomplished using the default
## GREAT algorithm
regions<-genomicRegions(geneDF$chr,geneDF$tss,
	geneDF$strand,1000,5000,100000)
## Assosiate Genes with peaks in env
## "name2" is the column name of geneList being returned
genes<-geneMatrix(env$over,env$reg[,contexts],
	regions,geneList,id="name2")
## Find the EZ IDs for use with davidWebService	
genesEZ<-local({
	fun<-function(x){
		symbolToENTREZID(unique(selectGenes(genes,x,df=FALSE)))
	}
	addNames(lapply(contexts,fun),contexts,list=TRUE)
})


```

___
### David from R
David provides a SOAP interface for its Gene Ontology(GO) functions. The R package to interface with the DAVID servers is RDAVIDWebService. Note as of writing this RDAVIDWebService requires JAVA8. See the **David** **Web** **Service** **WorkArround** section for details.

I also use xlsx to save the results in excell format.

Before Continuing

2. Visit the [David Webservice Registration](https://david.ncifcrf.gov/webservice/register.htm) page and sign up using your institution email
  * note that your gmail account is not sufficient to make an account
1. Make sure R knows to use Java 8
4. RDAVIDWebService can be installed from bioconductor 
   ```R
   source("http://bioconductor.org/biocLite.R")
   biocLite("RDAVIDWebService")

   ```
5. Generate the `genesEZ` variable as discussed in the **Gene** **Assosiation** section.


```R
library(RDAVIDWebService)
library(xlsx)

## email registered with DAVID
user<-"user@example.com"

contexts<-colnames(genes)
## Find the BP and MF GO terms for each context
GOres<-genGO(user,contexts,genesEz)

## Functions for swaping the EZ ids with the official gene symbols

#' Replace EZ ID
#' @description takes a table returned by david and swaps the EZ IDs
#' with official gene symbols
#' @param table David GO table
#' @param pvalue minimum cut off for keeping GO terms
#' @returns a table of GO terms with IDs replaced with names
#' @examples
#' genesEZ<-list(list("GATA","TAL1","RUNX1"))
#' contexts<-c("test")
#' user<-testuser@example.com
#' davidResults<-genGO(user,contexts,genesEZ)
#' ## Replace the ids with names and only keep the 
#' ## GO terms that have a p-value less than 1e-3
#' replaceEZID(davidResults[[1]][[1]],pvalue=1e-3)
#' @export
replaceEZID<-function(table,pvalue=1e-1){
	## Takes a string of EZ ids and splits it. Then the ids are 
	## swapped with gene names and returns a string of gene names
	ezidGene<-function(ezid){
		ezidSplit<-strsplit(ezid,", ")[[1]]
		geneName<-select(org.Hs.eg.db,ezidSplit, 
			"SYMBOL" ,"ENTREZID")[,2]
		paste(geneName,collapse=", ")
	}
    sigReg<-table[,5]<pvalue
    geneNames<-sapply(table[sigReg,6],ezidGene)
    stable<-table[sigReg,]
    stable[,6]<-geneNames
    stable
}

genesWNames<-lapply(contexts,
	function(x){ 
		list(BP=replaceEZID(GOres[[x]][[1]]),
			MF=replaceEZID(GOres[[x]][[2]]))})
			
## Save the results 
geneGOExcell<-createWorkbook()
geneGOExcell<-excellSheet(genes,genesWNames)
saveWorkbook(geneGOExcell,path.expand("~/Desktop/geneGO.xlsx"))

```

___
### David Web Service Workaround
As of writing this (aug 19 216) the david package requires JAVA8. They have also shifted from using a insicure HTTP to a TLS HTTPS. Make sure that your R environment is aware of the correct version of JAVA.

From the command line
```sh
R CMD javareconf -n
```

From within R
```R
library(rJava)
.jinit()
.jcall("java/lang/System", "S", "getProperty", "java.runtime.version")

```
If you are connected to the 1.8 runtime you should be good to go, otherwise update your version of JAVA before atempting to use davidWebService. 

___
### Notes
Specific tasks should have their own names. e.g. those that provide examples to gene assosiation or motif denovo. The examples should be self contained and when used transcribed into the specific *month*-*day*.r file appropriate for the day. There should be a refference included to any examples and as functions become codified add Roxygen documentation and integrate them into the CCCA package.


Following the reorganize commit things may be broken.

To fix switch:

```R
source("~/r-workspace/nov-functions.r")
source("~/r-workspace/nov-variables.r")
```

To:

```R
source("~/r-workspace/functions/nov-functions.r")
source("~/r-workspace/variables/nov-variables.r")
```

___
<center>&copy; Alexander Griffith, Ottawa, Canada, 2016</center>
