# R Workspace

## Development work done for CCCA

## Alexander Griffith MASc, UOttawa

### Version 0.1

* Functions are labled *month*-function.r
* Variables are labled *month*-variables.r

### Making the AFS

```R
## categories
cats<-c("jurk","eryt","cd34","cem_1")
## list of peak locations (in MACS' xls output)
peakfiles<-
rawDatafiles<-
```

```R
## Generate the AFS from peakfiles, uses cats to name the columns
## Restricts the combined peaks to a p-vale of 1E-20
## NOTE: these peaks have to be called in advance
AFS<-makeAFS(peakfiles,cats,pValue=20)
UDM<-peakDensity(AFS,rawDatafiles,n=4)
prc<-pca(UDM)
```

```R
env<-env(
	over=AFSm
	bed=AFS[,c("chr","start","end")]
	heights=UDM
	prc=prc)
env$reg<-with(env,{
	Leukemia<-CCCA:::normalize(prc$eigenVEctors[,1])<1
	Erythroid<-CCCA:::normalize(prc$eigenVEctors[,1])>1
	HSC<-CCCA:::normalize(prc$eigenVEctors[,2])<1
	cbind(Leukemia=Leukemia & !(Eythroid | HSC),
		Erythroid=Erythroid & !(Leukemia | HSC),
		HSC=HSC & !(Erythroud| Leukemia)
})
```

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
