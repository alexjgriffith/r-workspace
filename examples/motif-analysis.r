source("~/r-workspace/nov-functions.r")
source("~/r-workspace/dec-functions.r")
source("~/r-workspace/jan-functions.r")
source("~/r-workspace/feb-functions.r")
source("~/r-workspace/feb-variables.r")


qhwF<-function(fasta,r,p1,p2,pvalue,len=8,sd=0){
    print(paste("~/thesis-feb/",p1,"_",p2,"_pvalue=",pvalue,"_len=",len,"_sd=",sd,".motif",sep=""))
    homerWrapper(fasta,r[,p1],r[,p2],"~/Masters/mulcal/inst/lib/homer-4.7/bin/homer2",paste("~/thesis-feb/",p1,"_",p2,"_pvalue=",pvalue,"_len=",len,"_sd=",sd,".motif",sep=""),opts=paste("-S 25 -len ",len," > /tmp/homerTrash 2>&1 ",sep=""))
} 


qhwFALT<-function(fasta,r,p1,p2,pvalue,len=8,sd=0){
    print(paste("~/thesis-feb/",p1,"_",p2,"_pvalue=",pvalue,"_len=",len,"_sd=",sd,".motif",sep=""))
    homerWrapper(fasta,r[,p1],r[,p2],"~/Masters/mulcal/inst/lib/homer-4.7/bin/homer2",paste("~/thesis-feb/",p1,"_",p2,"_pvalue=",pvalue,"_len=",len,"_sd=",sd,"_alt.motif",sep=""),opts=paste("-S 25 -len ",len," > /tmp/homerTrash 2>&1 ",sep=""))
} 

env<-getPRC20()


callHomer<-function(env,pvalue,sd){
matrix<-with(env,{
    r<-cbind(ECFC=ECFC.alt,NONE=NONE,ALL=ALL,
             Erythroid=Erythroid.alt,Leukemia=Leukemia.alt,HSC=HSC.alt)
    for(b in c("NONE","ALL")){
        for(a in c("ECFC","HSC","Erythroid","Leukemia")){
            for(len in c(6,7,8))
                qhwF(fasta,r             
                    ,a,b,pvalue,len,paste0(sd,"_alt"))
        }
    }
})
}

callHomerALL<-function(env,pvalue,sd){
matrix<-with(env,{
    r<-cbind(ECFC=ECFC,NONE=NONE,ALL=ALL,
             Erythroid=Erythroid,Leukemia=Leukemia,HSC=HSC)
    for(b in c("NONE","ALL")){
        for(a in c("ECFC","HSC","Erythroid","Leukemia")){
            for(len in c(6,7,8))
                qhwF(fasta,r             
                    ,a,b,pvalue,len,paste0(sd,"_all"))
        }
    }
})
}



ahome<-function(){
    env<-getPRC20(0.25);
    env$fasta<-getSeq(BSgenome.Hsapiens.UCSC.hg19,env$bed$chr,start=env$bed$start+150,width=300)
    callHomer(env,21,0.25)
    callHomerALL(env,21,0.25)
    rm(env)
    env<-getPRC20(0.5);
    env$fasta<-getSeq(BSgenome.Hsapiens.UCSC.hg19,env$bed$chr,start=env$bed$start+150,width=300)
    callHomer(env,21,0.5)
    callHomerALL(env,21,0.5)
    rm(env)
    env<-getPRC20(2);
    env$fasta<-getSeq(BSgenome.Hsapiens.UCSC.hg19,env$bed$chr,start=env$bed$start+150,width=300)
    callHomer(env,21,2)
    callHomerALL(env,21,2)
    rm(env)

}


ahome<-function(){
    env<-getPRC20(2);
    env$fasta<-getSeq(BSgenome.Hsapiens.UCSC.hg19,env$bed$chr,start=env$bed$start+150,width=300)
    callHomer(env,20,2)
    callHomerALL(env,20,2)
    rm(env)
}

stampDF<-function(combs,pvalue,sd,dir="~/thesis-feb/")
    do.call(rbind,apply(combs,1,function(x,pvalue) data.frame(file=paste(dir,x[1],"_pvalue=",pvalue,"_len=",x[2],"_sd=",sd,".motif",sep=""),compatator=x[3],name=x[1],size=x[2],pvalue=pvalue,sd=sd),pvalue ))



combs<-rbind(
    stampDF(cbind(strsplit("Erythroid_NONE Erythroid_ALL Leukemia_NONE Leukemia_ALL HSC_NONE HSC_ALL ECFC_NONE ECFC_ALL"," ")[[1]],c(rep(6,8),rep(7,8),rep(8,8)),c(rep("TRANSFAC_Fams",24),rep("JASPAR_Fams",24))),pvalue=20,sd=1),      
       stampDF(cbind(strsplit("Erythroid_NONE Erythroid_ALL Leukemia_NONE Leukemia_ALL HSC_NONE HSC_ALL ECFC_NONE ECFC_ALL"," ")[[1]],c(rep(6,8),rep(7,8),rep(8,8)),c(rep("TRANSFAC_Fams",24),rep("JASPAR_Fams",24))),pvalue=20,sd=3,"~/thesis-feb/"),
    stampDF(cbind(strsplit("Erythroid_NONE Erythroid_ALL Leukemia_NONE Leukemia_ALL HSC_NONE HSC_ALL ECFC_NONE ECFC_ALL"," ")[[1]],c(rep(6,8),rep(7,8),rep(8,8)),c(rep("TRANSFAC_Fams",24),rep("JASPAR_Fams",24))),pvalue=5,sd=1),      
    stampDF(cbind(strsplit("Erythroid_NONE Erythroid_ALL Leukemia_NONE Leukemia_ALL HSC_NONE HSC_ALL ECFC_NONE ECFC_ALL"," ")[[1]],c(rep(6,8),rep(7,8),rep(8,8)),c(rep("TRANSFAC_Fams",24),rep("JASPAR_Fams",24))),pvalue=5,sd=3,"~/thesis-feb/")
    )


combs<-stampDF(cbind(strsplit("Erythroid_NONE Erythroid_ALL Leukemia_NONE Leukemia_ALL HSC_NONE HSC_ALL ECFC_NONE ECFC_ALL"," ")[[1]],c(rep(6,8),rep(7,8),rep(8,8)),c(rep("TRANSFAC_Fams",24),rep("JASPAR_Fams",24))),pvalue=20,sd=1)


for(i in combs$file){
    if(!file.exists(i))
        print(i)
}

combs<-stampDF(cbind(strsplit("Erythroid_NONE Erythroid_ALL Leukemia_NONE Leukemia_ALL HSC_NONE HSC_ALL ECFC_NONE ECFC_ALL"," ")[[1]],c(rep(6,8),rep(7,8),rep(8,8)),c(rep("TRANSFAC_Fams",24),rep("JASPAR_Fams",24))),pvalue=20,sd=3,"~/thesis-feb/")

i=1
x<-stampWrapper(as.character(combs$file[i]),combs[i,2],combs[i,3])

db<-multiStamp(combs)
write.table(db,"~/thesis-feb/motifAnnotations_20_3-alt.table",quote=FALSE,row.names=FALSE)


combs<-stampDF(cbind(strsplit("Erythroid_NONE Erythroid_ALL Leukemia_NONE Leukemia_ALL HSC_NONE HSC_ALL ECFC_NONE ECFC_ALL"," ")[[1]],c(rep(6,8),rep(7,8),rep(8,8)),c(rep("TRANSFAC_Fams",24),rep("JASPAR_Fams",24))),pvalue=20,sd="3_alt","~/thesis-feb/")

db<-multiStamp(combs)
write.table(db,"~/thesis-feb/motifAnnotations_20_3-alt.table",quote=FALSE,row.names=FALSE)


combs<-stampDF(cbind(strsplit("Erythroid_NONE Erythroid_ALL Leukemia_NONE Leukemia_ALL HSC_NONE HSC_ALL ECFC_NONE ECFC_ALL"," ")[[1]],c(rep(6,8),rep(7,8),rep(8,8)),c(rep("TRANSFAC_Fams",24),rep("JASPAR_Fams",24))),pvalue=20,sd="3","~/thesis-feb/")
db<-multiStamp(combs)
write.table(db,"~/thesis-feb/motifAnnotations_20_3.table",quote=FALSE,row.names=FALSE)
combs<-stampDF(cbind(strsplit("Erythroid_NONE Erythroid_ALL Leukemia_NONE Leukemia_ALL HSC_NONE HSC_ALL ECFC_NONE ECFC_ALL"," ")[[1]],c(rep(6,8),rep(7,8),rep(8,8)),c(rep("TRANSFAC_Fams",24),rep("JASPAR_Fams",24))),pvalue=20,sd="1","~/thesis-feb/")
db<-multiStamp(combs)
write.table(db,"~/thesis-feb/motifAnnotations_20_1.table",quote=FALSE,row.names=FALSE)
combs<-stampDF(cbind(strsplit("Erythroid_NONE Erythroid_ALL Leukemia_NONE Leukemia_ALL HSC_NONE HSC_ALL ECFC_NONE ECFC_ALL"," ")[[1]],c(rep(6,8),rep(7,8),rep(8,8)),c(rep("TRANSFAC_Fams",24),rep("JASPAR_Fams",24))),pvalue=20,sd="2","~/thesis-feb/")
db<-multiStamp(combs)
write.table(db,"~/thesis-feb/motifAnnotations_20_2.table",quote=FALSE,row.names=FALSE)
combs<-stampDF(cbind(strsplit("Erythroid_NONE Erythroid_ALL Leukemia_NONE Leukemia_ALL HSC_NONE HSC_ALL ECFC_NONE ECFC_ALL"," ")[[1]],c(rep(6,8),rep(7,8),rep(8,8)),c(rep("TRANSFAC_Fams",24),rep("JASPAR_Fams",24))),pvalue=20,sd="2_alt","~/thesis-feb/")
db<-multiStamp(combs)
write.table(db,"~/thesis-feb/motifAnnotations_20_2-alt.table",quote=FALSE,row.names=FALSE)



components<-c("dataset","genes","motif","rank","order","escore")

db[grep("smad",db$genes,ignore.case=TRUE),components]

db[grep("runx",db$genes,ignore.case=TRUE),components]

db[grep("fox",db$genes,ignore.case=TRUE),components]

db[grep("lmo",db$genes,ignore.case=TRUE),components]

db[grep("tal1",db$genes,ignore.case=TRUE),components]

db[grep("p53",db$genes,ignore.case=TRUE),components]

db[grep("egr",db$genes,ignore.case=TRUE),components]



get2<-function(clist){
    sapply(strsplit(as.character(clist),"_"),"[[",3)
}


db<-read.table("~/thesis-feb/motifAnnotations_20_2.table")

uin<-unique(get2(db[db$comparator=="TRANSFAC_Fams",]$genes))

gef<-data.frame(row.names=NULL,motif=uin,count=sapply(uin,function(i) length(grep(i,db$genes))))

ogef<-gef[order(gef[,2],decreasing=TRUE),]

data.frame(row.names=NULL,motif=sapply(ogef[,1],function(x)grep(x,db$genes,value=TRUE)[1]),count=ogef[,2])

Transcription factors of interest
E2F-1
EMP1
E47
SOX
TIG
LMO2
PAC4
HAC1
EGR
HAC1
YY1

ERF2
SP1
ABI4
E2F
PDR3
GAMYB
PU.1
HNF4

printSel<-function(x)
db[grep(x,db$genes,ignore.case=TRUE),components]

env<-getPRC(20,"combined")
a<-addRownames(do.call(rbind,lapply(seq(0.25,3,by=0.25),function(X) apply(getPRC20F(env,X)$reg,2,function(x) length(which(x))))),paste0(seq(0.25,3,by=0.25),"SD"))



cat(makeLatexTable(a))

db[db$motif=="GACAGC",components]


# use pvalue = 20 , sd = 2
"Smad3"
"Erythroid_ALL SMAD_SMAD3_M00701   GACAGC    4    10 7.5593e-06"
"Erythroid_ALL SMAD_SMAD3_M00701   GACAGC    4    10 7.5593e-06"
"Erythroid_NONE SMAD_SMAD3_M00701  GGCTGTC    1     5 7.6310e-07"

"HSC_ALL   SMAD_MAD_M01090   CCGWCG    1     1 6.2983e-04"
"HSC_ALL SMAD_SMAD3_M00701  GGKCTGT    1     2 3.0491e-05"

"Sox"
ECFC_ALL  HMG_SOX5_M00042   AACAAT    2     9 1.2208e-08
ECFC_ALL         HMG_Sox5   AACAAT    1     9 8.6265e-10
ECFC_NONE         HMG_Sox5 AAACAATG    4     6 1.0199e-06

ECFC_ALL HMG_SOX17_M01016   AACAAT    4     9 1.1858e-07
ECFC_NONE        HMG_Sox17 AAACAATG    5     6 9.3562e-06
Erythroid_NONE        HMG_Sox17 GCCTTGTC    1     3 4.6791e-04
Erythroid_ALL        HMG_Sox17 GGACAAGG    2     3 7.4812e-05

"Runx"
Leukemia_NONE        runt_AML1a_M00271   ACCACA    1     1 1.2650e-07
Leukemia_ALL        runt_AML1a_M00271   ACCACA    1     1 1.2595e-07
Leukemia_NONE        runt_AML1a_M00271   AACCRC    1     3 1.1475e-06
Leukemia_ALL          runt_AML_M00769   ACCACA    3     1 1.6408e-07
Leukemia_ALL          runt_AML_M00769  AACCACA    1    13 4.8159e-09
Leukemia_NONE         runt_AML1_M00751 GHTGTGGT    1     1 1.4063e-07

HSC_NONE               RUNT_RUNX1   CCACAG    1     2 3.4263e-04
HSC_ALL               RUNT_RUNX1   CGACAG    4     3 1.5114e-02


"GATA"
"Erythroid_NONE     CC_GATA-1_M00347   TTATCT    2     1 5.6036e-08"
"Erythroid_ALL     CC_GATA-1_M00347   CTTATC    2     1 2.1597e-08"
Erythroid_NONE     CC_GATA-1_M00347  CTTATCT    2     1 4.1304e-10
Erythroid_ALL     CC_GATA-1_M00347  CTTATCT    2     1 3.6813e-10
Erythroid_NONE     CC_GATA-1_M00347 SSCTTATC    2     1 2.1762e-08
Erythroid_ALL     CC_GATA-1_M00347 YCTTATCW    2     1 2.8077e-09

Erythroid_NONE     CC_GATA-2_M00348  CTTATCT    1     1 2.1594e-10
Erythroid_ALL     CC_GATA-2_M00348  CTTATCT    1     1 2.0670e-10

ECFC_NONE     CC_GATA-1_M00347   GATAAG    3     6 4.1107e-08
ECFC_ALL     CC_GATA-1_M00346   TGTTAT    1     5 1.0423e-07


"Fox"
ECFC_ALL  fork_FOXO3_M00477   TTGTTT    2     7 1.9661e-07
ECFC_ALL  fork_FOXO4_M00476  TGTTTTC    1     9 8.6309e-06
ECFC_NONE     FORKHEAD_FOXD1  TGTTTTC    1     2 1.5729e-06
ECFC_NONE    fork_FOX_M00809 AAACAATG    5     6 2.0321e-06
ECFC_NONE     FORKHEAD_Foxq1 AAACAATG    3     6 5.9407e-07

"YY1"
Erythroid_ALL CH_YY1_M00069 GGCTGCT    1     7 7.3594e-06

"EGR"
HSC_NONE   CH_Egr_M00807 GGGGGC    2     9 6.0177e-07
HSC_ALL   CH_Egr_M00807 CGCCCC    4     4 8.2525e-07


"LMO2"
Erythroid_ALL LIM_Lmo2_M00278   CTATCT    1     2 1.7875e-07
Erythroid_NONE LIM_Lmo2_M00278  YTATCTG    2     2 5.6578e-07

"HAC1"
ECFC_NONE bZIP_HAC1_M00730 GCCAGCG    1     6 2.0737e-08


"P54"
ECFC_NONE ETS_c-Ets-1_p54_M00032 GCAKCCGG    1    14 3.1171e-07
ECFC_ALL ETS_c-Ets-1_p54_M00032 GCATCCGG    1     9 4.0645e-08
ECFC_NONE ETS_c-Ets-1_p54_M00032   CTTCCG    3     2 2.4705e-07
ECFC_ALL ETS_c-Ets-1_p54_M00032   GGAAGC    3     3 7.2113e-06


"ebox"
Leukemia_ALL bHLH_Tal-1alpha-E47_M00066   CATCTG    3     9 5.7107e-07
Leukemia_ALL            bHLH_E47_M00071  ACACCTG    2     1 8.2121e-09
Leukemia_NONE            bHLH_E47_M00002 CGCAGGTG    3    10 3.7507e-09
Leukemia_ALL            bHLH_E47_M00002 GCAGGTGT    1     3 6.6893e-11
Leukemia_ALL bHLH_Tal-1alpha-E47_M00066 CAGCTGTT    2    12 1.0913e-07
Leukemia_ALL            bHLH_TAL1_M00993   CATCTG    1     9 9.8515e-08
Leukemia_NONE            bHLH_MyoD_M00184   CACCTG    1    12 2.4826e-08
Leukemia_ALL             bHLH_E12_M00693   CACCTG    2     8 3.1914e-08
Leukemia_NONE             bHLH_E2A_M00804 CGCAGGTG    4    10 7.1975e-08

HSC_NONE      bHLH_Hand1-E47_M00222  CAGAYGS    1     3 2.4332e-06
HSC_NONE   bHLH_TAL1-TCF3 CAGCTGBY    2    15 1.7642e-06
HSC_NONE             bHLH_HEB_M00698 CAGCTGBY    2    15 8.6559e-07


Erythroid_NONE        bHLH_myogenin_M00712 CAGCHGCC    2     7 5.8782e-07
Erythroid_NONE        bHLH_myogenin_M00712 TGCWGCTG    3     8 6.8565e-06


"erf2"
 Leukemia_NONE AP2_ERF2_M01057   CGCCGC    1    10 2.9894e-07
Leukemia_NONE AP2_ERF2_M01057   GGCGGC    1    11 8.9194e-09

ECFC_NONE AP2_ERF2_M01057   CGCCGC    1     9 1.1394e-06

HSC_NONE AP2_ERF2_M01057   GCGGCG    1    10 3.9818e-06
HSC_ALL AP2_ERF2_M01057  ACCGCYA    1     4 4.9065e-06


"ABI4"
Leukemia_NONE AP2_ABI4_M00958   GCACCS    1     5 9.7680e-06
Leukemia_ALL AP2_ABI4_M00958   CGCACC    1     4 8.6274e-06
Leukemia_NONE AP2_ABI4_M00958  GCGCCGC    1     8 1.5781e-06
Leukemia_ALL        AP2_ABI4 GCGGCACC    1    15 5.7651e-06

HSC_ALL AP2_ABI4_M00958 CGGCGGTG    1     5 1.1570e-08

Erythroid_NONE        AP2_ABI4 CGGAGCCC    1    11 8.7490e-06


"PDR3"
ECFC_NONE C6_PDR3_M00752 TTCCGCGG    1     2 5.1775e-11
HSC_NONE C6_PDR3_M00752 CACCGCGG    2     6 2.4218e-07

"GAMYB"
Leukemia_NONE  trp_GAMYB_M00345 AACCGCCG    1    17 5.8997e-09
Leukemia_NONE TRP-CLUSTER_GAMYB AACCGCCG    1    17 3.8521e-07


"PU.1"
HSC_ALL ETS_PU.1_M00658  CGTCCTC    1    10 9.3626e-07

ECFC_NONE ETS_PU.1_M00658  MCTTCCT    4     7 1.2660e-07
ECFC_ALL ETS_PU.1_M00658  GCTTCCT    3     3 4.7601e-08
ECFC_NONE ETS_PU.1_M00658 AGGAAGCS    3    17 3.3072e-07
ECFC_ALL ETS_PU.1_M00658 CGGAGGAA    1     6 8.8292e-07

"Nuclear Receptor"
ECFC_NONE        NUCLEAR_RECEPTOR_usp   TGACCC    1     8 9.0780e-08
ECFC_ALL     NUCLEAR_RECEPTOR_RORA_1   TGACCT    1     4 7.8015e-08
ECFC_ALL   NUCLEAR_RECEPTOR_RXRA-VDR  RGGTCAW    1     8 8.2547e-08
NUCLEAR_RECEPTOR_NR1H2-RXRA TGACCCTT    1    13 1.0971e-06
ECFC_ALL      NUCLEAR_RECEPTOR_NR2F1 DRAGGTCA    1     8 1.6971e-07

HSC_ALL              CH_NRSE_M00325  GGCGCTC    1     3 5.1484e-06
HSC_NONE              CH_NRSF_M01028 CACCGCGG    1     6 1.8598e-07


stampDF<-function(combs,pvalue,sd,dir="~/thesis-feb/")
    do.call(rbind,apply(combs,1,function(x,pvalue) data.frame(file=paste(dir,x[1],"_pvalue=",pvalue,"_len=",x[2],"_sd=",sd,".motif",sep=""),compatator=x[3],name=x[1],size=x[2],pvalue=pvalue,sd=sd),pvalue ))


combs<-stampDF(cbind(strsplit("Erythroid_NONE Erythroid_ALL Leukemia_NONE Leukemia_ALL HSC_NONE HSC_ALL ECFC_NONE ECFC_ALL"," ")[[1]],c(rep(6,8),rep(7,8),rep(8,8)),c(rep("TRANSFAC_Fams",24),rep("JASPAR_Fams",24))),pvalue=21,sd="0.5_all","~/feb-homer/")


for(i in combs$file){
    if(file.exists(i))
        print(i)
}


combs<-stampDF(cbind(strsplit("Erythroid_NONE Erythroid_ALL Leukemia_NONE Leukemia_ALL HSC_NONE HSC_ALL ECFC_NONE ECFC_ALL"," ")[[1]],c(rep(6,8),rep(7,8),rep(8,8)),c(rep("TRANSFAC_Fams",24),rep("JASPAR_Fams",24))),pvalue=21,sd="0.5_alt","~/feb-homer/")
db<-multiStamp(combs)
write.table(db,"~/feb-homer/motifAnnotations_20_0.5.table",quote=FALSE,row.names=FALSE)
combs<-stampDF(cbind(strsplit("Erythroid_NONE Erythroid_ALL Leukemia_NONE Leukemia_ALL HSC_NONE HSC_ALL ECFC_NONE ECFC_ALL"," ")[[1]],c(rep(6,8),rep(7,8),rep(8,8)),c(rep("TRANSFAC_Fams",24),rep("JASPAR_Fams",24))),pvalue=21,sd="0.25_alt","~/feb-homer/")
db<-multiStamp(combs)
write.table(db,"~/feb-homer/motifAnnotations_20_0.5.table",quote=FALSE,row.names=FALSE)
