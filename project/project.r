excellSheet<-function(matrix, GO){
    ework<-createWorkbook()
    s<-createSheet(ework,"Genes")
    for(i in seq(dim(matrix)[2])){
        addDataFrame(selectGenes(matrix,colnames(matrix)[i]),col.names=TRUE,row.names=FALSE,sheet=s,startColumn = i)
    }
    for(i in colnames(matrix)){
        if(i %in% contexts){
        s<-createSheet(ework,paste0("BP-",i))
        r<-addDataFrame(GO[[i]][[1]],col.names=TRUE,row.names=FALSE,sheet=s)
        print(paste0("MF-",i))
        d<-createSheet(ework,paste0("MF-",i))
        r<-addDataFrame(GO[[i]][[2]],col.names=TRUE,row.names=FALSE,sheet=d)
    }
    }
    ework    
}

makePeakFiles<-function(categories,mock="combined")paste("~/Dropbox/Data/august_peaks/",categories,"~",mock,"_mock_peaks.xls",sep="")


swapFunM<-stringToSwap("tall_p1 P2B tall_p2_1 P1 tall_p2_2 P1_R tall_p3_1 P2A tall_p3_2 P2A_R jurk_sandar_1 Jurk1 jurk_sandar Jurk1_R jurk_1 Jurk2 jurk Jurk2 jurk_2 Jurk2_R rpmi_1 RPMI rpmi_2 RPMI_R cem_1 CEM cem_2 CEM_R ecfc-tsa ECFC ecfc.tsa ECFC ecfc ECFC meka MEKA cd133 CD133 cd34 CD34_1 cd34_new CD34_2 eryt_1 ProEB_A_1 eryt_2 ProEB_A_1_R eryt ProEB_A_1 eryt_f ProEB_F eryt_a ProEB_A_2 k562_1 K562 k562_2 K562_R")
swapFunB<-stringToSwap("rpmi Leukemia tall_p3 Leukemia tall_p2 Leukemia tall_p1 Leukemia tall_p2_1 Leukemia tall_p2_2 Leukemia tall_p3_1 Leukemia tall_p3_2 Leukemia jurk_sandar_1 Leukemia jurk_sandar Leukemia jurk Leukemia jurk_1 Leukemia jurk_2 Leukemia rpmi_1 Leukemia rpmi_2 Leukemia cem_1 Leukemia cem_2 Leukemia cem Leukemia ecfc-tsa ECFC ecfc_old ECFC ecfc.tsa ECFC ecfc ECFC meka MEKA cd133 HSC cd34 HSC cd34_new HSC eryt Erythroid eryt_1 Erythroid eryt_2 Erythroid eryt_f Erythroid eryt_a Erythroid k562 Erythroid k562_1 Erythroid k562_2 Erythroid")
swapFun<-stringToSwap("rpmi RPMI tall_p3 Prima5 tall_p2 Prima2 tall_p1 Prima5 tall_p2_1 Prima2 tall_p2_2 Prima2 tall_p3_1 Prima5 tall_p3_2 Prima5 jurk_sandar_1 Jurkat jurk_sandar Jurkat jurk Jurkat jurk_1 Jurkat jurk_2 Jurkat rpmi_1 RPMI rpmi_2 RPMI cem_1 CEM cem_2 CEM cem CEM ecfc-tsa ECFC ecfc_old ECFC ecfc.tsa ECFC ecfc ECFC meka MEKA cd133 HSC cd34 HSC cd34_new HSC eryt Erythroid eryt_1 Erythroid eryt_2 Erythroid eryt_f Erythroid eryt_a Erythroid k562 K562 k562_1 K562 k562_2 K562")
swapFunC<-stringToSwap("Leukemia blue Erythroid red HSC orange ECFC green MEKA brown")
mockToNorm<-stringToSwap("cd133_mock cd133 cd34_mock cd34 cd34_new_mock cd34_new cem_mock cem_1 ecfc-tsa_mock ecfc-tsa eryt_a_mock eryt_a eryt_f_mock eryt_f jurk_mock jurk jurk_sandar_mock jurk_sandar k562_mock_1 k562_1 k562_mock_2 k562_2 meka_mock meka rpmi_mock_1 rpmi_1 tall_p1_mock tall_p1 tall_p2_mock tall_p2_1 tall_p3_mock tall_p2_1")
## combines Meka and HSC
swapFunD<-stringToSwap("rpmi Leukemia tall_p3 Leukemia tall_p2 Leukemia tall_p1 Leukemia tall_p2_1 Leukemia tall_p2_2 Leukemia tall_p3_1 Leukemia tall_p3_2 Leukemia jurk_sandar_1 Leukemia jurk_sandar Leukemia jurk Leukemia rpmi_1 Leukemia rpmi_2 Leukemia cem_1 Leukemia cem_2 Leukemia cem Leukemia ecfc-tsa ECFC ecfc_old ECFC ecfc.tsa ECFC ecfc ECFC meka HSC cd133 HSC cd34 HSC cd34_new HSC eryt Erythroid eryt_f Erythroid eryt_a Erythroid k562 Erythroid k562_1 Erythroid k562_2 Erythroid")


readAFS<-function(pvalue=5,control="single",directory="~/Dropbox/Data/AFS/",forward=22){
    fname<-paste0(directory,forward,"_pvalue_",pvalue,"_control_",control,".txt")
    fname
    read.table(fname,header=T)
}

readUDM<-function(pvalue=5,tr="treatment",control="single",directory="~/Dropbox/Data/UDM/",forward=22){
    fname<-paste0(directory,forward,"_",tr,"_pvalue_",pvalue,"_control_",control,".txt")
    fname
    read.table(fname,header=T)
}

getPRC<-function(pvalue=5,control="single",treatment="treatment",forward="22"){
    categories<-readCategories(paste0("~/Dropbox/Data/categories/",forward,"-categories.txt"))
    over<-orderBed(readAFS(pvalue,control,forward=forward))
    bed<-over[,1:3]
    heights<-readUDM(pvalue,treatment,control,forward=forward)
    prc<-pca(heights)    
    list(over=over,bed=bed,heights=heights,prc=prc,categories=categories)
}

getPRC20<-function(n=1,inner=0.25){
    out<-getPRC(20,"combined")
    add<-with(out,{
    ECFC<-CCCA:::normalize(prc$eigenVectors[,4])<(-(n))&CCCA:::normalize(prc$eigenVectors[,2])<(-(n))
    ECFCA<-CCCA:::normalize(prc$eigenVectors[,2])<(-n)
    ECFCB<-CCCA:::normalize(prc$eigenVectors[,4])<(-n)
    HSC<-CCCA:::normalize(prc$eigenVectors[,4])>(n)
    Leukemia<-CCCA:::normalize(prc$eigenVectors[,1])<(-n)
    Erythroid<-CCCA:::normalize(prc$eigenVectors[,1])>(n)
    ALL<-rep(TRUE,length(Erythroid))
    NONE<-CCCA:::normalize(prc$eigenVectors[,4])>(-inner)&
        CCCA:::normalize(prc$eigenVectors[,4])<(inner) &
            CCCA:::normalize(prc$eigenVectors[,1])>(-inner) &
                CCCA:::normalize(prc$eigenVectors[,1])<(inner) &
                    CCCA:::normalize(prc$eigenVectors[,2])>(-inner)     
    ECFC.alt<-ECFC & !(Erythroid | Leukemia | HSC)
    HSC.alt<-HSC & !(Erythroid | Leukemia | ECFC)
    Leukemia.alt<-Leukemia & !(HSC | Erythroid | ECFC)
    Erythroid.alt<-Erythroid & !(HSC | Leukemia | ECFC)
    reg<-addColnames(cbind(Erythroid.alt,Leukemia.alt,HSC.alt,ECFC.alt,NONE,ALL),c("Erythroid","Leukemia","HSC","ECFC","NONE","ALL"))
    list(reg=reg,ALL=ALL,NONE=NONE,ECFC=ECFC,HSC=HSC,Leukemia=Leukemia,Erythroid=Erythroid,
         ECFC.alt=ECFC.alt,HSC.alt=HSC.alt,Leukemia.alt=Leukemia.alt,Erythroid.alt=Erythroid.alt,
         ECFCB=ECFCB,ECFCA=ECFCA)

})
    append(out,add)
}


getPRC20F<-function(out,n=1,inner=0.25){
    add<-with(out,{
    ECFC<-function(prc=prc,n=n,inner=inner)
        {normalize(prc$eigenVectors[,4])<(-(n^2/2))&
             normalize(prc$eigenVectors[,2])<(-(n/2))}
    
    HSC<-function(prc=prc, n=n,inner=inner){
        normalize(prc$eigenVectors[,4])>(n)
    }
    Leukemia<-function(prc=prc2, n=n,inner=inner){
        normalize(prc$eigenVectors[,1])<(-n)
    }
    Erythroid<-function(prc=prc, n=n,inner=inner){
        normalize(prc$eigenVectors[,1])>(n)
    }
    ALL<-function(prc=prc, n=n,inner=inner){
        rep(TRUE,length(Erythroid))
    }
    NONE<-function(prc, n=n,inner=inner){
        normalize(prc$eigenVectors[,4])>(-inner)&
        normalize(prc$eigenVectors[,4])<(inner) &
            normalize(prc$eigenVectors[,1])>(-inner) &
                normalize(prc$eigenVectors[,1])<(inner) &
                    normalize(prc$eigenVectors[,2])>(-inner)
    }
    ECFC.alt<-function(prc=prc, n=n,inner=inner){
        ECFC(prc,n,inner) &
            !(Erythroid(prc,n,inner) | Leukemia(prc,n,inner) | HSC(prc,n,inner))
    }
    HSC.alt<-function(prc=prc, n=n,inner=inner){
        HSC(prc,n,inner) &
            !(Erythroid(prc,n,inner) | Leukemia(prc,n,inner) | ECFC(prc,n,inner))
    }
    Leukemia.alt<-function(prc=prc, n=n,inner=inner){
        Leukemia(prc,n,inner) &
            !(HSC(prc,n,inner) | Erythroid(prc,n,inner) | ECFC(prc,n,inner))
    }
    Erythroid.alt<-function(prc=prc, n=n,inner=inner){
        Erythroid(prc,n,inner) &
            !(HSC(prc,n,inner) | Leukemia(prc,n,inner) | ECFC(prc,n,inner))
    }
    reg<-addColnames(cbind(Erythroid.alt(prc,n,inner),
                           Leukemia.alt(prc,n,inner),
                           HSC.alt(prc,n,inner),
                           ECFC.alt(prc,n,inner),
                           NONE(prc,n,inner),
                           ALL(prc,n,inner)),c("Erythroid","Leukemia","HSC","ECFC","NONE","ALL"))
    list(reg=reg)
})
    append(out,add)
}
