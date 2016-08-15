swapFunB<-stringToSwap("rpmi Leukemia tall_p3 Leukemia tall_p2 Leukemia tall_p1 Leukemia tall_p2_1 Leukemia tall_p2_2 Leukemia tall_p3_1 Leukemia tall_p3_2 Leukemia jurk_sandar_1 Leukemia jurk_sandar Leukemia jurk Leukemia jurk_1 Leukemia jurk_2 Leukemia rpmi_1 Leukemia rpmi_2 Leukemia cem_1 Leukemia cem_2 Leukemia cem Leukemia ecfc-tsa ECFC ecfc_old ECFC ecfc.tsa ECFC ecfc ECFC meka MEKA cd133 HSC cd34 HSC cd34_new HSC eryt Erythroid eryt_1 Erythroid eryt_2 Erythroid eryt_f Erythroid eryt_a Erythroid k562 Erythroid k562_1 Erythroid k562_2 Erythroid")
swapFun<-stringToSwap("rpmi RPMI tall_p3 Prima5 tall_p2 Prima2 tall_p1 Prima5 tall_p2_1 Prima2 tall_p2_2 Prima2 tall_p3_1 Prima5 tall_p3_2 Prima5 jurk_sandar_1 Jurkat jurk_sandar Jurkat jurk Jurkat jurk_1 Jurkat jurk_2 Jurkat rpmi_1 RPMI rpmi_2 RPMI cem_1 CEM cem_2 CEM cem CEM ecfc-tsa ECFC ecfc_old ECFC ecfc.tsa ECFC ecfc ECFC meka MEKA cd133 HSC cd34 HSC cd34_new HSC eryt Erythroid eryt_1 Erythroid eryt_2 Erythroid eryt_f Erythroid eryt_a Erythroid k562 K562 k562_1 K562 k562_2 K562")
swapFunC<-stringToSwap("Leukemia blue Erythroid red HSC orange ECFC green MEKA brown")
mockToNorm<-stringToSwap("cd133_mock cd133 cd34_mock cd34 cd34_new_mock cd34_new cem_mock cem_1 ecfc-tsa_mock ecfc-tsa eryt_a_mock eryt_a eryt_f_mock eryt_f jurk_mock jurk jurk_sandar_mock jurk_sandar k562_mock_1 k562_1 k562_mock_2 k562_2 meka_mock meka rpmi_mock_1 rpmi_1 tall_p1_mock tall_p1 tall_p2_mock tall_p2_1 tall_p3_mock tall_p2_1")
## combines Meka and HSC
swapFunD<-stringToSwap("rpmi Leukemia tall_p3 Leukemia tall_p2 Leukemia tall_p1 Leukemia tall_p2_1 Leukemia tall_p2_2 Leukemia tall_p3_1 Leukemia tall_p3_2 Leukemia jurk_sandar_1 Leukemia jurk_sandar Leukemia jurk Leukemia rpmi_1 Leukemia rpmi_2 Leukemia cem_1 Leukemia cem_2 Leukemia cem Leukemia ecfc-tsa ECFC ecfc_old ECFC ecfc.tsa ECFC ecfc ECFC meka HSC cd133 HSC cd34 HSC cd34_new HSC eryt Erythroid eryt_f Erythroid eryt_a Erythroid k562 Erythroid k562_1 Erythroid k562_2 Erythroid")

swapFunE<-stringToSwap("rpmi Leukemia tall_p3 LeukemiaP tall_p2 LeukemiaP tall_p1 LeukemiaP tall_p2_1 LeukemiaP tall_p2_2 LeukemiaP tall_p3_1 LeukemiaP tall_p3_2 LeukemiaP jurk_sandar_1 Leukemia jurk_sandar Leukemia jurk Leukemia rpmi_1 Leukemia rpmi_2 Leukemia cem_1 Leukemia cem_2 Leukemia cem Leukemia ecfc-tsa ECFC ecfc_old ECFC ecfc.tsa ECFC ecfc ECFC meka MEKA cd133 HSC cd34 HSC cd34_new HSC eryt ErythroidP eryt_f ErythroidP eryt_a ErythroidP k562 Erythroid k562_1 Erythroid k562_2 Erythroid")
swapFunC2<-stringToSwap("LeukemiaP turquoise ErythroidP tomato Leukemia blue Erythroid red HSC orange ECFC green MEKA brown")

p1<-plotPCMat2D(pca2Matr(env$prc),c("PC1","PC4"),categories,swapFun,swapFunB,swapFunC)
ggsave("~/Dropbox/pca_ggplot_20_1_4.png",p1)

png("~/Dropbox/pca_rplot_20_1_4.png")
plotPCs(env$prc$eigenVectors,c(1,4),env$prc$normData,categories)
dev.off()

p1<-plotPCMat2D(pca2Matr(env$prc),c("PC1","PC2"),categories,swapFun,swapFunB,swapFunC)
ggsave("~/Dropbox/pca_ggplot_20_1_2.png",p1)
png("~/Dropbox/pca_rplot_20_1_2.png")
plotPCs(env$prc$eigenVectors,c(1,2),env$prc$normData,categories)
dev.off()

tev<-pca2Matr(env$prc)
t2<-addColnames(cbind(tev[,1],tev[,2]*-1),c("PC1","PC2"))


p1<-plotPCMat2D(t2,c("PC1","PC2"),categories,swapFun,swapFunE,swapFunC2)
ggsave("~/Dropbox/pca_ggplot_20_1_2_rev_P.png",p1)


p1<-plotPCMat2D(pca2Matr(env$prc),c("PC1","PC4"),categories,swapFun,swapFunE,swapFunC2)
ggsave("~/Dropbox/pca_ggplot_20_1_4_P.png",p1)


png("~/Dropbox/pca_rplot_20_1_2.png")
plotPCs(env$prc$eigenVectors,c(1,2),env$prc$normData,categories)
dev.off()

regions<-makeGeneRegions(
    "/home/agriffith/Dropbox/UTX-Alex/jan/hg19.RefSeqGenes.csv",rna=rna)
geneList<-makeGeneList("/home/agriffith/Dropbox/UTX-Alex/jan/hg19.RefSeqGenes.csv",rna=rna)

genesGNN<-geneMatrix(env$over,env$reg,regions,geneList)


