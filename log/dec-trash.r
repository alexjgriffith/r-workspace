#### dec-1
# debugging
sapply(sapply(unlist(loadPWM("~/thesis-november/PC1-1_PC1+1_pvalue=5_len=6.motif")[,1]),function(x) substr(x,2,nchar(x))),function(x,y) x%in%y, valueanot[pvalueanot$filename=="~/thesis-november/PC1-1_PC1+1_pvalue=5_len=6.motif"&pvalueanot$comparator=="TRANSFAC_Fams"&pvalueanot$rank==1,"motif"])

