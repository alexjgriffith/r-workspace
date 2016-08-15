categories<-as.character(unlist(read.table("~/Dropbox/Data/categories/22-categories.txt")))

contexts<-unique(swapFunD(categories))
