t<-function(){x=0;function(){y<-x;x<<-x+1;return(paste0(y,"> "))}}
makeActiveBinding(".prompt",t(),as.environment("package:base"))

repos=structure(c(CRAN="http://cran.utstat.utoronto.ca/"))
