saveEnv<-function(env,folder,rdata=TRUE,overwrite=TRUE){
    if(rdata){
        save(env,folder)
    }
    else{
        if(dir.exists(folder) ){
            if(overwrite==FALSE){                        
                stop(paste0(folder," exists and saveEnv overwrite == FALSE."))
            }
        }
        else{
            dir.create(folder)
        }
        elements<-c("over","bed","categories","heights")
        saveF<-function(n)
            write.table(env[[n]],paste0(folder,"/",n),row.names=FALSE,quote=FALSE)
        mapply(saveF, elements)
        print(names(env))
        write.table(elements,paste0(folder,"/","elements"),
                    row.names=FALSE,quote=FALSE,col.names=FALSE)
    }
}

loadEnv<-function(folder,rdata=TRUE){
    if(rdata){
        load(folder,lenv)
    }
    else{
        elements<-list("over","bed","categories","heights")
        env<-lapply(elements,function(x)read.table(paste0(folder,"/",x),header=T))
        names(env)<-elements
        prc<-pca(env$heights)
        env<-c(env,list(prc=prc))
    }
    env
}
