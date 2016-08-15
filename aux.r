addColnames<-function(matrix,colnames){
    colnames(matrix)<-colnames
    matrix
}

addRownames<-function(matrix,rownames){
    rownames(matrix)<-rownames
    matrix
}

addNames<-function(matrix,colnames,rownames=colnames,list=NULL){
    if(is.null(list)){
        colnames(matrix)<-colnames
        rownames(matrix)<-rownames
    }
    else{
        names(matrix)<-colnames
    }
    matrix        
}
