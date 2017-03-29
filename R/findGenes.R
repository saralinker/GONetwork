findGenes <- function(M,term = "nuclear"){
  tmp <- M[,grep(term,colnames(M),ignore.case=TRUE)]
  if(!is.null(dim(tmp))){
    tmp2 <- tmp[rowSums(tmp) > 0,]
    return(rownames(tmp2))
    
  }else{
    #tmp2 <- tmp[rowSums(tmp) > 0,]
    #tmp2 <- as.list(tmp2)
    #names(tmp2) <- rownames(M[tmp> 0,])
    return(rownames(M[tmp> 0,]))
  }
  
}