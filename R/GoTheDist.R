GoTheDist <- function(M,Min = 4,minparents =   8, MinRow = 0, Threads = NULL){
  require(parallel)
  M2 <- M[,(as.vector(apply(M,2,function(x) sum(x > 0))) > Min)]
  parent <- as.numeric(apply(X=M2,MARGIN=2,max))
  M2 <- M2[,parent >= minparents]
  M2 <- M2[rowSums(M2) > MinRow,]
  if (!is.null(Threads)){
    cl <- makeCluster(getOption("cl.cores",  Threads))
    clusterExport(cl=cl, varlist=c("M2","InnerFun","cosine"))
    tmp <- parLapply(cl, X = c(1:nrow(M2)),fun = OuterFun)
    #tmp <- lapply( X = c(1:nrow(M2)),FUN = OuterFun)
    D2 <- matrix(unlist(tmp),nrow = nrow(M2), ncol = nrow(M2))
    rownames(D2) <- colnames(D2) <-  rownames(M2)
    stopCluster(cl)
  } else{
    D2 <- as.matrix(cosine(as.matrix(t(M2)))) #distance measure is cosine distance
  }
  D2[is.na(D2)] <- 0
  #Col <- colSums(D2)
  #Row <- rowSums(D2)
  rownames(D2) <- colnames(D2) <-  rownames(M2)
  return(D2)
}