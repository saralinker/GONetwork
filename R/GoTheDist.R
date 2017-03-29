#' Calculate the distance matrix between genes based on GO annotation
#'
#' .
#'
#' @param M matrix of GO terms for each gene 
#' @param Min Cutoff to exclude low frequency GO terms. Minimum number of genes that contain the GO term. default = 4.
#' @param minparents Cutoff to exclude non-informative GO terms (ex: nucleus vs. Pol-II transcription factor). defualt = 8.
#' @param MinRow 
#' 
#' @return matrix in a cytoscape-friendly format
#'
#' @export


GoTheDist <- function(M,Min = 4,minparents = 8, MinRow = 0, Threads = NULL ,method = "cosine"){
  requireNamespace("parallel", quietly = TRUE)
  requireNamespace("lsa", quietly = TRUE)
  M2 <- M[,(as.vector(apply(M,2,function(x) sum(x > 0))) > Min)]
  parent <- as.numeric(apply(X=M2,MARGIN=2,max))
  M2 <- M2[,parent >= minparents]
  M2 <- M2[rowSums(M2) > MinRow,]
  
  if(!is.null(Threads) & method != "cosine"){
    warning("Parallelization only currently available for cosine distances.")
  }
  
  if (!is.null(Threads) & method == "cosine"){
    cl <- makeCluster(getOption("cl.cores",  Threads))
    clusterExport(cl=cl, varlist=c("M2","InnerFun","cosine"))
    tmp <- parLapply(cl, X = c(1:nrow(M2)),fun = OuterFun)
    #tmp <- lapply( X = c(1:nrow(M2)),FUN = OuterFun)
    D2 <- matrix(unlist(tmp),nrow = nrow(M2), ncol = nrow(M2))
    rownames(D2) <- colnames(D2) <-  rownames(M2)
    stopCluster(cl)
  } else{
       if(method=="cosine"){
        D2 <- as.matrix(cosine(as.matrix(t(M2)))) #distance measure is cosine distance
        
      } else{
        D2 <- as.matrix(dist(as.matrix(M2),method = method))
      }
    }
  D2[is.na(D2)] <- 0
  rownames(D2) <- colnames(D2) <-  rownames(M2)
  return(D2)
}


