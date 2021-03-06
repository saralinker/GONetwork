#' AssignCluster
#'
#' assign clusters for each gene based on GO distance. By default, it incorporates this "group" into the cytoscape file
#'
#' @param tab, output from cyto()
#' @param return_in_cytotable, returns the group numbers automatically into the original table provided from cyto(). Alternatively, you can set this to FALSE and simply return the group for each gene.
#' @param K, number of groups to subset genes into
#' @param rmv_stragglers, run a second round of clustering after removing the unclustered genes default = TRUE
#' 
#' @return .
#'
#' @export
#' 

AssignCluster <- function(tab,D, K = 10, rmv_stragglers = TRUE, cutoff = 0.3){
  if(rmv_stragglers == TRUE){
    k <-cutree(hclust(dist(D)),k = K) 
    g <- names(k[k>1])
    D2 <- D[g,g]
    k <-cutree(hclust(dist(D2)),k = K) 
    
    tmp <- melt(t(D2))
    tmp2 <- tmp[tmp$value > cutoff & tmp$value !=1,]
    tmp2$value <- 1 - tmp2$value
    colnames(tmp2) <- c("origin","destination","distance")
    tab <- tmp2
  }else{
    k <-cutree(hclust(dist(D)),k = K) 
  }
  tab$group <- as.numeric(k[match(as.character(tab$origin), as.vector(names(k)))])
  return(tab)
}

