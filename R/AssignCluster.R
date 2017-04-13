#' AssignCluster
#'
#' assign clusters for each gene based on GO distance. By default, it incorporates this "group" into the cytoscape file
#'
#' @param tab, output from cyto()
#' @param cutoff, the minimum distance between nodes to be considered a cluster
#' @param return_in_cytotable, returns the group numbers automatically into the original table provided from cyto(). Alternatively, you can set this to FALSE and simply return the group for each gene.
#' 
#' @return .
#'
#' @export
#' 

AssignCluster <- function(tab, cutoff = 0.1, return_in_cytotable = FALSE){
  cutoff <- 1 - cutoff
  tab$origin <- as.character(tab$origin)
  tab$destination <- as.character(tab$destination)
  allgenes <- unique(c(tab$origin, tab$destination))
  k.out <- data.frame(genes = allgenes, k = NA)
  j <- 0
  for (g in allgenes){
    j <- j + 1
    a <- tab[tab$origin == g & tab$distance < cutoff | tab$destination == g & tab$distance < cutoff,]
    a <- unique(c(a$origin, a$destination))
    k.tmp <- k.out[match(a, k.out$genes),"k"]
    if(sum(is.na(k.tmp)) == length(a)){
      k.out[match(a, k.out$genes),"k"] <- j
    }else{
      for(old.k in unique(k.tmp[!is.na(k.tmp)])){
        k.out[k.out$k == old.k & !is.na(k.out$k),"k"] <- j
        k.out[match(a, k.out$genes),"k"] <- j
      }
    }
  }
  k.out$k <- as.numeric(as.factor(k.out$k))
  if(return_in_cytotable == TRUE){
    tab$group <- k.out[match(tab$destination, as.character(k.out$genes)),"k" ]
    return(tab)
  }else{return(k.out)}
}