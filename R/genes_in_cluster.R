#' genes_in_cluster
#' return a list of all genes assigned to a given cluster
#'
#' 
#'
#' @param cluster, the cluster number 
#' @param tab, cluster matrix
#' 
#' @return list of genes belonging to a given cluster
#'
#' @export
#' 

genes_in_cluster <- function(cluster = 1, tab = NULL){
  return(unique(as.character(tab[tab$group == cluster, "origin"])))
}