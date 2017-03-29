#' Find GO terms linking together a cluster of genes
#'
#' Returns a vector of GO terms associated with a group of genes
#'
#' @param M 
#' @param genes
#' @param proportion.shared
#'
#'
#' @export
findTerms <- function(M,genes = NULL, proportion.shared = 1){
  
  terms <- M[genes,]
  if(length(genes) > 1){
    terms <- names(terms[(rowSums(terms) / length(genes)) >= proportion.shared,])
  }else{
    terms <- names(terms[terms > 0])
  }
  
  return(terms)
}