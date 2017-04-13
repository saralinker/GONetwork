#' Find GO terms linking together a cluster of genes
#'
#' Returns a vector of GO terms associated with a group of genes
#'
#' @param M  Matrix of GO terms for each gene
#' @param genes a vector of gene(s) to retrieve the associated terms for.
#' @param proportion.shared if more than one gene is used as input, GO terms will only be returned if at least \textbf{proportion.shared} genes share the annotation. default = 1
#'
#' @return vector of GO terms
#'
#' @export
findTerms <- function(M,genes = NULL, proportion.shared = 1){
  
  terms <- M[genes,]
  if(length(genes) > 1){
    terms[terms > 0] <- 1
    terms <- colnames(terms[,(colSums(terms) / length(genes)) >= proportion.shared])
  }else{
    terms <- names(terms[terms > 0])
  }
  
  return(terms)
}