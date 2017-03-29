#' Finds exact GO term annotation associated with a general description
#'
#' Akin to grep for the GO term annotations
#'
#' @param M matrix of GO terms for each gene 
#' @param term the general term that is being searched
#' 
#' @return vector of exact terms that match the given term (ie: grep)
#'
#' @export


exactTerms <- function(M,term = "nuclear"){
  x <- colnames(M)
  tmp <- x[grep(term,x,ignore.case=TRUE)]
  return(tmp)
}