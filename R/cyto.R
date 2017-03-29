#' Convert distance matrix to a cytoscape formatted table
#'
#' .
#'
#' @param D matrix of gene distances based on GO terms
#' @param cutoff, minimum distance between nodes
#' 
#' @return matrix in a cytoscape-friendly format
#'
#' @export


cyto <- function(D,cutoff = 0.5){
  tmp <- melt(t(D))
  tmp2 <- tmp[tmp$value > cutoff & tmp$value !=1,]
  tmp2$value <- 1 - tmp2$value
  colnames(tmp2) <- c("origin","destination","distance")
  return(tmp2)
}