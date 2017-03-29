#' Outer function for an internal parallelizable cosine distance calculation
#'
#' .
#'
#' @param j row in M2
#' @param x x-value in cosine
#' 
#' @return cosine(x,y)
#'
#' @export


OuterFun <- function(i){
  x <- as.numeric(M2[i,])
  tmp <- lapply(1:nrow(M2), InnerFun,x=x)
}