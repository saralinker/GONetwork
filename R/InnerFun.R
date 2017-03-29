#' Inner function for an internal parallelizable cosine distance calculation
#'
#' .
#'
#' @param j row in M2
#' @param x x-value in cosine
#' 
#' @return cosine(x,y)
#'
#' @export


InnerFun <- function(j,x){
  y <- as.numeric(M2[j,])
  return(cosine(x,y))
}