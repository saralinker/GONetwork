#' FillM
#'
#' internal function required for getGo
#'
#' @param j row in M2
#' @param x x-value in cosine
#' 
#' @return cosine(x,y)
#'
#' @export


FillM <- function(g){
  x <- as.numeric(M2[i,])
  tmp <- lapply(1:nrow(M2), InnerFun,x=x)
}