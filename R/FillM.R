#' FillM
#'
#' internal function required for getGo
#'
#' @param g gene name
#' 
#' @return .
#'
#' @export


FillM <- function(g,ord,res){
  a <- rep(0,length(ord))
  go <- as.character(res[res$external_gene_name == g,"name_1006"])
  group <- match(go,ord)
  a[group] <- 1
  return(a)
}