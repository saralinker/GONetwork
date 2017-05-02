#' Convert distance matrix to a cytoscape formatted table
#'
#' .
#'
#' @param D, matrix of gene distances based on GO terms
#' @param cutoff, minimum distance between nodes
#' @param remove.dups, default = FALSE This will remove self-loops. This is an advantage for visulazing in cytoscape. If you continue on to calling groups make sure this is set to FALSE because cytoscape will only color nodes that are in the origin column, thereby missing the remaining nodes.
#' 
#' 
#' @return matrix in a cytoscape-friendly format
#'
#' @export


cyto <- function(D,  cutoff = 0.5, remove.dups = FALSE){
  require("reshape", quietly = TRUE)
  tmp <- melt(t(D))
  tmp2 <- tmp[tmp$value > cutoff & tmp$value !=1,]
  tmp2$value <- 1 - tmp2$value
  colnames(tmp2) <- c("origin","destination","distance")
  if(remove.dups == TRUE){
    tmp2$dup <-  0
    for(i in 1:nrow(tmp2)){
      a <- as.character(tmp2[i,1])
      b <- as.character(tmp2[i,2])
      dup <- which(tmp2[,1] == b & tmp2[,2] == a)
      if(tmp2[dup,"dup"] != "original"){
        tmp2[dup,"dup"] <- "duplicated"
        tmp2[i,"dup"] <- "original"
      }
    }
    tmp2 <- tmp2[tmp2$dup == "original",c(1:3)]
  }
  
  return(tmp2)
}
