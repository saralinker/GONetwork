#' Find genes associated with a given GO term
#'
#' Returns a vector of genes associated with a given GO term
#'
#' @param M  Matrix of GO terms for each gene
#' @param term The GO term to call
#' 
#' @return vector of gene names
#'
#' @export


findGenes <- function(M,term = "nuclear"){
  tmp <- M[,grep(term,colnames(M),ignore.case=TRUE)]
  if(!is.null(dim(tmp))){
    tmp2 <- tmp[rowSums(tmp) > 0,]
    return(rownames(tmp2))
    
  }else{
    return(rownames(M[tmp > 0,]))
  }
}