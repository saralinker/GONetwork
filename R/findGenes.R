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


findGenes <- function(M,tab, group = NULL, term = "nuclear"){
  if(sum(colnames(M) == term) == 1){
    tmp <- M[,term]
  }else{
    tmp <- M[,grep(term,colnames(M),ignore.case=TRUE)]  
  }
  if(!is.null(dim(tmp))){
    tmp2 <- tmp[rowSums(tmp) > 0,]
    df <- data.frame(gene = rownames(tmp2),group = rep(0, times = nrow(tmp2)))
  }else{
    df <- data.frame(gene = names(tmp[tmp > 0]), group = rep(0, length(names(tmp[tmp > 0]))))
  }
  for (g in as.character(unique(tab$origin))){
    df[df$gene == g,"group" ] <- unique(as.character(tab[tab$origin == g, "group"]))
  }
  df <- df[order(df$group, decreasing=TRUE),]
  if (!is.null(group)){
    df <- as.character(df[df$group == group,]$gene)
  }
  return(df)
}