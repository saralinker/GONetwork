#' Calculate the GO term matrix without hierarchical information
#'
#' .
#'
#' @param genes genes to query GO terms for 
#' @param species species to query from. Currently only supports human and mouse
#' @param preMinCol default = 0
#' @param preMinCol default = 0
#' 
#' @return Matrix of GO terms per gene
#'
#' @export


getGo_bin <- function(genes, species = "mouse", preMinCol = 0, preMinRow = 0){
  requireNamespace("biomaRt", quietly = TRUE)
  require(GO.db, quietly = TRUE)
  if(species == "human"){
    ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
  }else if (species == "mouse"){
    ensembl <- useMart("ensembl",dataset="mmusculus_gene_ensembl")
  }else{
    warning("species not interpretable, try using human or mouse")
  }
  filter <- c("external_gene_name") #TRY: listFilters(ensembl)
  attrib <- c("external_gene_name","name_1006","go_id","namespace_1003")
  res = getBM(attributes=attrib,filters=filter,values=genes,mart=ensembl)
  
  ##########################
  ######## Create a {0,1} table of term belongingness per gene
  ##########################
  tab <- table(res$name_1006)
  if (sum(names(tab) == "") > 0){
    tab <- tab[-(c(which(names(tab) == "")))]
  }
  M <- as.data.frame(matrix(data = 0,nrow=length(genes),ncol=length(tab)))
  rownames(M) <- genes
  colnames(M) <- names(tab)
  res <- res[res$name_1006 != "",]
  for (g in genes){
    group <- res[res$external_gene_name == g,"name_1006"]
    #group <- group[-c(which(group == ""))]
    M[g,group] <- 1
  }
  Col <- colSums(M)
  Row <- rowSums(M)
  M <- M[as.vector(Row) > preMinRow ,as.vector(Col) > preMinCol]
  return(M)
}