#' Calculate the GO term matrix using hierarchical knowledge
#'
#' .
#'
#' @param genes genes to query GO terms for 
#' @param species species to query from. Currently only supports human and mouse
#' @param preMinCol default = 0
#' @param preMinCol default = 0
#' @param maxthreads default = 3
#' 
#' @return Matrix of GO terms per gene
#'
#' @export



getGo <- function(genes, species = "mouse", preMinCol = 0, preMinRow = 0, maxthreads = 3){
  require("parallel", quietly = TRUE)
  load("R/sysdata.rda")
  source("R/FillM.R")
  
  
  if(species == "human"){
    gaf <- human_gaf
  }else if (species == "mouse"){
    gaf <- mouse_gaf
  }else{
    warning("species other than human and mouse are not yet supported")
  }
  
  res <- na.exclude(human_gaf[matches(genes, as.character(human_gaf$V3))$y,c(3,5)])
  colnames(res) <- c("external_gene_name","name_1006")
  ##########################
  ######## Create a {0,1} table of term belongingness per gene
  ##########################
  tab <- table(res$name_1006)
  #tab <- tab[-(c(which(names(tab) == "")))]
  ord <- names(tab)
  res <- res[res$name_1006 != "",]
  
  ####### 
  # Fill in Matrix (binary step); parallel step
  ####### 
  cl <- makeCluster(getOption("cl.cores",  maxthreads))
  clusterExport(cl=cl, varlist=c("genes","FillM","ord","res"))
  tmp <- parLapply(cl, X = genes, fun = FillM, ord = ord, res= res)
  M <- t(matrix(unlist(tmp),ncol = length(genes), nrow = length(ord)))
  rownames(M) <- genes
  colnames(M) <- ord
  stopCluster(cl)
  
  #######################
  ## Multiply each column by its respective number of parents
  #######################
  M2 <- M
  parents <- GO_space[match(colnames(M), as.character(GO_space$id)),"parents"]
  parents[is.na(parents)] <- 0
  M2 <- t(t(M2) * parents)
  Col <- colSums(M2)
  Row <- rowSums(M2)
  M2 <- M2[as.vector(Row) > preMinRow ,as.vector(Col) > preMinCol]
  
  return(M2)
}