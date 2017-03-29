#' Calculate the GO term matrix using hierarchical knowledge
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

getGo <- function(genes, species = "mouse", preMinCol = 0, preMinRow = 0){
  require("biomaRt", quietly = TRUE)
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
  tab <- tab[-(c(which(names(tab) == "")))]
  M <- as.data.frame(matrix(data = 0,nrow=length(genes),ncol=length(tab)))
  rownames(M) <- genes
  colnames(M) <- names(tab)
  res <- res[res$name_1006 != "",]
  for (g in genes){
    group <- res[res$external_gene_name == g,"name_1006"]
    #group <- group[-c(which(group == ""))]
    M[g,group] <- 1
  }
  ##########################
  ######## Create a table of term parents and children
  ##########################
  terms <- data.frame(id = res$go_id , space = res$namespace_1003, name = res$name_1006)
  terms <- unique(terms[terms$id != "",])
  #BP
  bp <- terms[terms$space == "biological_process",]
  bp2 <- as.vector(lapply(as.character(bp$id),function(x) (length(GOBPANCESTOR[[x]]) )))
  bp$parents <- bp2
  #CC
  cc <- terms[terms$space == "cellular_component",]
  cc2 <- as.vector(lapply(as.character(cc$id),function(x) (length(GOCCANCESTOR[[x]]) )))
  cc$parents <- cc2
  #MF
  mf <- terms[terms$space == "molecular_function",]
  mf2 <- as.vector(lapply(as.character(mf$id),function(x) (length(GOMFANCESTOR[[x]]) )))
  mf$parents <- mf2
  
  terms <- as.data.frame(rbind(bp,cc,mf))
  #######################
  ## Multiply each column by its respective number of parents
  #######################
  M2 <- M
  for (i in 1:ncol(M)){
    M2[,colnames(M)[i]] <- M[,colnames(M)[i]] * unlist(terms[terms$name == colnames(M)[i],"parents"])
  }
  M2[is.na(M2)] <- 0
  Col <- colSums(M2)
  Row <- rowSums(M2)
  M2 <- M2[as.vector(Row) > preMinRow ,as.vector(Col) > preMinCol]
  
  return(M2)
}