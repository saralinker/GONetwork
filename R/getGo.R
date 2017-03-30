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
  require("biomaRt", quietly = TRUE)
  require("parallel", quietly = TRUE)
  load("R/sysdata.rda")
  
  
  if(species == "human"){
    ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
  }else if (species == "mouse"){
    ensembl <- useMart("ensembl",dataset="mmusculus_gene_ensembl")
  }else{
    warning("species other than human and mouse are not yet supported")
  }
    
    filter <- c("external_gene_name") #TRY: listFilters(ensembl)
    attrib <- c("external_gene_name","name_1006","go_id","namespace_1003")
    res = getBM(attributes=attrib,filters=filter,values=genes,mart=ensembl)
  
  
  ##########################
  ######## Create a {0,1} table of term belongingness per gene
  ##########################
  tab <- table(res$name_1006)
  tab <- tab[-(c(which(names(tab) == "")))]
  ord <- names(tab)
  res <- res[res$name_1006 != "",]
  
  ####### 
  # Fill in Matrix (binary step); parallel step
  ####### 
  cl <- makeCluster(getOption("cl.cores",  maxthreads))
  clusterExport(cl=cl, varlist=c("genes","FillM","ord","res"))
  tmp <- parLapply(cl, X = genes, fun = FillM, ord = ord, res= res)
  M <- matrix(unlist(tmp),nrow = length(genes), ncol = length(ord))
  rownames(M) <- genes
  colnames(M) <- ord
  stopCluster(cl)
  
  
  ##########################
  ######## Create a table of term parents and children; not parallelized because the GO.db class is not subsettable
  ##########################
  terms <- data.frame(id = res$go_id , space = res$namespace_1003, name = res$name_1006)
  terms <- unique(terms[terms$id != "",])
  

  

  

  #BP
  bp <- terms[terms$space == "biological_process",]
  bp$parents <- bpanc[as.character(bp$id),"bp"]
  #CC
  cc <- terms[terms$space == "cellular_component",]
  cc$parents <- ccanc[as.character(cc$id),"bp"]
  #MF
  mf <- terms[terms$space == "molecular_function",]
  mf$parents <- mfanc[as.character(mf$id),"bp"]
  
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