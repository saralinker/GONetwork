#############################################
## GoNetwork_slinker, created 11/20/2015
## Purpose: Generate a network of how genes are linked based on GO terminology
## Requirements: reshape, biomaRt, lsa, GO.db
## User Input: genes, that's it!
## Other functions that can be added if you ask (just no time right now)
##    - adding in important variables so you can color/size/etc... based on FC or p-value in cytoscape
#############################################
# LOAD Public functions --------------------------------------------------------
requireNamespace("reshape", quietly = TRUE)
requireNamespace("biomaRt", quietly = TRUE)
requireNamespace("lsa", quietly = TRUE)
#requireNamespace("GO.db", quietly = TRUE)
#Initialize environment
gonet.env <- new.env()
# LOAD Homemade functions -------------------------------------------------
getGo <- function(genes, species = "mouse", preMinCol = 0, preMinRow = 0){
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
getGo_bin <- function(genes, species = "mouse", preMinCol = 0, preMinRow = 0){
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
OuterFun <- function(i){
  x <- as.numeric(M2[i,])
  tmp <- lapply(1:nrow(M2), InnerFun,x=x)
}
InnerFun <- function(j,x){
  y <- as.numeric(M2[j,])
  return(cosine(x,y))
}
GoTheDist <- function(M,Min = 4,minparents =   8, MinRow = 0, Threads = NULL){
  require(parallel)
  M2 <- M[,(as.vector(apply(M,2,function(x) sum(x > 0))) > Min)]
  parent <- as.numeric(apply(X=M2,MARGIN=2,max))
  M2 <- M2[,parent >= minparents]
  M2 <- M2[rowSums(M2) > MinRow,]
  if (!is.null(Threads)){
    cl <- makeCluster(getOption("cl.cores",  Threads))
    clusterExport(cl=cl, varlist=c("M2","InnerFun","cosine"))
    tmp <- parLapply(cl, X = c(1:nrow(M2)),fun = OuterFun)
    #tmp <- lapply( X = c(1:nrow(M2)),FUN = OuterFun)
    D2 <- matrix(unlist(tmp),nrow = nrow(M2), ncol = nrow(M2))
    rownames(D2) <- colnames(D2) <-  rownames(M2)
    stopCluster(cl)
  } else{
    D2 <- as.matrix(cosine(as.matrix(t(M2)))) #distance measure is cosine distance
  }
  D2[is.na(D2)] <- 0
  #Col <- colSums(D2)
  #Row <- rowSums(D2)
  rownames(D2) <- colnames(D2) <-  rownames(M2)
    return(D2)
}
Names <- function(x){
  paste(min(as.character(x[c(1,2)])),
        max(as.character(x[c(1,2)])),
        sep = "."
        )
}

cyto <- function(D,cutoff = 0.5){

  tmp <- melt(t(D))
  tmp2 <- tmp[tmp$value > cutoff & tmp$value !=1,]
  tmp2$value <- 1 - tmp2$value
  #tmp2$combname <- as.vector(unlist(apply(tmp2, 1, Names)))
  #tmp2$X3 <- do.call("rbind",strsplit(tmp2$combname,".",fixed =TRUE))[,1]
  #tmp2$X4 <- do.call("rbind",strsplit(tmp2$combname,".",fixed =TRUE))[,2]
  #tmp3 <- unique(tmp2[,c(3,5,6)])
  colnames(tmp2) <- c("origin","destination","distance")
  return(tmp2)
}
exactTerms <- function(x,term = "nuclear"){
  tmp <- x[grep(term,x,ignore.case=TRUE)]
  return(tmp)
}
findGenes <- function(M,term = "nuclear"){
  tmp <- M[,grep(term,colnames(M),ignore.case=TRUE)]
  if(!is.null(dim(tmp))){
    tmp2 <- tmp[rowSums(tmp) > 0,]
    return(rownames(tmp2))

  }else{
    #tmp2 <- tmp[rowSums(tmp) > 0,]
    #tmp2 <- as.list(tmp2)
    #names(tmp2) <- rownames(M[tmp> 0,])
    return(rownames(M[tmp> 0,]))
    }

}
findTerms <- function(M,geneA,geneB=NA,geneC=NA,geneD=NA){
  genecount <- 1
  a <- M[geneA,]
  a <- names(a[which(a > 0)])
  if(!is.na(geneB)){
    genecount  = genecount + 1
    b <- M[geneB,]
    a <- c(a,names(b[which(b > 0)]))
  }
  if(!is.na(geneC)){
    genecount  = genecount + 1
    b <- M[geneC,]
    a <- c(a,names(b[which(b > 0)]))
  }
  if(!is.na(geneD)){
    genecount  = genecount + 1
    b <- M[geneD,]
    a <- c(a,names(b[which(b > 0)]))
  }
  a <- table(a)
  return(names(a[a==genecount]))
}



