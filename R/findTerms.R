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