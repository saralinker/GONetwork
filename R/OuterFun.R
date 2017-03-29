OuterFun <- function(i){
  x <- as.numeric(M2[i,])
  tmp <- lapply(1:nrow(M2), InnerFun,x=x)
}