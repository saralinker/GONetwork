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