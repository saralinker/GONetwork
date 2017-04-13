#' calculate the relative abundance of each group within a cluster compared to chance
#'
#' 
#'
#' @param k.out, output from AssignClusters(return_in_cytotable = FALSE)
#' @param group, a data.frame dictating which genes are in which test group. The first column is the gene name and the second column is the group name. Currently this test only supports two groups. For example, genes that are up-regulated versus genes that are down-regulated
#' 
#' @return results of a fisher.test() for each cluster
#'
#' @export

testClusters <- function(k.out, group){
  k.test <- data.frame()
  a <- table(group[,2]) 
  if(length(a) > 2){
    warning("testClusters currently only supports a comparison of 2 groups")
  }else{
    for(i in unique(na.exclude(k.out$k))){
      k_genes <- as.character(k.out[k.out$k == i,"genes"])
      k_direction <- table((group[match(k_genes, as.character(group[,1])),2]))
      df <- data.frame(exp = as.numeric(a), ob = as.numeric(k_direction),row.names = names(a))
      t <- fisher.test(df)
      k.test[i,"k"] <- i
      k.test[i,names(a)[1]] <- as.numeric(k_direction[1])
      k.test[i,names(a)[2]] <- as.numeric(k_direction[2])
      k.test[i,"fisher.oddsratio"] <- as.numeric(t$estimate)
      k.test[i,"fisher.p"] <- as.numeric(t$p.value)
    }
    k.test$padj <- p.adjust(k.test$fisher.p)
  }
  return(k.test)
}



