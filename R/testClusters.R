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
  op <- options(warn = (-1))
  k.test <- data.frame()
  a <- table(group[, 2])
  group[,2] <- as.factor(group[,2])
  
  for (i in unique(na.exclude(k.out$k))) {
    k_genes <- as.character(na.exclude(k.out[k.out$k == i, "genes"]))
    k_direction <- table((group[match(k_genes, as.character(group[,1])), 2]))
    propdif <- length(k_genes) / sum(a)
    df <- data.frame(exp = as.numeric(a) , ob = as.numeric(k_direction), 
                     row.names = names(a))
    t <- fisher.test(df)
    k.test[i, "k"] <- i
    for(n in 1:length(names(a))){
      k.test[i, names(a)[n]] <- as.numeric(k_direction[n])
    }
    k.test[i, "chi_statistic"] <- as.numeric(t$statistic)
    k.test[i, "chi_p"] <- as.numeric(t$p.value)
  }
  k.test$padj <- p.adjust(k.test$chi_p)
  options(op)
  return(k.test)
}



