#' calculate the relative abundance of each group within a cluster compared to chance
#'
#' 
#'
#' @param k.out, output from cyto
#' @param group, a data.frame dictating which genes are in which test group. The first column is the gene name and the second column is the group name. Currently this test only supports two groups. For example, genes that are up-regulated versus genes that are down-regulated
#' 
#' @return results of a fisher.test() for each cluster
#'
#' @export

testClusters <- function(tab, group){
  k.test <- data.frame()
  a <- table(group[, "group"])
  group[,"group"] <- as.factor(group[,"group"])
  if(length(unique(group[,"group"])) > 2){
    warning("Using more than 2 groups in proportion test. Testing the proportion of the first group...")
  }
  STOP = FALSE
  if(sum(colnames(group) == "genes") == 0){
    warning('A column of genes, named "genes" is required')
    STOP = TRUE
  }
  if(sum(colnames(group) == "group") == 0){
    warning('A column of group identifiers, named "group" is required')
    STOP = TRUE
  }
  if (STOP == TRUE){
    stop()
  }
  op <- options(warn = (-1))
  k.out <- data.frame()
  i <- 0
  for(g in unique(as.character(tab$origin))){
    i <- i + 1
    k.out[i,"genes"] <- g
    k.out[i,"k"] <- unique(tab[tab$origin == g, "group"])
  }
  
  for (i in unique(na.exclude(k.out$k))) {
    k_genes <- as.character(na.exclude(k.out[k.out$k == i, "genes"]))
    k_direction <- table((group[match(k_genes, as.character(group[,"genes"])), 2]))
    t <- prop.test(x = k_direction[1],n = sum(k_direction),p = a[1]/sum(a))
    k.test[i, "group"] <- i
    for(n in 1:length(names(a))){
      k.test[i, paste("obs.",names(a)[n],sep = "")] <- as.numeric(k_direction[n])
    }
    for(n in 1:length(names(a))){
      k.test[i,paste("exp.",names(a)[n],sep = "")] <- (as.numeric(a[n]) / sum(a)) * sum(k_direction)
    }
    k.test[i, "odds_ratio"] <- as.numeric(t$statistic)
    k.test[i, "prop_test_p"] <- as.numeric(t$p.value)
  }
  k.test$padj <- p.adjust(k.test$prop_test_p)
  k.test <- k.test[order(k.test$prop_test_p),]
  k.test$group <- factor(k.test$group, k.test$group[1:nrow(k.test)])
  options(warn = 0)
  return(k.test)
}



