#' loop_clusters
#' loop through each cluster to test proportions
#'
#' 
#'
#' @param cluster,  cluster id
#' @param prop, proportion of genes required to share the term; default = 0.3
#' 
#' @return result of proportions test for all terms in a clusters
#'
#' @export

loop_clusters <- function(CLUSTER, prop = 0.1,M = M, TAB  = tab, return.max = 5){
  #find all terms associated with group
  c.genes <- genes_in_cluster(CLUSTER, TAB)
  ft <- findTerms(M, genes = c.genes,proportion.shared = prop)
  if(return.max < length(ft)){
    ft.ancestor <- lapply(ft, FUN = function(x){max(M[,x])})
    #proportion
    terms <- M[c.genes,ft] 
    terms[terms > 0] <- 1
    ft.prop <- (colSums(terms) / length(genes)) 
    #score
    ft.score <- length(ft) - rank(as.numeric(ft.ancestor) * ft.prop)
    ft <- ft[ft.score <= return.max ]
  }
  #find genes associated with each term 
  if(length(ft) > 0){
    term_genes <- data.frame(do.call("rbind", lapply(X = ft, FUN = test_global_cluster_proportion, cluster = CLUSTER, M = M, tab  = TAB,cl.genes = c.genes)))
    dimnames(term_genes) <- list(ft, c("in.group.n","in.group.x","out.group.x","p","ancestor" ))
    term_genes <- na.exclude(term_genes)
    term_genes$padj <- p.adjust(term_genes$p)
    term_genes <- term_genes[term_genes$padj < 0.05,]
    term_genes <- term_genes[order(term_genes$p),]
  }else{
    term_genes <- NA
  }
  return(term_genes)
}

