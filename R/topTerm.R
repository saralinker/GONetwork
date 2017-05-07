#' topTerm
#' determine the most descriptive term for each cluster
#'
#' 
#'
#' @param tab, output from AssignCluster
#' @param M, matrix containing hierarchy information for all genes
#' 
#' @return top term for each cluster assigned with AssignCluster
#'
#' @export

setClass("TermList",slots = c(sig.clusters = "data.frame", terms.in.matrix = "list") )

topTerm <- function(tab,M){
  if(sum(colnames(tab) == "group") == 0){
    warning("Need to assign clusters with AssignCluster() before determining the top term per cluster.")
    stop()
  }
  topterms <- lapply(X = unique(tab$group),FUN = loop_clusters, M = M, TAB = tab)
  names(topterms) <- unique(tab$group)
  sig <- unlist(lapply(names(topterms), FUN = function(x){nrow(topterms[[x]])}))
  sig <- data.frame(cluster = names(topterms[sig > 0]), term_count = sig[sig > 0])
  output.terms <- new(Class = "TermList", sig.clusters = sig, terms.in.matrix = topterms)
  return(output.terms)
}

