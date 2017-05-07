#' test_global_cluster_proportion
#' test the likelihood of a term belonging to a cluster
#'
#' 
#'
#' @param t, term
#' 
#' @return result of proportions test for a given term
#'
#' @export


test_global_cluster_proportion <- function(M, tab, t, cluster,cl.genes = c.genes){
    op <- options(warn = (-1))
    fg <- findGenes(M = M[unique(as.character(tab$origin)),], tab = tab, term = t)
    fg.t <- table(fg$group)
    in.group.x <- as.numeric(fg.t[as.character(cluster)]) 
    in.group.n <- length(cl.genes)
    out.group.x <- sum(as.numeric(fg.t[-which(names(fg.t) == as.character(cluster))])) 
    out.group.n <- length(unique(tab$origin)) - length(cl.genes)
    a <- in.group.n
    b <- in.group.x
    d <- out.group.x
    e <- prop.test(x = c(in.group.x, out.group.x), n = c(in.group.n, out.group.n))
    if(e$estimate[1] > e$estimate[2]){
      tmp <- as.numeric(unique(M[,t]))
      f <- tmp[tmp>0]
    return(c(a,b,d,as.numeric(e$p.value),f))
    }else{return(c(NA,NA,NA,NA,NA))
      options(warn = 0)
    }
}