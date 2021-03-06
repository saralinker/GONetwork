% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testClusters.R
\name{testClusters}
\alias{testClusters}
\title{calculate the relative abundance of each group within a cluster compared to chance}
\usage{
testClusters(tab, group)
}
\arguments{
\item{group, }{a data.frame dictating which genes are in which test group. The first column is the gene name and the second column is the group name. Currently this test only supports two groups. For example, genes that are up-regulated versus genes that are down-regulated}

\item{k.out, }{output from cyto}
}
\value{
results of a fisher.test() for each cluster
}
\description{
calculate the relative abundance of each group within a cluster compared to chance
}
