% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findTerms.R
\name{findTerms}
\alias{findTerms}
\title{Find GO terms linking together a cluster of genes}
\usage{
findTerms(M, genes = NULL, proportion.shared = 1)
}
\arguments{
\item{M}{Matrix of GO terms for each gene}

\item{genes}{a vector of gene(s) to retrieve the associated terms for.}

\item{proportion.shared}{if more than one gene is used as input, GO terms will only be returned if at least \textbf{proportion.shared} genes share the annotation. default = 1}
}
\value{
vector of GO terms
}
\description{
Returns a vector of GO terms associated with a group of genes
}
