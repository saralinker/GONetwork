\name{findTerms}
\alias{findTerms}
\title{Find GO terms linking together a cluster of genes}
\usage{
findTerms(M,genes = NULL, proportion.shared = 1)
}
\arguments{
\item{M}{Matrix of GO terms for each gene}
\item{genes}{a vector of gene(s) to retrieve the associated terms for.}
\item{proportion.shared}{if more than one gene is used as input, GO terms will only be returned if at least \textbf{proportion.shared} genes share the annotation. default = 1 }
}
\description{
If a single gene, all associated GO terms will be reported. If multiple genes, GO terms will only be returned if at least \textbf{proportion.shared} genes share the annotation.
}

\details{

The DESCRIPTION file:
\packageDESCRIPTION{GONetwork}
\packageIndices{GONetwork}
~~ An overview of how to use the package, including the ~~
~~ most important functions ~~
}
