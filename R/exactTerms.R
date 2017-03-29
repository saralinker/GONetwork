exactTerms <- function(x,term = "nuclear"){
  tmp <- x[grep(term,x,ignore.case=TRUE)]
  return(tmp)
}