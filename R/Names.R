Names <- function(x){
  paste(min(as.character(x[c(1,2)])),
        max(as.character(x[c(1,2)])),
        sep = "."
  )
}