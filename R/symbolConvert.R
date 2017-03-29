#' symbolConvert
#'
#' .
#'
#' @param x .
#' @param in.type .
#' @param out.type .
#' 
#' @return value
#'
#' @export
 


symbolConvert <- function(x, in.type, out.type){
  
  if(is.vector(x)==TRUE ) {
    x=as.data.frame(x)
  }
  #manually set objects for when you are debugging, ignore otherwise  
  #in.type="entrezgene"
  #out.type="hgnc_symbol"
  #x=adhdq
  
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  conversions <- getBM(filters= in.type, attributes= c(in.type,out.type),values=x,mart= mart)
  #conversions<-conversions[match(x,conversions[,1]),]
  
  row1<-as.character(x[,1])
  row2<-as.character(conversions[,1])
  nas<-is.na(match(row1,row2))
  notfound<-x[grep("TRUE",nas),]
  
  if(length(notfound)>=1) {
    warning(paste("the following entries were not found:",paste(notfound,sep=",",collapse=" , ")))
  }
  
  c<-conversions[match(row1,row2),]
  c<-c[is.na(c[,1])==FALSE,]
  c<-c[,2]
  
  return(c)
  
}
