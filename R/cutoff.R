#' cutoff
#'
#' .
#'
#' @param data .
#' @param method .
#' @param stringency .
#' @param Min .
#' @param minparents .
#' 
#' @return value
#'
#' @export

cutoff<- function(data,method,stringency,Min,minparents){
  requireNamespace("modes", quietly = TRUE)
  d.initial<- GoTheDist(data,method=paste(method),Min=Min,minparents=minparents)
  
  q<-quantile(d.initial,c(.05,.25,.5,.75,.95))
  
  if(method=="cosine"){
    qq<-match(T,q!=0)
    cut<- distandcyto(data, method=paste(method),low.cutoff=q[qq],high.cutoff=Inf,Min=Min,minparents=minparents )
  }
  else{
    cut<- distandcyto(data, method=paste(method),low.cutoff=-Inf,high.cutoff=q[3],Min=Min,minparents=minparents )
  }
  
  dist<- cut[,3]
  peaks<- amps(dist)$Peaks
  sec_peak<- peaks[2,1]
  s_t_midpoint<- mean(c(peaks[3,1],peaks[2,1]))
  third_peak<-peaks[3,1]
  three_four_midpoint<- mean(c(peaks[3,1],peaks[4,1]))
  fourth_peak<-peaks[4,1]
  if(stringency==2){
    return(sec_peak)
  }
  if(stringency==2.5){
    return(s_t_midpoint)
  }
  if(stringency==3){
    return(third_peak)
  }
  
  if(stringency==3.5){
    return(three_four_midpoint)
  }
  if(stringency==4){
    return(fourth_peak)
  }
  
}
