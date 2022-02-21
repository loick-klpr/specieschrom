#' Mean of the n highest values of a vector
#' Estimates the mean of the n% of the samples in x with the highest abundance
#'
#' @param x a vector with the values to average
#' @param n an integer corresponding to the percentage of values to use
#'
#' @return mean of the n% of the samples in x with the highest values
#' @export

nanmean4<-function(x,n){
  p<-length(x)

  p1<-ceiling((n*p[1])/100)

  ind<-sort(x,decreasing = TRUE)
  y<-mean(ind[1:p1],na.rm = TRUE)
}

