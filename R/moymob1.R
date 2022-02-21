#' m order simple moving average
#'
#' @param x a matrix with in lines the observations and in column the variables
#' @param m an integer corresponding to the window size (2m+1)
#'
#' @return smoothed matrix
#' @export

moymob1<-function(x,m){
  n<-length(x)
  fenetre<-2*m+1
  y_2<-matrix(nrow = n,ncol = 1)

  for (i in (m+1):(n-m)) {
    y_2[i]<-mean(x[(i-m):(i+m)],na.rm = TRUE)
  }

  m2<-m

  while (m2!=1) {
    m2<-m2-1
    for (i in (m2+1):m) {
      y_2[i]<-mean(x[(i-m2):(i+m2)],na.rm = TRUE)
    }
  }

  y_2[1]<-mean(x[1:2])

  m2=m
  while (m2!=1) {
    m2<-m2-1
    for (i in (n-m+1):(n-m2)) {
      y_2[i]<-mean(x[(i-m2):(i+m2)],na.rm = TRUE)
    }
  }

  y_2[n]<-mean(x[(n-1):n])
  y_2<-y_2
}
