#' Nested function for niche overlapping estimation with 'combina_niche3.R'.
#'
#' @description Select the couple of species for niche overlapping (index D) estimation.
#'
#' @param sp_chr a three dimensional matrix with the species chromatograms (alpha category by p environmental variables by species)
#' @param Thres_T an integer corresponding to the threshold of minimal abundance in a category for the niche breadth estimation
#'
#' @return index D
#' @export

niche_difer_sp<-function(sp_chr,Thres_T){

  n<-dim(sp_chr)
  y2<-matrix(nrow = n[length(n)],ncol = n[length(n)])

  if (length(n)==2) {
    for (i in 1:(n[length(n)]-1)) {
      for (j in ((i+1):n[length(n)])) {
        y2[i,j]<-niche_difer2(as.matrix(sp_chr[,i]),as.matrix(sp_chr[,j]),Thres_T)
      }
    }

  }else{
    for (i in 1:(n[length(n)]-1)) {
      for (j in ((i+1):n[length(n)])) {
        y2[i,j]<-niche_difer2(as.matrix(sp_chr[,,i]),as.matrix(sp_chr[,,j]),Thres_T)
      }
    }
  }

  return(y2)
}
