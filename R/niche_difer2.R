#' Nested function for niche overlapping estimation with 'combina_niche3.R' and 'niche_differ_sp.R'.
#'
#' @description Estimate the niche overlapping (index D) between two species. Assume that the niche is rectangular.
#' Columns (i.e. environmental variables) in 'spe_chr1' and 'sp_chr2' have to be the same.
#'
#' @param sp_chr1 a matrix (niche categories by environmental variables) for species 1
#' @param sp_chr2 a matrix (niche categories by environmental variables) for species 2
#' @param Thres_T an integer corresponding to the threshold of minimal abundance in a category for the niche breadth estimation
#'
#' @return index D

niche_difer2<-function(sp_chr1,sp_chr2,Thres_T){

  if (Thres_T==0) {
    sp_chr1[sp_chr1>Thres_T]<-1
    sp_chr1[sp_chr1<=Thres_T]<-0

    sp_chr2[sp_chr2>Thres_T]<-1
    sp_chr2[sp_chr2<=Thres_T]<-0

  }else{
    sp_chr1[sp_chr1>=Thres_T]<-1
    sp_chr1[sp_chr1<Thres_T]<-0

    sp_chr2[sp_chr2>=Thres_T]<-1
    sp_chr2[sp_chr2<Thres_T]<-0
  }

  n<-dim(sp_chr1)
  for (i in 1:n[length(n)]) {
    f1<-which(sp_chr1[,i]==1)
    f2<-which(sp_chr2[,i]==1)
    sp_chr1[f1[1]:f1[length(f1)],i]<-1
    sp_chr2[f2[1]:f2[length(f2)],i]<-1
    rm(f1,f2)
  }

  y<-matrix(nrow = n[1],ncol = n[2])
  a1<-sp_chr1+sp_chr2
  f1<-which(a1==2)

  a2<-matrix(nrow = n[1],ncol = n[2])
  a2[f1]<-1

  z<-prod(colSums(a2,na.rm=TRUE),na.rm=TRUE)
  a3<-prod(colSums(sp_chr1,na.rm=TRUE),na.rm=TRUE)
  a4<-prod(colSums(sp_chr2,na.rm=TRUE),na.rm=TRUE)

  D<-(100*z)/(a3+a4-z)
  return(D)
}
