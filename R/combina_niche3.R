#' Estimate the degree of niche overlapping among species from their chromatograms
#'
#' @param sp_chr a matrix with the species chromatograms (categories by environmental variables by species). Outputs of 'chromato_env16.R'
#' @param T an integer corresponding to the threshold of minimal abundance in a category for niche breadth estimation
#'
#' @return results1 a matrix with the mean degree of niche overlapping. The first column displays the number of dimensions
#' considered simultaneously, columns 2 to 10 display the combinations of dimensions.
#' The last column displays index D associated with the combination of environmental dimensions.
#' D=0 when species niches are fully different and D=100 when species niches are identical;
#' the higher the number of dimensions, the lower the value of index D.
#' Only the combinations of environmental variables that minimise values of index D are displayed.
#' @return results2 a vector with the mean degree of niche overlapping when each dimension is considered alone.
#' @return y2 a matrix with the degree of niche overlapping species by species when all the dimensions are considered.
#' @export

combina_niche3<-function(sp_chr,T){

  n<-dim(sp_chr)
  v<-(1:n[2])
  nb<-0

  for (i in 1:n[2]) {
    temp<-t(utils::combn(v,i))
    ntemp<-dim(temp)
    nb<-ntemp[1]+nb
    rm(temp,ntemp)
  }


  point<-matrix(nrow = nb,ncol = 1)
  point2<-matrix(nrow = nb,ncol = n[2])
  y<-array(data=NA,dim = c(n[3],n[3],nb))
  moyenne<-matrix(nrow = nb,ncol = 1)

  cp<-0
  for (i in 1:n[2]) {
    temp<-t(utils::combn(v,i))
    ntemp<-dim(temp)

    for (j in 1:ntemp[1]) {
      cp<-cp+1
      y[,,cp]<-niche_difer_sp(sp_chr[,temp[j,],],T)
      point[cp]<-i
      point2[cp,1:i]<-temp[j,1:i]

      test<-y[,,cp]
      test2<-upper.tri(test,diag = FALSE)
      moyenne[cp]<-mean(test[test2==TRUE],na.rm = TRUE)
      rm(test,test2)
    }
    rm(temp,ntemp)
  }

  results1<-matrix(nrow = n[2],ncol = n[2]+2)

  for (i in 1:n[2]) {
    f<-which(point==i)
    temp1<-moyenne[f]
    mini<-min(temp1,na.rm = TRUE)
    pos<-which(temp1==mini)

    if (length(pos)>1 & length(unique(temp1))==1) {
      results1[i,]<-t(c(i,point2[f[pos[1]],],mini))
    }else if (length(pos)==1) {
      results1[i,]<-t(c(i,point2[f[pos],],mini))
    }else{
      warning('Problem with minimal overlapping values')
    }

    rm(f,temp1,mini,pos)
  }

  Index_D<-list(results1,y[,,cp],moyenne[1:n[2]])
  names(Index_D)<-c('combi_dim','sp_by_sp','dim_alone')
  return(Index_D)
}