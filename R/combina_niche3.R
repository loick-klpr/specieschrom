#' Estimate the degree of niche overlapping among species from their chromatograms.
#'
#' @description Function that returns the lowest degree of niche overlapping (index D) among all the species when 1 to p environmental variables
#' are considered. Also return the degree of niche overlapping (index D) among each couple of species when all the environmental
#' variables are considered simultaneously and the lowest degree of niche overlapping (index D) when each environmental variable is considered alone.
#'
#' @param sp_chr a matrix with the species chromatograms (alpha categories by p-environmental variables by species). Outputs of 'chromato_env16.R'
#' @param T an integer corresponding to the threshold of minimal abundance in a category for niche breadth estimation
#'
#' @return Return a list composed of three matrices:
#' @return combi_dim, which contains the mean degree of niche overlapping. In combi_dim, the first column displays the number of dimensions
#' considered simultaneously and the last column displays the index D associated with the combination of dimensions.
#' Columns between the first and the last display the combinations of dimensions considered.
#'
#' @return sp_by_sp, a matrix with the degree of niche overlapping (index D) species by species when all the dimensions are considered.
#'
#' @return dim_alone, a vector with the mean degree of niche overlapping (index D) when each dimension is considered alone.
#'
#' @return D=0 when species niches are fully different and D=100 when species niches are identical;
#' the higher the number of dimensions, the lower the value of index D.
#'
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
