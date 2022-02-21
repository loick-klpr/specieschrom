#' Display the multidimensional ecological niche of a species with a species chromatogram
#'
#' @param z a matrix with n samples by p environmental variables (i.e. the value of each environmental varibale in each sample)
#' @param y a vector with the abundance of a species in the n samples
#' @param alpha an integer corresponding to the number of category along each environmental variable
#' @param m an integer corresponding to the lowest number of samples needed in a category in order to have an estimation of the mean abundance
#' @param k an integer corresponding to the percentage of samples with the highest abundance values to use to estimate the mean abundance in a given category
#' @param order_smth an integer corresponding the order of the simple moving average applied along each niche dimension
#'
#' @return chr2 a matrix corresponding to the species chromatogram (alpha categories by p environmental variables)
#' @export

chromato_env16<-function(z,y,alpha,m,k,order_smth){
  n<-dim(z)
  xst<-matrix(nrow = n[1],ncol = n[2])

  for (i in 1:n[2]){
    xst[,i]<-(z[,i]-min(z[,i],na.rm = TRUE))/(max(z[,i],na.rm = TRUE)-min(z[,i],na.rm = TRUE))
  }

  catego<-seq(0, 1, by = (1/alpha))
  z1<-catego[1:length(catego)-1]
  z2<-catego[2:length(catego)]

  chr1<-matrix(nrow = length(z1), ncol = n[2])
  nbval<-matrix(nrow = length(z1),ncol = n[2])

  for (i in 1:length(z1)) {
    for (j in 1:n[2]) {

      if (i==length(z1)){
        f<-which(xst[,j]>=z1[i] & xst[,j]<=z2[i])
      }else{
        f<-which(xst[,j]>=z1[i] & xst[,j]<z2[i])
      }

      if (length(xst[f,j])>=m){
        chr1[i,j]<-nanmean4(y[f],k)
        nbval[i,j]<-length(xst[f,j])
      }
      rm(f)
    }
  }

  chr2<-matrix(nrow = length(z1), ncol = n[2])

  for (i in 1:n[2]) {
    tps<-moymob1(chr1[,i],order_smth)
    chr2[,i]<-tps/max(tps,na.rm = TRUE)
  }

  y2<-as.data.frame(chr2)
  category<-factor(1:alpha)
  y3<-cbind(y2,category)
  y4<-reshape2::melt(y3,id.vars="category")

  chrom_plot<-ggplot2::ggplot(y4,ggplot2::aes_string(x='variable', y='category', colour=0, fill='value')) + ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colours=colorRamps::blue2green2red(400)) + ggplot2::scale_color_gradientn(colours="#000000") +
    ggplot2::guides(colour="none") +
    ggplot2::labs(fill = "Standardised\nabundance") + ggplot2::xlab("Environmental variables")
  print(chrom_plot)

  return(chr2)
}

