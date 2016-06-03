#! /usr/bin/Rscript
#
### Load necessary libraries (all CRAN libraries can be acquired using install.packages)
##
#
#

library(compiler)
library(AnalyzeFMRI)
library(MESS)
library(abind)
library(fda)
library(fields)
library(Rcpp)
library(RcppEigen)
library(speedglm)
library(pracma)
library(lattice)

### Images are sampled at .2 per second
### Time offset 340*170*750 >>> 150*70*170

enableJIT(3)

### Working Directories containing NII versions of my files
setwd("/data/home/uqajank2/data/NIF053-ZFISH-SCAPE/SCAPE-MNC2/20150628_5dpf_H2BS_CON_LR_F11/03-chunk/")

### Declare what files to load in
### Remove time points that do not look stationary
file_list <- list.files()
file_list <- file_list[grep('.nii',file_list)]
file_number <- as.numeric(substr(file_list,7,12))
# file_list <- file_list[which(file_number>2000 & file_number<5000)] # F11
file_list <- file_list[which(file_number>500 & file_number<=3300)] # F03 and F04
file_number <- as.numeric(substr(file_list,7,12)) - 500

### Declare number of time slices left and initialize list
N <- length(file_number)
full_array <- list()

print(paste("Found ", N, " time slices"))

### Read in 75th Z slice
pb <- txtProgressBar(min = 0, max = N, style = 3)
for (ii in 1:N) {
  file_name <- substring(file_list[ii],1,12)
  nii_name <- paste(file_name,'.nii',sep='')
  NIFTI <- f.read.nifti.slice(nii_name,75,1)
  full_array[[ii]] <- NIFTI
  setTxtProgressBar(pb, ii)
  #print(c(ii,ii/N))
}

### Combine all of the time slices
full_array <- abind(full_array,along=3)

### Store data into matrix instead of array
print("Re-org to matrix from array")
count <- 0
full_mat <- matrix(NA,prod(170*70),N)
pb <- txtProgressBar(min = 0, max = 170*70, style = 3)
for (ii in 1:170) {
  for (jj in 1:70) {
    count <- count + 1
    full_mat[count,] <- full_array[ii,jj,]
    setTxtProgressBar(pb, count)
    #print(count)
  }
}

### Detrend data
print("Detrend data")
pb <- txtProgressBar(min = 0, max = 170*70, style = 3)
for (ii in 1:(170*70)) {
  full_mat[ii,] <- full_mat[ii,]-speedlm.fit(full_mat[ii,],cbind(1,file_number))$coefficients[2]*file_number
  # full_mat[ii,] <- detrend(full_mat[ii,])
  setTxtProgressBar(pb, ii)
  #print(ii)
}

### Declare number of bases to use
Basis_number <- 200

### Generate Fourier Basis/B-Spline
Basis <- create.bspline.basis(c(0,max(file_number)), nbasis=Basis_number)
BS <- eval.basis(file_number,Basis)

### Fit B-spline to all time series
print("Bspline fit data")
coeff_mat <- matrix(NA,170*70,Basis_number)
count <- 0
pb <- txtProgressBar(min = 0, max = 170*70, style = 3)
for (ii in 1:(170*70)) {
  LM <- fastLmPure(BS,full_mat[ii,])
  # LM <- speedlm.fit(full_mat[ii,],BS)
  coeff_mat[ii,] <- LM$coefficients
  setTxtProgressBar(pb, ii)
  #print(ii)
}

### Function for computing BIC of K-means
kmeansBIC <- function(fit,data) {
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  PI = table(fit$cluster)/n
  log_densities = -as.matrix(pdist2(data,fit$center))^2/2 - (m/2)*log(2*pi) - log(m)/2
  inner = sweep(log_densities,2,log(PI),'+')
  max_val = apply(inner,1,max)
  inner_minus = sweep(inner,1,max_val,'-')
  log_like_comps = max_val+log(rowSums(exp(inner_minus)))
  log_like = sum(log_like_comps)
  BIC = -2*log_like + log(n)*(m*k + k - 1)
  BIC
}

### Scale the coeffient matrix for K-means
scale_mat <- scale(coeff_mat)

### Compute BIC over a range of K
BIC_val <- c()
for (kk in 10:20) {
  BIC_val[kk] <- kmeansBIC(kmeans(scale_mat,kk,100,10),scale_mat)
  print(c(kk,BIC_val[kk]))
}

### Get the optimal K and computer clustering under optimal K
comp <- which.min(BIC_val)
clustering <- kmeans(scale_mat,comp,100,20)

### Plot Clustering Image
image_hold <- matrix(NA,170,70)
count <- 0
for (ii in 1:170) {
  for (jj in 1:70) {
    count <- count + 1
    image_hold[ii,jj] <- clustering$cluster[count]
  }
}
image.plot(image_hold,col=tim.colors(comp))

### Plot first 3 principal components
princip <- princomp(scale_mat,cor=T)
plot(as.data.frame(princip$scores[,c(1:3)]),col=tim.colors(comp)[clustering$cluster])

### Plot the Basis function coefficients for each cluster
plot(c(0,Basis_number),c(min(clustering$centers),max(clustering$centers)),type='n')
for (ii in 1:comp) {
  lines(1:Basis_number,clustering$centers[ii,],col=tim.colors(comp)[ii])
}

### Plot a specific cluster overlayed on a mean image
which_clust <- 6
mean_image <- matrix(NA,170,70)
count <- 0
for (ii in 1:170) {
  for (jj in 1:70) {
    count <- count + 1
    mean_image[ii,jj] <- mean(full_mat[count,])
  }
}
image.plot(1:170,1:70,mean_image,col=heat.colors(100))
count <- 0
for (ii in 1:170) {
  for (jj in 1:70) {
    count <- count + 1
    if (clustering$cluster[count]==which_clust) {
      points(ii,jj,pch=15,cex=1,col='black')
    }
  }
}

### Obtain the cluster mean time series
centers <- clustering$centers
centers <- sweep(centers,2,attributes(scale_mat)$'scaled:scale','*')
centers <- sweep(centers,2,attributes(scale_mat)$'scaled:center','+')
pred_range <- eval.basis(1:max(file_number),Basis)
PRED <- matrix(NA,dim(centers)[1],dim(pred_range)[1])
for (ii in 1:dim(centers)[1]) {
  PRED[ii,] <- apply(pred_range,1,function(x) {sum(x*centers[ii,])})
}


### Plot any time series with its optimal cluster mean
number <- 10
clustering$cluster[number]
plot(full_mat[number,]~file_number,type='l')
lines(PRED[clustering$cluster[number],],col='red')

### Plot dissimilarity matrix, based on correlation
dissimilarity <- 1-cor(as.data.frame(t(PRED)),method='spearman')
distance <- as.dist(dissimilarity)
image.plot(1:comp,1:comp,cor(as.data.frame(t(PRED))[,hclust(distance)$order]))

### Plot joint trajectory of any two mean curves
xyplot(PRED[1,]~PRED[14,],col=tim.colors(max(file_number)),type='b',cex=2,pch=19)

### Plot period space spectrogram of any mean curve
SS <- spectrum(PRED[10,],c(11,11))
plot(log(0.2/SS$freq,2),SS$spec,type='l')
0.2/SS$freq[which.max(SS$spec)]
