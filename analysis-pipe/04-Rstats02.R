#! /usr/bin/Rscript

### Load necessary libraries (all CRAN libraries can be acquired using install.packages)
library(compiler)
library(AnalyzeFMRI)
library(MASS)
library(abind)
library(fda)
library(fields)
library(Rcpp)
library(RcppEigen)
library(speedglm)
library(pracma)
library(lattice)
library(tclust)
library(dynamicTreeCut)

### Enable Just In Time Compiling, for ????MAYBE???? more speed.
enableJIT(3)

# SPLINE FITTING --------------------------------------

### Fit Splines to the time series in each of the slice files (ss loops over slices)
for (ss in 1:150) {
  ### Set where the slices are
  setwd("~/Dropbox/Collaborations/Jeremy Ullmann/Zebrafish_260216/for-hien/20150628_5dpf_H2BS_CON_LR_F11/03-chunk/NII")  
  ### Declare what files to load in
  ### Remove time points that do not look stationary
  file_list <- list.files()
  file_list <- file_list[grep('.nii',file_list)]
  file_number <- as.numeric(substr(file_list,7,12))
  
  ### These lines are kind of data dependent. I've determined what works for F11, F03, and F04.
  file_list <- file_list[which(file_number>2000 & file_number<5000)] # F11
  # file_list <- file_list[which(file_number>500 & file_number<=3300)] # F03 and F04
  file_number <- as.numeric(substr(file_list,7,12)) - 2000 # - 500 is F03 and F04 dependent, for F11 it would be - 2000.
  
  file_number.old <- file_number
  max_file <- max(file_number)
  file_number <- (file_number-1)*340 + ss
  
  ### Declare number of time slices left and initialize list
  N <- length(file_number)
  full_array <- list()
  
  ### Read in 75th Z slice
  for (ii in 1:N) {
    file_name <- substring(file_list[ii],1,12)
    nii_name <- paste(file_name,'.nii',sep='')
    NIFTI <- f.read.nifti.slice(nii_name,ss,1)
    full_array[[ii]] <- NIFTI
    print(c('reading',ss,ii,ii/N))
  }
  
  ### Combine all of the time slices
  full_array <- abind(full_array,along=3)
  
  ### Store data into matrix instead of array
  count <- 0
  full_mat <- matrix(NA,prod(170*70),N)
  for (ii in 1:170) {
    for (jj in 1:70) {
      count <- count + 1
      full_mat[count,] <- full_array[ii,jj,]
      print(c('storing',ss,count))
    }
  }
  
  ### Detrend data
  for (ii in 1:(170*70)) {
    full_mat[ii,] <- full_mat[ii,]-speedlm.fit(full_mat[ii,],cbind(1,file_number))$coefficients[2]*file_number
    # full_mat[ii,] <- detrend(full_mat[ii,])
    print(c('detrending',ss,ii))
  }
  
  ### Declare number of bases to use
  Basis_number <- 200
  Basis <- create.bspline.basis(c(0,(max_file-1)*340+150),
                                nbasis=Basis_number)
  BS <- eval.basis(file_number,Basis)

  ### Fit B-spline to all time series
  FD <- smooth.basisPar(file_number,t(full_mat),Basis)
  coeff_mat <- as.matrix(t(FD$fd$coefs))
  
  ### Set file location for where you want to save the results
  setwd("~/Dropbox/Collaborations/Jeremy Ullmann/Zebrafish_260216/for-hien/20150628_5dpf_H2BS_CON_LR_F11/03-chunk")  
  
  ### Save the results in the format coeff_mat_(SLICE NUMBER).rdata
  save(coeff_mat,file=paste('coeff_mat_',ss,'.rdata',sep=''))
}


# TRIMMED K-MEANS CLUSTERING ------------------------

### Set file location where the spline fitting results are saved
setwd("~/Dropbox/Collaborations/Jeremy Ullmann/Zebrafish_260216/for-hien/20150628_5dpf_H2BS_CON_LR_F11/03-chunk")  

### Load all of the coefficients and put them into a BIG matrix
load(paste('coeff_mat_',1,'.rdata',sep=''))
big_mat <- coeff_mat
for (ss in 2:150) {
  load(paste('coeff_mat_',ss,'.rdata',sep=''))
  big_mat <- rbind(big_mat,coeff_mat)
  print(ss)
}

### Function for computing BIC of K-means
tmeansBIC <- function(fit,data) {
  m = nrow(fit$centers)
  n = length(fit$cluster)
  k = ncol(fit$centers)
  PI = rep(1/k,k)
  log_densities = -as.matrix(pdist2(data,t(fit$center)))^2/2 - (m/2)*log(2*pi) - log(m)/2
  inner = sweep(log_densities,2,log(PI),'+')
  max_val = apply(inner,1,max)
  inner_minus = sweep(inner,1,max_val,'-')
  log_like_comps = max_val+log(rowSums(exp(inner_minus)))
  log_like = sum(log_like_comps)
  BIC = -2*log_like + log(n)*(m*k + k - 1)
  BIC
}

### Scale the coeffient matrix for K-means
scale_mat <- scale(big_mat)
# Remove big_mat for memory saving
rm(big_mat)

### Set a location for where BIC results should be saved
setwd("~/Dropbox/Collaborations/Jeremy Ullmann/Zebrafish_260216/for-hien/20150628_5dpf_H2BS_CON_LR_F11")  

### Compute BIC over a range of K (here 50--100)
BIC_val <- c()
for (kk in 50:100) {
  # Conduct timing while computing BIC values
  TIME <- proc.time()
  BIC_val[kk] <- tmeansBIC(tkmeans(x=scale_mat,k=kk,alpha=.9,nstart=3,iter.max=10),scale_mat)
  print(proc.time() - TIME)
  print(c(kk,BIC_val[kk]))
  
  # Save the results
  save(BIC_val,file='BIC_values.rdata')
}

### Get the optimal K and computer clustering under optimal K
comp <- which.min(BIC_val)

### Cluster using the optimal value for K
clustering <- tkmeans(x=scale_mat,k=comp,alpha=.9,nstart=5,iter.max=10)

### Function for allocating observations to cluster from a tkmeans clustering
tmeansClust <- function(fit,data) {
  apply(as.matrix(pdist2(data,t(fit$center))),1,which.min)
}

### Get a clustering using the tkmeans clustering form variable "clustering"
clustering_cluster <- tmeansClust(clustering,scale_mat)

### Produce Volume with cluster labels coordinates are (z,x,y)
image_hold <- array(NA,c(150,170,70))
count <- 0
for (ss in 1:150) {
  for (ii in 1:170) {
    for (jj in 1:70) {
      count <- count + 1
      image_hold[ss,ii,jj] <- clustering_cluster[count]
    }
  } 
}
# Plot the volume at the 75th slice
image.plot(image_hold[75,,],col=tim.colors(comp))

### Obtain the cluster mean time series
# Reload Graphing Parameters (These values are specific to F03 and F04)
file_number <- file_number.old
max_file <- max(file_number)
file_number <- (file_number-1)*340 + ss
Basis_number <- 200
Basis <- create.bspline.basis(c(0,(max_file-1)*340+150),
                              nbasis=Basis_number)
BS <- eval.basis(file_number,Basis)
# Compute the mean time series
centers <- t(clustering$centers)
centers <- sweep(centers,2,attributes(scale_mat)$'scaled:scale','*')
centers <- sweep(centers,2,attributes(scale_mat)$'scaled:center','+')
pred_range <- eval.basis(seq(1,((max_file-1)*340+150),1000),Basis)
# Each row is a mean time series over the pred_range values
PRED <- matrix(NA,dim(centers)[1],dim(pred_range)[1])
for (ii in 1:dim(centers)[1]) {
  PRED[ii,] <- apply(pred_range,1,function(x) {sum(x*centers[ii,])})
}
# The time average values of each series
MEAN_PRED <- rowMeans(PRED)


### Make a set of functions for evaluating splines and convolutions of splines
Spline_function <- function(x,cc) {sum(eval.basis(x,Basis)*centers[cc,])}
S1 <- function(x) Spline_function(x,1) 
S2 <- function(x) Spline_function(x,2) 
S_Prod <- function(x) S1(x)*S2(x)
# INTEGRAL <- integrate(Vectorize(S_Prod),1,(max_file-1)*340+150)$value

### Compute the Covariance of each mean function 
COVAR <- c()
for (cc in 1:comp) {
  S1 <- function(x) Spline_function(x,cc) 
  S2 <- function(x) Spline_function(x,cc) 
  S_Prod <- function(x) S1(x)*S2(x)
  INTEGRAL <- quadv(S_Prod,1,(max_file-1)*340+150)$Q
  COVAR[cc] <- INTEGRAL
}

### Compute the Correlation between pairs of mean functions
CORR <- matrix(NA,comp,comp)
for (c1 in 1:comp) {
  for (c2 in c1:comp) {
    S1 <- function(x) Spline_function(x,c1) 
    S2 <- function(x) Spline_function(x,c2) 
    S_Prod <- function(x) S1(x)*S2(x)
    INTEGRAL <- quadv(S_Prod,1,(max_file-1)*340+150)$Q
    CORR[c1,c2] <- INTEGRAL/sqrt(COVAR[c1]*COVAR[c2])
  }
}
# Plot Correlation matrix
image.plot(CORR)


### Create a Mean value image (using the predicted mean values from MEAN_PRED)
which_clust <- 6
mean_image <- array(NA,c(150,170,70))
count <- 0
for (ss in 1:150)
{
  for (ii in 1:170) {
    for (jj in 1:70) {
      count <- count + 1
      mean_image[ss,ii,jj] <- MEAN_PRED[clustering_cluster[count]]
    }
  }
}
# Plot the mean image at slice 90
image.plot(1:170,1:70,mean_image[90,,],col=grey.colors(100,0,1))


# Graphing --------------------------------------------
# Make New directory and set to new directory for graphs
dir.create('Graphs_and_results')
setwd('./Graphs_and_results')

### First set of graphs
# Graph clustering on every 5th slice
for (slice in seq(5,150,by=5)) {
  pdf(paste('Clustering_Slice_',slice,'.pdf',sep=''),paper='a4r')
  # pdf(paste('Clustering_Slice_',slice,'.pdf',sep=''),paper='a4r')
  ### Plot Clustering Image
  image_hold <- array(NA,c(150,170,70))
  count <- 0
  for (ss in 1:150) {
    for (ii in 1:170) {
      for (jj in 1:70) {
        count <- count + 1
        image_hold[ss,ii,jj] <- clustering_cluster[count]
      }
    } 
  }
  
  image.plot(image_hold[slice,,],col=tim.colors(comp))
  dev.off()
}

# Graph Mean Slices on every 5th slice
for (slice in seq(5,150,by=5)) {
  pdf(paste('Mean_Slice_',slice,'.pdf',sep=''),paper='a4r')
  # jpeg(paste('Mean_Slice_',slice,'.jpeg',sep=''))
  ### Plot Mean Image
  mean_image <- array(NA,c(150,170,70))
  count <- 0
  for (ss in 1:150)
  {
    for (ii in 1:170) {
      for (jj in 1:70) {
        count <- count + 1
        mean_image[ss,ii,jj] <- MEAN_PRED[clustering_cluster[count]]
      }
    }
  }
  image.plot(1:170,1:70,mean_image[slice,,],col=grey.colors(100,0,1))
  dev.off()
  
}

# Graph the location of clusters on slice 75
for (cc in 1:comp) {
  pdf(paste('Location_on_Slice_75_Cluster_',cc,'.pdf',sep=''),paper='a4r')
  # jpeg(paste('Location_on_Slice_75_Cluster_',cc,'.jpeg',sep=''))
  image.plot(1:170,1:70,mean_image[75,,],col=grey.colors(100,0,1))
  count <- 0
  for (ss in 1:150) {
    for (ii in 1:170) {
      for (jj in 1:70) {
        count <- count + 1
        if (clustering_cluster[count]==cc & ss==75) {
          points(ii,jj,pch=15,cex=1,col='green')
        }
      }
    }  
  }
  dev.off()
}

# Graph the location of clusters on slice 50
for (cc in 1:comp) {
  pdf(paste('Location_on_Slice_50_Cluster_',cc,'.pdf',sep=''),paper='a4r')
  # jpeg(paste('Location_on_Slice_75_Cluster_',cc,'.jpeg',sep=''))
  image.plot(1:170,1:70,mean_image[50,,],col=grey.colors(100,0,1))
  count <- 0
  for (ss in 1:150) {
    for (ii in 1:170) {
      for (jj in 1:70) {
        count <- count + 1
        if (clustering_cluster[count]==cc & ss==50) {
          points(ii,jj,pch=15,cex=1,col='green')
        }
      }
    }  
  }
  dev.off()
}

# Graph the location of clusters on slice 100
for (cc in 1:comp) {
  pdf(paste('Location_on_Slice_100_Cluster_',cc,'.pdf',sep=''),paper='a4r')
  # jpeg(paste('Location_on_Slice_75_Cluster_',cc,'.jpeg',sep=''))
  image.plot(1:170,1:70,mean_image[100,,],col=grey.colors(100,0,1))
  count <- 0
  for (ss in 1:150) {
    for (ii in 1:170) {
      for (jj in 1:70) {
        count <- count + 1
        if (clustering_cluster[count]==cc & ss==100) {
          points(ii,jj,pch=15,cex=1,col='green')
        }
      }
    }  
  }
  dev.off()
}

# Graph the Mean functions
for (cc in 1:comp) {
  pdf(paste('Cluster_Mean_Function_',cc,'.pdf',sep=''),paper='a4r')
  # pdf(paste('Cluster_Mean_Function_',cc,'.pdf',sep=''),paper='a4r')
  plot(seq(1,((max_file-1)*340+150),1000),PRED[cc,],type='l',xlab='Time',ylab='signal',main=cc)
  dev.off()
}

# Plot the Correlation matrix
pdf('Correlation_matrix.pdf',paper='a4r')
image.plot(1:comp,1:comp,CORR)
dev.off()


# Hierachical Clustering ------------------------------
# Make a distance metric
DIST <- as.dist(1-t(CORR))
# Conduct hierachical clustering
HCLUST <- hclust(DIST)

# Make tree cut using dynamicTreeCut
CUT <- cutreeDynamicTree(HCLUST,minModuleSize=1)

# Plot dendrogram
pdf('Dendrogram_clusters.pdf',width=30,height=10)
plot(HCLUST)
rect.hclust(HCLUST,max(CUT),border=rainbow(max(CUT)))
dev.off()

### Plot the 10 HCLust Cluster Means
for (ii in 1:max(CUT)) {
  pdf(paste('HCLUST_',ii,'.pdf',sep=''),paper='a4r')
  plot(c(1,((max_file-1)*340+150)),c(min(PRED),max(PRED)),
       type='n',xlab='Time',ylab='signal',main=ii)
  for (ss in 1:comp) {
    if (CUT[ss]==ii) {
      lines(seq(1,((max_file-1)*340+150),1000),PRED[ss,],col=tim.colors(comp)[ss])
    }
  }
  dev.off()
}

# Plot Frequency Histogram
pdf('Frequency_of_clusters.pdf',paper='a4r')
plot(table(clustering_cluster),xlab='cluster',ylab='Frequency')
dev.off()

# Plot K means by HCLUST reference
pdf('Cluster_by_HCLUST.pdf',paper='a4r')
plot(CUT,col=tim.colors(86),xlim=c(0,comp+1),ylim=c(0,8),xlab='Cluster',ylab='HCLUST')
grid()
dev.off()
