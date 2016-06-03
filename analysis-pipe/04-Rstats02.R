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
library(NbClust)
library(capushe)

### Enable Just In Time Compiling, for ????MAYBE???? more speed.
enableJIT(3)

### Parse CMD Line Arguments
args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]

### These lines are kind of data dependent. I've determined what works for F11, F03, and F04.
## FROM and TO determine the start and end points of the trimming
FROM <- as.numeric(args[2])
TO <- as.numeric(args[3])

mask_fn <- args[4]

print(c('1indir', indir))
print(c('2FROM', FROM))
print(c('3TO', TO))
print(c('4mask', mask_fn))

#
#X_SIZE <- 70
#Y_SIZE <- 170
#Z_SIZE <- 150

# new crop sizes BAD post-doc... Swap X and Y later...
X_SIZE <- 130
Y_SIZE <- 640
Z_SIZE <- 300


### Set working directory for computation
# BAD post-doc, do not do this. (this made things fail later)
# see changes to list.files below also
# setwd(indir)

# SPLINE FITTING --------------------------------------

### Fit Splines to the time series in each of the slice files (ss loops over slices)
for (ss in 1:Z_SIZE) {
  print(c('Doing ', ss))
  ### Get already existing outputs
  Completed_output <- list.files(path=indir, pattern='coeff_mat', full.names=FALSE)
  #Completed_output <- Completed_output[grep('coeff_mat',Completed_output)]
  ### Skip things already in Completed_output

  if (paste('coeff_mat_',ss,'.rdata',sep='')%in%Completed_output){
    print(paste('coeff_mat_', ss, '.rdata', ' exists, loading', sep=''))
  }
  else{
    ### Declare what files to load in
    ### Remove time points that do not look stationary
    file_list <- list.files(path=indir, pattern='R-out-.*.nii', full.names=FALSE)
    #file_list <- file_list[grep('.nii',file_list)]
    file_number <- as.numeric(substr(file_list,7,12))

    file_list <- file_list[which(file_number>FROM & file_number<=TO)]
    file_number <- as.numeric(substr(file_list,7,12)) - FROM

    file_number.old <- file_number
    max_file <- max(file_number)
    file_number <- (file_number-1)*Z_SIZE + ss

    ### Declare number of time slices left and initialize list
    N <- length(file_number)
    full_array <- list()

    ### Read in Z slices
    for (ii in 1:N) {
      file_name <- substring(file_list[ii],1,12)
      nii_name <- paste(indir, '/', file_name,'.nii',sep='')
      print(c('reading',nii_name, ss,ii,ii/N))
      NIFTI <- f.read.nifti.slice(nii_name,ss,1)
      full_array[[ii]] <- NIFTI
    }


    ### Combine all of the time slices
    full_array <- abind(full_array,along=3)


    ### Store data into matrix instead of array
    count <- 0
    full_mat <- matrix(NA,prod(Y_SIZE*X_SIZE),N)
    for (ii in 1:Y_SIZE) {
      for (jj in 1:X_SIZE) {
        count <- count + 1
        full_mat[count,] <- full_array[ii,jj,]
        print(c('storing',ss,count))
      }
    }


    ### Detrend data
    print('detrending')
    for (ii in 1:(Y_SIZE*X_SIZE)) {
      full_mat[ii,] <- full_mat[ii,]-speedlm.fit(full_mat[ii,],cbind(1,file_number))$coefficients[2]*file_number
      # full_mat[ii,] <- detrend(full_mat[ii,])
      # print(c('detrending',ss,ii))
    }


    ### Declare number of bases to use
    Basis_number <- 100
    Basis <- create.bspline.basis(c(0,(max_file-1)*Z_SIZE+Z_SIZE),
                                  nbasis=Basis_number)
    BS <- eval.basis(file_number,Basis)


    ### Fit B-spline to all time series
    print(c('Spline fitting',ss))
    FD <- smooth.basisPar(file_number,t(full_mat),Basis)
    coeff_mat <- as.matrix(t(FD$fd$coefs))


    ### Save the results in the format coeff_mat_(SLICE NUMBER).rdata
    print(paste('Saving as: ', indir, '/coeff_mat_', ss, '.rdata', sep=''))
    save(coeff_mat,file=paste(indir,'/coeff_mat_',ss,'.rdata',sep=''))
  }
}


# TRIMMED K-MEANS CLUSTERING ------------------------
### Load Mask File
print(paste("Loading mask ", mask_fn, sep=""))
MASK <- f.read.nifti.volume(mask_fn)
MASK <- MASK[,,,1]
### Dummy Mask
D_Mask <- array(NA,c(Z_SIZE,Y_SIZE,X_SIZE))
for (ss in 1:Z_SIZE) {
  for (ii in 1:Y_SIZE) {
    for (jj in 1:X_SIZE) {
      D_Mask[ss,ii,jj] <- MASK[ii,jj,ss]>=.99995
    }
  }
}
ssDM <- sum(D_Mask)
### Make big matrix to store all series voxels that are unmasked
Basis_number <- 100
big_mat <- matrix(NA,ssDM,Basis_number)
Count <- 0
for (ss in 1:Z_SIZE) {
  load(paste(indir,'/coeff_mat_',ss,'.rdata',sep=''))
  InCount <- 0
  for (ii in 1:Y_SIZE) {
    for (jj in 1:X_SIZE) {
      InCount <- InCount + 1
      if(D_Mask[ss,ii,jj]) {
        Count <- Count + 1
        big_mat[Count,] <- coeff_mat[InCount,]
      }
    }
  }
  print(c(ss))
}
### Scale the coeffient matrix for K-means
scale_mat <- scale(big_mat)
# Remove big_mat for memory saving
rm(big_mat)

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

### Compute BIC over a range of K (here 50--100)
BIC_val <- c()
TIME_STORE <- c()
for (kk in 2:50) {
  # Conduct timing while computing BIC values
  TIME <- proc.time()
  BIC_val[kk] <- tmeansBIC(tkmeans(x=scale_mat,k=kk,alpha=.9,nstart=1,iter.max=20),scale_mat)
  TIME_DUMMY <- proc.time() - TIME
  print(TIME_STORE)
  TIME_STORE[kk] <- TIME_DUMMY[1] + TIME_DUMMY[2]
  print(c(kk,BIC_val[kk]))

  # Save the results
  save(TIME_STORE,file='Time_store.rdata')
  save(BIC_val,file='BIC_values.rdata')
}

### Get the optimal K and computer clustering under optimal K
n <- dim(scale_mat)[1]
m <- Basis_number
neg_like <- BIC_val - log(n)*(m*(1:length(BIC_val)))
log_like <- neg_like/(2)
ave_log_like <- log_like/n
names_vec <- 1:length(BIC_val)
complexity_h <- shape_h <- (m*(1:length(BIC_val)))
SHDATA <- cbind(names_vec,shape_h,complexity_h,ave_log_like)
DD <- DDSE(SHDATA)
comp <- as.numeric(attributes(DD)$model)

### Cluster using the optimal value for K
clustering <- tkmeans(x=scale_mat,k=comp,alpha=.9,nstart=5,iter.max=20)

### Function for allocating observations to cluster from a tkmeans clustering
tmeansClust <- function(fit,data) {
  apply(as.matrix(pdist2(data,t(fit$center))),1,which.min)
}

### Get a clustering using the tkmeans clustering form variable "clustering"
clustering_cluster <- tmeansClust(clustering,scale_mat)

### Produce Volume with cluster labels coordinates are (z,x,y)
image_hold <- array(NA,c(Z_SIZE,Y_SIZE,X_SIZE))
Count <- 0
for (ss in 1:Z_SIZE) {
  for (ii in 1:Y_SIZE) {
    for (jj in 1:X_SIZE) {
      if (D_Mask[ss,ii,jj]) {
        Count <- Count + 1
        image_hold[ss,ii,jj] <- clustering_cluster[Count]
      }
    }
  }
}
image.plot(1:Y_SIZE,1:X_SIZE,image_hold[Z_SIZE,,])

### Obtain the cluster mean time series
# Reload Graphing Parameters (These values are specific to F03 and F04)
file_number <- file_number.old
max_file <- max(file_number)
file_number <- (file_number-1)*Z_SIZE + ss
Basis_number <- Basis_number
Basis <- create.bspline.basis(c(0,(max_file-1)*Z_SIZE+Z_SIZE),
                              nbasis=Basis_number)
BS <- eval.basis(file_number,Basis)
# Compute the mean time series
centers <- t(clustering$centers)
centers <- sweep(centers,2,attributes(scale_mat)$'scaled:scale','*')
centers <- sweep(centers,2,attributes(scale_mat)$'scaled:center','+')
pred_range <- eval.basis(seq(1,((max_file-1)*Z_SIZE+Z_SIZE),1000),Basis)
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
# INTEGRAL <- integrate(Vectorize(S_Prod),1,(max_file-1)*Z_SIZE+Z_SIZE)$value

### Compute the Covariance of each mean function
COVAR <- c()
for (cc in 1:comp) {
  S1 <- function(x) Spline_function(x,cc)
  S2 <- function(x) Spline_function(x,cc)
  S_Prod <- function(x) S1(x)*S2(x)
  INTEGRAL <- quadv(S_Prod,1,(max_file-1)*Z_SIZE+Z_SIZE)$Q
  COVAR[cc] <- INTEGRAL
}

### Compute the Correlation between pairs of mean functions
CORR <- matrix(NA,comp,comp)
for (c1 in 1:comp) {
  for (c2 in c1:comp) {
    S1 <- function(x) Spline_function(x,c1)
    S2 <- function(x) Spline_function(x,c2)
    S_Prod <- function(x) S1(x)*S2(x)
    QUAD <- quadv(S_Prod,1,(max_file-1)*Z_SIZE+Z_SIZE)
    INTEGRAL <- QUAD$Q
    INT_OLD <- INTEGRAL
    PREC <- QUAD$estim.prec
    CORR[c1,c2] <- INTEGRAL/sqrt(COVAR[c1]*COVAR[c2])
    while( xor(CORR[c1,c2] > 1, CORR[c1,c2] < -1) ) {
      if (CORR[c1,c2] > 1) {INTEGRAL <- INTEGRAL - PREC*INT_OLD}
      if (CORR[c1,c2] < -1) {INTEGRAL <- INTEGRAL + PREC*INT_OLD}
      CORR[c1,c2] <- INTEGRAL/sqrt(COVAR[c1]*COVAR[c2])
    }
  }
}



### Create a Mean value image (using the predicted mean values from MEAN_PRED)
which_clust <- 6
mean_image <- array(NA,c(Z_SIZE,Y_SIZE,X_SIZE))
count <- 0
for (ss in 1:Z_SIZE)
{
  for (ii in 1:Y_SIZE) {
    for (jj in 1:X_SIZE) {
      if (D_Mask[ss,ii,jj]) {
        count <- count + 1
        mean_image[ss,ii,jj] <- MEAN_PRED[clustering_cluster[count]]
      }
    }
  }
}

### Save All Intermediate Results
save(clustering_cluster,file=paste(indir,'/clustering.rdata',sep=''))
save(mean_image,file=paste(indir,'/mean_image.rdata',sep=''))
save(PRED,file=paste(indir,'/predicted_means.rdata',sep=''))
save(CORR,file=paste(indir,'/correlation_matrix.rdata',sep=''))
save(BIC_val,file=paste(indir,'/BIC_values.rdata',sep=''))
save(TIME_STORE,file=paste(indir,'/Time_Store.rdata',sep=''))

# Graphing --------------------------------------------

### First set of graphs
# Graph clustering on every 5th slice
for (slice in seq(5,Z_SIZE,by=5)) {
  pdf(paste(indir,'/Clustering_Slice_',slice,'.pdf',sep=''),paper='a4r')
  ### Plot Clustering Image
  image_hold <- array(NA,c(Z_SIZE,Y_SIZE,X_SIZE))
  count <- 0
  for (ss in 1:Z_SIZE) {
    for (ii in 1:Y_SIZE) {
      for (jj in 1:X_SIZE) {
        if (D_Mask[ss,ii,jj]) {
          count <- count + 1
          image_hold[ss,ii,jj] <- clustering_cluster[count]
        }
      }
    }
  }


  image.plot(image_hold[slice,,],col=tim.colors(comp))
  dev.off()
}

# Graph Mean Slices on every 5th slice
for (slice in seq(5,Z_SIZE,by=5)) {
  pdf(paste(indir,'/Mean_Slice_',slice,'.pdf',sep=''),paper='a4r')
  ### Plot Mean Image
  mean_image <- array(NA,c(Z_SIZE,Y_SIZE,X_SIZE))
  count <- 0
  for (ss in 1:Z_SIZE)
  {
    for (ii in 1:Y_SIZE) {
      for (jj in 1:X_SIZE) {
        if (D_Mask[ss,ii,jj]) {
          count <- count + 1
          mean_image[ss,ii,jj] <- MEAN_PRED[clustering_cluster[count]]
        }
      }
    }
  }
  image.plot(1:Y_SIZE,1:X_SIZE,mean_image[slice,,],col=grey.colors(100,0,1))
  dev.off()


}

# Graph the location of clusters on slice 75
for (cc in 1:comp) {
  pdf(paste(indir,'/Location_on_Slice_75_Cluster_',cc,'.pdf',sep=''),paper='a4r')
  image.plot(1:Y_SIZE,1:X_SIZE,mean_image[75,,],col=grey.colors(100,0,1))
  count <- 0
  for (ss in 1:Z_SIZE) {
    for (ii in 1:Y_SIZE) {
      for (jj in 1:X_SIZE) {
        if (D_Mask[ss,ii,jj]) {
          count <- count + 1
          if (clustering_cluster[count]==cc & ss==75) {
            points(ii,jj,pch=15,cex=1,col='green')
          }
        }
      }
    }
  }
  dev.off()
}

# Graph the location of clusters on slice 50
for (cc in 1:comp) {
  pdf(paste(indir,'/Location_on_Slice_50_Cluster_',cc,'.pdf',sep=''),paper='a4r')
  image.plot(1:Y_SIZE,1:X_SIZE,mean_image[50,,],col=grey.colors(100,0,1))
  count <- 0
  for (ss in 1:Z_SIZE) {
    for (ii in 1:Y_SIZE) {
      for (jj in 1:X_SIZE) {
        if (D_Mask[ss,ii,jj]) {
          count <- count + 1
          if (clustering_cluster[count]==cc & ss==50) {
            points(ii,jj,pch=15,cex=1,col='green')
          }
        }
      }
    }
  }
  dev.off()
}

# Graph the location of clusters on slice 100
for (cc in 1:comp) {
  pdf(paste(indir,'/Location_on_Slice_100_Cluster_',cc,'.pdf',sep=''),paper='a4r')
  image.plot(1:Y_SIZE,1:X_SIZE,mean_image[100,,],col=grey.colors(100,0,1))
  count <- 0
  for (ss in 1:Z_SIZE) {
    for (ii in 1:Y_SIZE) {
      for (jj in 1:X_SIZE) {
        if (D_Mask[ss,ii,jj]) {
          count <- count + 1
          if (clustering_cluster[count]==cc & ss==100) {
            points(ii,jj,pch=15,cex=1,col='green')
          }
        }
      }
    }
  }
  dev.off()
}

# Graph the Mean functions
for (cc in 1:comp) {
  pdf(paste(indir,'/Cluster_Mean_Function_',cc,'.pdf',sep=''),paper='a4r')
  plot(seq(1,((max_file-1)*Z_SIZE+Z_SIZE),1000),PRED[cc,],type='l',xlab='Time',ylab='signal',main=cc)
  dev.off()
}

# Hierachical Clustering ------------------------------
# Plot the Correlation matrix
pdf(paste(indir,'/Correlation_matrix.pdf',sep=''),paper='a4r')
image.plot(1:comp,1:comp,CORR)
dev.off()

# Make a distance metric
DIST <- as.dist(1-t(CORR))
HCLUST <- hclust(DIST,method='average')

# Make tree cut using Dunn Index
NB <- NbClust(diss=DIST,distance=NULL,method='average',index='silhouette',max.nc=ceiling(comp-2))
CUT <- NB$Best.partition

# Plot dendrogram
pdf(paste(indir,'/Dendrogram_clusters.pdf',sep=''),width=30,height=10)
plot(HCLUST,xlab='')
rect.hclust(HCLUST,max(CUT),border=rainbow(max(CUT)))
dev.off()

### Plot the 10 HCLust Cluster Means
for (ii in 1:max(CUT)) {
  pdf(paste(indir,'/HCLUST_',ii,'.pdf',sep=''),paper='a4r')
  plot(c(1,((max_file-1)*Z_SIZE+Z_SIZE)),c(min(PRED),max(PRED)),
       type='n',xlab='Time',ylab='signal',main=ii)
  for (ss in 1:comp) {
    if (CUT[ss]==ii) {
      lines(seq(1,((max_file-1)*Z_SIZE+Z_SIZE),1000),PRED[ss,],col=tim.colors(comp)[ss])
    }
  }
  dev.off()
}

# Plot Frequency Histogram
pdf(paste(indir,'/Frequency_of_clusters.pdf',sep=''),paper='a4r')
plot(table(clustering_cluster),xlab='cluster',ylab='Frequency')
dev.off()

# Plot K means by HCLUST reference
pdf(paste(indir,'/Cluster_by_HCLUST.pdf',sep=''),paper='a4r')
plot(CUT,col=tim.colors(comp),xlim=c(0,comp+1),ylim=c(0,8),xlab='Cluster',ylab='HCLUST')
grid()
dev.off()


# Graph clustering based on cuts on every 5th slice
for (slice in seq(5,Z_SIZE,by=5)) {
  pdf(paste(indir,'/CUT_Slice_',slice,'.pdf',sep=''),paper='a4r')
  ### Plot Clustering Image
  image_hold <- array(NA,c(Z_SIZE,Y_SIZE,X_SIZE))
  count <- 0
  for (ss in 1:Z_SIZE) {
    for (ii in 1:Y_SIZE) {
      for (jj in 1:X_SIZE) {
        if (D_Mask[ss,ii,jj]) {
          count <- count + 1
          image_hold[ss,ii,jj] <- CUT[clustering_cluster[count]]
        }
      }
    }
  }


  image.plot(image_hold[slice,,],col=tim.colors(max(CUT)))
  dev.off()
}
