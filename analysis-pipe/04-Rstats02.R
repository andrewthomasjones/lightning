### Load necessary libraries (all CRAN libraries can be acquired using install.packages)
library(compiler)
library(AnalyzeFMRI)
library(ggplot2)
library(reshape2)
library(MASS)
library(abind)
library(fda)
library(fields)
library(speedglm)
library(pracma)
library(lattice)
library(tclust)
#library(NbClust)
library(capushe)

#C++ stuff
library(Rcpp)
library(RcppEigen)
library(RcppArmadillo)
#needs to be c++ 11 because I used the tuple class
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
if(file.exists('./bin/tkmeans.cpp')){sourceCpp('./bin/tkmeans.cpp')} #for testing
if(file.exists('./tkmeans.cpp')){sourceCpp('./tkmeans.cpp')} #actual file structure

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

# #
# X_SIZE <- 70
# Y_SIZE <- 170
# Z_SIZE <- 150

# new crop sizes BAD post-doc... Swap X and Y later...
#this will be approximately 18.6 GB in RAM with 100 basis function params

X_SIZE <- 130
Y_SIZE <- 640
Z_SIZE <- 300


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
print(paste("Trying to allocate matrix of size (approx)", round(X_SIZE*Y_SIZE*Z_SIZE*Basis_number*8/(1024^3),2), " GB..."))
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
  print(paste("Loading slice", c(ss), "of", Z_SIZE))
}
### Scale the coeffient matrix for K-means
#big_mat <- scale(big_mat)
## Remove big_mat for memory saving
#rm(big_mat)

#scales inplace, returns means and sd for later use
mean_sd_from_unscaled<-scale_lowmem(big_mat)
save(mean_sd_from_unscaled, file = paste0(indir,"/scalingparams.rdata"))
print("Matrix successfully loaded and rescaled.")

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


#probably put a better check here
if(!file.exists(paste0(indir,'/clustering.rdata')) | !file.exists(paste0(indir,'/centers.rdata')) ){
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
 

  # ### Function for computing BIC of K-means
  # tmeansBIC <- function(fit,data) {
  #   m = nrow(fit$centers)
  #   n = length(fit$cluster)
  #   k = ncol(fit$centers)
  #   PI = rep(1/k,k)
  #   log_densities = -as.matrix(pdist2(data,t(fit$center)))^2/2 - (m/2)*log(2*pi) - log(m)/2
  #   inner = sweep(log_densities,2,log(PI),'+')
  #   max_val = apply(inner,1,max)
  #   inner_minus = sweep(inner,1,max_val,'-')
  #   log_like_comps = max_val+log(rowSums(exp(inner_minus)))
  #   log_like = sum(log_like_comps)
  #   BIC = -2*log_like + log(n)*(m*k + k - 1)
  #   BIC
  # }
  if(file.exists(paste0(indir,'/centers.rdata'))){load(file=paste0(indir,'/centers.rdata'))}
  
  if(!exists('clustering')){
    print("Doing clustering")
    ### Compute BIC over a range of K (here 50--100)
    BIC_val <- c()
    TIME_STORE <- c()
    
    max_clust<-50 #speed up by making this look at smaller
    
    # load any partial results already saved
    if(file.exists(paste0(indir,'/Time_store.rdata'))){load(file=paste0(indir,'/Time_store.rdata'))}
    if(file.exists(paste0(indir,'/BIC_values.rdata'))){load(file=paste0(indir,'/BIC_values.rdata'))}
  
    
    start_kk = max(length(BIC_val)+1,2)
    
    
  
    if(start_kk <= max_clust){ #check if havent gone too far
      print(paste("Any previous results saved, begining from number of clusters =", start_kk,". Search continues up to", max_clust))
      
      for (kk in (start_kk:max_clust)) {
        # Conduct timing while computing BIC values
        TIME <- proc.time()
        BIC_val[kk] <- BIC_lowmem(tkmeans_lowmem(big_mat,kk,.9,1,20),big_mat)
        TIME_DUMMY <- proc.time() - TIME
        print(TIME_STORE)
        TIME_STORE[kk] <- TIME_DUMMY[1] + TIME_DUMMY[2]
        print(c(kk,BIC_val[kk]))
  
        # Save the results
        save(TIME_STORE,file=paste0(indir,'/Time_store.rdata'))
        save(BIC_val,file=paste0(indir,'/BIC_values.rdata'))
      }
    }
    ### Get the optimal K and computer clustering under optimal K
    n <- dim(big_mat)[1]
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
    #clustering <- tkmeans(x=big_mat,k=comp,alpha=.9,nstart=5,iter.max=20)
    n_starts  = 5
    starts_list<-list()
    BIC_list<-array(0, n_starts)
    
    print(paste("Running multiple (", n_starts, ") starts on optimal clustering number,", comp))  
    for(j in 1:n_starts){
      
      starts_list[[j]] <- tkmeans_lowmem(big_mat, comp,.9,1,20) #only one start
      
      BIC_list[j]<-BIC_lowmem(starts_list[[j]] ,big_mat)
      
    }
    
    clustering<-starts_list[[which.min(BIC_list)]]
    #save centres
    save(clustering,file=paste0(indir,'/centers.rdata'))
  }else{
    print("Clustering already done. Loading saved data.")
    load(paste0(indir,'/centers.rdata'))
    comp<-dim(clustering)[1]
    
  }
  ### Function for allocating observations to cluster from a tkmeans clustering
  # tmeansClust <- function(fit,data) {
  #   apply(as.matrix(pdist2(data,t(fit$center))),1,which.min)
  # }
  
  print("Doing some more analysis on clustering results...")
  ### Get a clustering using the tkmeans clustering form variable "clustering"
  clustering_cluster <- tmeansClust_lowmem(big_mat,clustering)
  ## save cluster allocations
  save(clustering_cluster,file=paste(indir,'/clustering.rdata',sep=''))
}

#reload params
load(file = paste0(indir,"/scalingparams.rdata"))
load(file=paste0(indir,'/centers.rdata'))
load(file=paste(indir,'/clustering.rdata',sep=''))
load(file=paste(indir,'/image_hold.rdata',sep=''))
load(file=paste(indir,'/predicted_means.rdata',sep=''))
load(file=paste(indir,'/correlation_matrix.rdata',sep=''))
load(file=paste(indir,'/mean_image.rdata',sep=''))

comp<-dim(clustering)[1]

clustering_cluster <- tmeansClust_lowmem(big_mat,clustering)
save(clustering_cluster, file=paste0(indir,"/cluster_allocations.rdata"))

#print(dim(clustering_cluster))
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

save(image_hold,file=paste(indir,'/image_hold.rdata',sep=''))
image_hold[is.na(image_hold)]<-0
f.write.nifti(image_hold,file=paste0(indir,'/clusters.nii'), nii=TRUE)

image.plot(1:Y_SIZE,1:X_SIZE,image_hold[Z_SIZE,,])
# 
# ### Obtain the cluster mean time series
# # Reload Graphing Parameters (These values are specific to F03 and F04)
 file_number <- file_number.old
 max_file <- max(file_number)
file_number <- (file_number-1)*Z_SIZE + ss
Basis_number <- Basis_number
Basis <- create.bspline.basis(c(0,(max_file-1)*Z_SIZE+Z_SIZE),nbasis=Basis_number)
BS <- eval.basis(file_number,Basis)
# Compute the mean time series
centers <- clustering
#centers <- sweep(centers,2,attributes(big_mat)$'scaled:scale','*')
#centers <- sweep(centers,2,attributes(big_mat)$'scaled:center','+')








#reload scaling if need be
if(!exists("mean_sd_from_unscaled")){
  if(file.exists(paste0(indir,"/scalingparams.rdata"))){
    load(paste0(indir,"/scalingparams.rdata"))
  }
}

centers <- sweep(centers,2,mean_sd_from_unscaled[1,],'*')
centers <- sweep(centers,2,mean_sd_from_unscaled[2,],'+')



pred_range <- eval.basis(seq(1,((max_file-1)*Z_SIZE+Z_SIZE),1000),Basis)

# Each row is a mean time series over the pred_range values
PRED <- matrix(NA,dim(centers)[1],dim(pred_range)[1])
for (ii in 1:dim(centers)[1]) {
  PRED[ii,] <- apply(pred_range,1,function(x) {sum(x*centers[ii,])})
}
# The time average values of each series
MEAN_PRED <- rowMeans(PRED)
save(PRED,file=paste(indir,'/predicted_means.rdata',sep=''))

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

save(CORR,file=paste(indir,'/correlation_matrix.rdata',sep=''))

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

save(mean_image,file=paste(indir,'/mean_image.rdata',sep=''))
mean_image[is.na(mean_image)]<-0
f.write.nifti(mean_image,file=paste0(indir,'/mean_image.nii'), nii=TRUE)


print("Analysis done and saved. Producing and saving plots.")
# Graphing --------------------------------------------

# ### First set of graphs
# # Graph clustering on every 5th slice
# for (slice in seq(5,Z_SIZE,by=5)) {
#   pdf(paste(indir,'/Clustering_Slice_',slice,'.pdf',sep=''),paper='a4r')
#   ### Plot Clustering Image
#   image_hold <- array(NA,c(Z_SIZE,Y_SIZE,X_SIZE))
#   count <- 0
#   for (ss in 1:Z_SIZE) {
#     for (ii in 1:Y_SIZE) {
#       for (jj in 1:X_SIZE) {
#         if (D_Mask[ss,ii,jj]) {
#           count <- count + 1
#           image_hold[ss,ii,jj] <- clustering_cluster[count]
#         }
#       }
#     }
#   }
# 
#   if(!all(is.na(image_hold[slice,,]))){
#     image.plot(image_hold[slice,,],col=tim.colors(comp))
#     dev.off()
#   }
# }
# 
# # Graph Mean Slices on every 5th slice
# for (slice in seq(5,Z_SIZE,by=5)) {
#   pdf(paste(indir,'/Mean_Slice_',slice,'.pdf',sep=''),paper='a4r')
#   ### Plot Mean Image
#   mean_image <- array(NA,c(Z_SIZE,Y_SIZE,X_SIZE))
#   count <- 0
#   for (ss in 1:Z_SIZE)
#   {
#     for (ii in 1:Y_SIZE) {
#       for (jj in 1:X_SIZE) {
#         if (D_Mask[ss,ii,jj]) {
#           count <- count + 1
#           mean_image[ss,ii,jj] <- MEAN_PRED[clustering_cluster[count]]
#         }
#       }
#     }
#   }
#   if(!all(is.na(mean_image[slice,,]))){
#     image.plot(1:Y_SIZE,1:X_SIZE,mean_image[slice,,],col=grey.colors(100,0,1))
#   dev.off()
#   }
# 
# }


#replace here using ggplot2 face command




# 
# if(Z_SIZE >= 100){
#   # Graph the location of clusters on slice 75
#   for (cc in 1:comp) {
#     pdf(paste(indir,'/Location_on_Slice_75_Cluster_',cc,'.pdf',sep=''),paper='a4r')
#     if(!all(is.na(mean_image[75,,]))){
#       image.plot(1:Y_SIZE,1:X_SIZE,mean_image[75,,],col=grey.colors(100,0,1))
#       count <- 0
#       for (ss in 1:Z_SIZE) {
#         for (ii in 1:Y_SIZE) {
#           for (jj in 1:X_SIZE) {
#             if (D_Mask[ss,ii,jj]) {
#               count <- count + 1
#               if (clustering_cluster[count]==cc & ss==75) {
#                 points(ii,jj,pch=15,cex=1,col='green')
#               }
#             }
#           }
#         }
#       }
#       dev.off()
#     }
#   }
#   
#   # Graph the location of clusters on slice 50
#   for (cc in 1:comp) {
#     pdf(paste(indir,'/Location_on_Slice_50_Cluster_',cc,'.pdf',sep=''),paper='a4r')
#     if(!all(is.na(mean_image[50,,]))){
#       image.plot(1:Y_SIZE,1:X_SIZE,mean_image[50,,],col=grey.colors(100,0,1))
#       count <- 0
#       for (ss in 1:Z_SIZE) {
#         for (ii in 1:Y_SIZE) {
#           for (jj in 1:X_SIZE) {
#             if (D_Mask[ss,ii,jj]) {
#               count <- count + 1
#               if (clustering_cluster[count]==cc & ss==50) {
#                 points(ii,jj,pch=15,cex=1,col='green')
#               }
#             }
#           }
#         }
#       }
#       dev.off()
#     }
#   }
#   
#   # Graph the location of clusters on slice 100
#   for (cc in 1:comp) {
#     pdf(paste(indir,'/Location_on_Slice_100_Cluster_',cc,'.pdf',sep=''),paper='a4r')
#     if(!all(is.na(mean_image[100,,]))){
#       image.plot(1:Y_SIZE,1:X_SIZE,mean_image[100,,],col=grey.colors(100,0,1))
#       count <- 0
#       for (ss in 1:Z_SIZE) {
#         for (ii in 1:Y_SIZE) {
#           for (jj in 1:X_SIZE) {
#             if (D_Mask[ss,ii,jj]) {
#               count <- count + 1
#               if (clustering_cluster[count]==cc & ss==100) {
#                 points(ii,jj,pch=15,cex=1,col='green')
#               }
#             }
#           }
#         }
#       }
#       dev.off()
#     }
#   }
# 
# }

#also graph here  with gggplot2


# # Graph the Mean functions
# for (cc in 1:comp) {
#   pdf(paste(indir,'/Cluster_Mean_Function_',cc,'.pdf',sep=''),paper='a4r')
#   plot(seq(1,((max_file-1)*Z_SIZE+Z_SIZE),1000),PRED[cc,],type='l',xlab='Time',ylab='signal',main=cc)
#   dev.off()
# }

##R version issue here
# pred.df<-data.frame(t(PRED))
# names(pred.df) <- paste0("Cluster ", 1:comp)
# pred.df$Time<-seq(0,dim(PRED)[2]-1,1)
# pred.df.flat<-melt(pred.df, id='Time')
# pred.df.flat$Signal<-pred.df.flat$value
# 
# pdf(paste(indir,'/All_Cluster_Mean_Function.pdf',sep=''),paper='a4')
# ggplot(data=pred.df.flat)+geom_line(aes(y=Signal, x=Time))+facet_wrap(~variable)+theme_bw()
# dev.off()

# Plot the Correlation matrix
pdf(paste(indir,'/Correlation_matrix.pdf',sep=''),paper='a4r')
image.plot(1:comp,1:comp,CORR)
dev.off()


# Hierachical Clustering ------------------------------
# # Make a distance metric
DIST <- as.dist(1-t(CORR))
HCLUST <- hclust(DIST,method='average')
pdf(paste(indir,'/Correlation_matrix.pdf',sep=''),paper='a4r')
plot(HCLUST)
dev.off()
print("04-Rstats02.R complete.")
# # 
# # # Make tree cut using Dunn Index
# #NB <- NbClust( diss=DIST,distance=NULL,method='average',index='silhouette',max.nc=ceiling(comp-2))
# NB<-NbClust(data = 1-t(CORR), diss=DIST,distance="NULL",method='average',index='silhouette', max.nc=ceiling(comp-2))
# CUT <- NB$Best.partition
# 
# # Plot dendrogram
# pdf(paste(indir,'/Dendrogram_clusters.pdf',sep=''),width=30,height=10)
# plot(HCLUST,xlab='')
# rect.hclust(HCLUST,max(CUT),border=rainbow(max(CUT)))
# dev.off()

# ### Plot the 10 HCLust Cluster Means
# for (ii in 1:max(CUT)) {
#   pdf(paste(indir,'/HCLUST_',ii,'.pdf',sep=''),paper='a4r')
#   plot(c(1,((max_file-1)*Z_SIZE+Z_SIZE)),c(min(PRED),max(PRED)),
#        type='n',xlab='Time',ylab='signal',main=ii)
#   for (ss in 1:comp) {
#     if (CUT[ss]==ii) {
#       lines(seq(1,((max_file-1)*Z_SIZE+Z_SIZE),1000),PRED[ss,],col=tim.colors(comp)[ss])
#     }
#   }
#   dev.off()
# }
# 
# # Plot Frequency Histogram
# pdf(paste(indir,'/Frequency_of_clusters.pdf',sep=''),paper='a4r')
# plot(table(clustering_cluster),xlab='cluster',ylab='Frequency')
# dev.off()

# # Plot K means by HCLUST reference
# pdf(paste(indir,'/Cluster_by_HCLUST.pdf',sep=''),paper='a4r')
# plot(CUT,col=tim.colors(comp),xlim=c(0,comp+1),ylim=c(0,8),xlab='Cluster',ylab='HCLUST')
# grid()
# dev.off()


# # Graph clustering based on cuts on every 5th slice
# for (slice in seq(5,Z_SIZE,by=5)) {
#   pdf(paste(indir,'/CUT_Slice_',slice,'.pdf',sep=''),paper='a4r')
#   ### Plot Clustering Image
#   image_hold <- array(NA,c(Z_SIZE,Y_SIZE,X_SIZE))
#   count <- 0
#   for (ss in 1:Z_SIZE) {
#     for (ii in 1:Y_SIZE) {
#       for (jj in 1:X_SIZE) {
#         if (D_Mask[ss,ii,jj]) {
#           count <- count + 1
#           image_hold[ss,ii,jj] <- CUT[clustering_cluster[count]]
#         }
#       }
#     }
#   }
# 
#   if(!all(is.na(image_hold[slice,,]))){
#     image.plot(image_hold[slice,,],col=tim.colors(max(CUT)))
#     dev.off()
#   }
# }
