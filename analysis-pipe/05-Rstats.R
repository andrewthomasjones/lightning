### Load necessary libraries (all CRAN libraries can be acquired using install.packages)
library(compiler)
library(AnalyzeFMRI)

library(MASS)
library(abind)
library(fda)
library(fields)
library(speedglm)
library(pracma)
library(lattice)
library(tclust)
library(NbClust)
library(capushe)

### Parse CMD Line Arguments
args <- commandArgs(trailingOnly = TRUE)

indir <- args[1]
mask_fn <- args[2]

print(c('1indir', indir))

X_SIZE <- 130
Y_SIZE <- 640
Z_SIZE <- 300

load(file = paste0(indir,"Scalingparams.rdata"))
load(file=paste0(indir,'/centers.rdata'))
load(file=paste(indir,'/clustering.rdata',sep=''))
load(file=paste(indir,'/image_hold.rdata',sep=''))
load(file=paste0(indir,"Scalingparams.rdata"))
load(file=paste(indir,'/predicted_means.rdata',sep=''))
load(file=paste(indir,'/correlation_matrix.rdata',sep=''))
load(file=paste(indir,'/mean_image.rdata',sep=''))


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


comp<-dim(clustering)[1]


clustering_cluster <- tmeansClust_lowmem(big_mat,clustering)    


# Define a function that obtains the Moore neighborhood
# of distance 'dist' for any coordinate (x_co,y_co)
Moore_nbh <- function(input,dist,x_co,y_co,z_co) {
  C_xyz <- input[x_co, y_co, z_co]
  
  x_size = dim(input)[1]
  y_size = dim(input)[2]
  z_size = dim(input)[3]
  
  #print(paste0("(",x_size,",", y_size,",", z_size,")"))
  
  output <- input[max(1,x_co-dist):min(x_size,x_co+dist),max(1,y_co-dist):min(y_size,y_co+dist), max(1,z_co-dist):min(z_size,z_co+dist)]
  output <- c(output)
  return(output)
}

# Define a function that obtains the von Neumann
# neighborhood of distance 'dist' for any coordinate
# (x_co,y_co).
vN_nbh <- function(input,dist,x_co,y_co, z_co) {
  C_xyz <- input[x_co, y_co, z_co]
  
  x_size = dim(input)[1]
  y_size = dim(input)[2]
  z_size = dim(input)[3]
  output <- c()
  
  count <- 1
  list_x <- max(1,x_co-dist):min(x_size ,x_co+dist)
  list_y <- max(1,y_co-dist):min(y_size ,y_co+dist)
  list_z <- max(1,z_co-dist):min(z_size ,z_co+dist)
  
  for (x in list_x) {
    for (y in list_y) {
      for (z in list_z) {
        centre_dist = abs(x-x_co) + abs(y-y_co) + abs(z-z_co)
        if (centre_dist <= dist) {
          output[count] <- input[x,y,z]
          count <- count + 1
        }
      }
    }
  }
  
  return(output)
}

# Create a variable for the function evaluates eta_fun and
# put into the function, the average count from each of the
# 4 classes, within the Moore neighborhood of distance 1
# around each coordinate (x,y).
eta_fun <- function(input,dist, weighting) {
  count <- 1
  classes <- max(input,na.rm=TRUE)
  
  n<-length(input)
  
  x_size = dim(input)[1]
  y_size = dim(input)[2]
  z_size = dim(input)[3]
  
  output <- matrix(NA,n, classes)
  
  for (x in 1:x_size) {
    for (y in 1:y_size) {
      for(z in 1:z_size) {
        
        for(i in 1:classes){
            output[count,i] <- weighting[i]*mean(Moore_nbh(input,dist,x,y, z)==i,na.rm=T) 
        }
        
        count <- count + 1
      }
    }
  }
  return(output)
}  
  
mrf_val <- eta_fun(image_hold,1)
# Construct an MRF from the image and neighbourhood function
# eta_fun using the multinomial regression function multinom.

df <- data.frame(y=c(image_hold),x=mrf_val)
MRF_1 <- multinom(y~., data=df) 

# Obtain the parameter estimates from the fitted MRF model.
summary(MRF_1)

# Compute the PLIC value of the fitted MRF model.
BIC(MRF_1)

# Plot an image of the obtained estimated signals from the
# MRF model.
smooth_im_1 <- (predict(MRF_1 ))
#table(smooth_im_1 )
#table(image_hold)


testout1<-image_hold
testout1[is.na(testout1)]<-0
testout2<-array(smooth_im_1, dim = dim(testout1))

#also stuff is getting jumbled here

f.write.nifti(testout1 ,file='./testorig.nii', nii=TRUE)
f.write.nifti(testout2 ,file='./testorig2.nii', nii=TRUE)


