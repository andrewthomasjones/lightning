### Load necessary libraries (all CRAN libraries can be acquired using install.packages)
library(compiler)
library(AnalyzeFMRI)
#library(ggplot2)
#library(reshape2)
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
#if(file.exists('./temp/tkmeans.cpp')){sourceCpp('./temp/tkmeans.cpp')} #for testing
#if(file.exists('./tkmeans.cpp')){sourceCpp('./tkmeans.cpp')} #actual file structure
#sourceCpp('/data/nif02/uqajon14/temp/tkmeans.cpp')
### Enable Just In Time Compiling, for ????MAYBE???? more speed.
enableJIT(3)

### Parse CMD Line Arguments
args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
outdir <- args[2]

### These lines are kind of data dependent. I've determined what works for F11, F03, and F04.
## FROM and TO determine the start and end points of the trimming
FROM <- as.numeric(args[3])
TO <- as.numeric(args[4])

mask_fn <- args[5]

print(c('1indir', indir))
print(c('2outdir', outdir))
print(c('3FROM', FROM))
print(c('4TO', TO))
print(c('5mask', mask_fn))

# #
# X_SIZE <- 70
# Y_SIZE <- 170
# Z_SIZE <- 150

# new crop sizes BAD post-doc... Swap X and Y later...
#this will be approximately 18.6 GB in RAM with 100 basis function params

X_SIZE <- 130
Y_SIZE <- 640
Z_SIZE <- 300

# X_SIZE <- 100
# Y_SIZE <- 250
# Z_SIZE <- 175


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
#if(!file.exists(paste0(outdir,'/clustering.rdata')) | !file.exists(paste0(outdir,'/centers.rdata')) ){
  
  #print("Pre-saved clustering doesn't exist. Beginning clustering process..")
  ### Set working directory for computation
  # BAD post-doc, do not do this. (this made things fail later)
  # see changes to list.files below also
  # setwd(indir)
  
  # SPLINE FITTING --------------------------------------
  
z_set = c(100,150,200,250)
  ### Fit Splines to the time series in each of the slice files (ss loops over slices)
  #for (ss in 1:Z_SIZE) {
for (ss in z_set) {
    print(c('Doing ', ss))
 
      ### Declare what files to load in
      ### Remove time points that do not look stationary
      file_list <- list.files(path=indir, pattern='R-out-.*.nii', full.names=FALSE)
      #file_list <- file_list[grep('.nii',file_list)]
      file_number <- as.numeric(substr(file_list,7,12))
      
      if(FROM==0){# what do when from and to are both zero basically
        if(TO==0){
          TO_2 = min(3500,max(file_number))
        }else{TO_2 = TO}
        
        file_list <- file_list[which(file_number>=1000& file_number<=TO_2)] #this throws out first image no matter what, maybe change
      }else{
          file_list <- file_list[which(file_number>FROM & file_number<=TO)]
      }
      
      file_number <- as.numeric(substr(file_list,7,12)) - FROM
  
      file_number.old <- file_number
      max_file <- max(file_number)
      file_number <- (file_number-1)*Z_SIZE + ss
  
      ### Declare number of time slices left and initialize list
      N <- length(file_number)
      print(paste(N, "timeslices to read."))
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
      write.csv(full_mat,file=paste(outdir,'/time_series',ss,'.csv',sep=''))
  
      ### Detrend data
      #print('detrending')
      #for (ii in 1:(Y_SIZE*X_SIZE)) {
       # full_mat[ii,] <- full_mat[ii,]-speedlm.fit(full_mat[ii,],cbind(1,file_number))$coefficients[2]*file_number
        # full_mat[ii,] <- detrend(full_mat[ii,])
        # print(c('detrending',ss,ii))
      #}
  
  
      ### Declare number of bases to use
      #Basis_number <- 100
      #Basis <- create.bspline.basis(c(0,(max_file-1)*Z_SIZE+Z_SIZE),
                                    #nbasis=Basis_number)
      #BS <- eval.basis(file_number,Basis)
  
  
      ### Fit B-spline to all time series
      #print(c('Spline fitting',ss))
      #FD <- smooth.basisPar(file_number,t(full_mat),Basis)
      #coeff_mat <- as.matrix(t(FD$fd$coefs))
  
  
      ### Save the results in the format coeff_mat_(SLICE NUMBER).rdata
      #print(paste('Saving as: ', outdir, '/coeff_mat_', ss, '.rdata', sep=''))
      
    }
  
 # }
  

