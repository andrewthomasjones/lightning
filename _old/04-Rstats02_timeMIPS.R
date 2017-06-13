
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
#library(Rcpp)
#library(RcppEigen)
#library(RcppArmadillo)

library(lowmemtkmeans)
#needs to be c++ 11 because I used the tuple class
#Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
#if(file.exists('./temp/tkmeans.cpp')){sourceCpp('./temp/tkmeans.cpp')} #for testing
#if(file.exists('./tkmeans.cpp')){sourceCpp('./tkmeans.cpp')} #actual file structure
#sourceCpp('./tkmeans.cpp')
#sourceCpp('/data/nif02/uqajon14/temp/tkmeans.cpp')
### Enable Just In Time Compiling, for ????MAYBE???? more speed.
#enableJIT(3)

### Parse CMD Line Arguments
# args <- commandArgs(trailingOnly = TRUE)
# indir <- args[1]
# outdir <- args[2]
# FROM <- as.numeric(args[3])
# TO <- as.numeric(args[4])
# mask_fn <- args[5]

#NAME<-'/media/andrew/Port/20150628_5dpf_H2BS_CON_LR_F03/'
NAME<-'/media/andrew/Port/20150629_6dpf_H2BS_PTZ_LR_F1_1/'

indir<-  paste0(NAME,"03-chunk02_old")
outdir <- paste0(NAME,"TIMEMIPS")
FROM <- 55
TO <- 2200
mask_fn  <-paste0(indir, "/mask.nii")

if(!dir.exists(outdir)){
  dir.create(outdir)
}

print(c('1indir', indir))
print(c('2outdir', outdir))
print(c('3FROM', FROM))
print(c('4TO', TO))
print(c('5mask', mask_fn))

# #
# X_SIZE <- 10
# Y_SIZE <- 10
# Z_SIZE <- 10
# X_SIZE <- 40
# Y_SIZE <- 100
# Z_SIZE <- 50

# X_START <- 50
# Y_START <- 300
# Z_START <- 150

X_START <- 0
Y_START <- 0
Z_START <- 0

X_SIZE <- 130
Y_SIZE <- 640
Z_SIZE <- 300

# new crop sizes BAD post-doc... Swap X and Y later...
#this will be approximately 18.6 GB in RAM with 100 basis function params

# X_SIZE <- 100
# Y_SIZE <- 250
# Z_SIZE <- 175


# ### Load Mask File
# print(paste("Loading mask ", mask_fn, sep=""))
# MASK <- f.read.nifti.volume(mask_fn)
# MASK <- MASK[,,,1]
# ### Dummy Mask
# D_Mask <- array(NA,c(Z_SIZE,Y_SIZE,X_SIZE))
# 
# for (ss in Z_START:(Z_SIZE+Z_START)) {
#   for (ii in Y_START:(Y_SIZE+Y_START)) {
#     for (jj in X_START:(X_SIZE+X_START))  {
#       D_Mask[ss-Z_START,ii-Y_START,jj-X_START] <- MASK[ii, jj, ss]>=.99995
#     }
#   }
# }
# ssDM <- sum(D_Mask)


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

snap_mips <- array(NA,c(Z_SIZE,Y_SIZE,X_SIZE)) 
snap_mips_detrend <- array(NA,c(Z_SIZE,Y_SIZE,X_SIZE))
progress<-0 

if(file.exists(paste(outdir,'/progress.rdata',sep=''))){
  load(file=paste(outdir,'/progress.rdata',sep=''))
  load(file=paste(outdir,'/snap_mips.rdata',sep=''))
  load(file=paste(outdir,'/snap_mips_detrend.rdata',sep=''))
  }

start_val=progress+1
  ### Fit Splines to the time series in each of the slice files (ss loops over slices)
  for (ss in start_val:Z_SIZE) {
    print(c('Doing ', ss))
    
    ### Get already existing outputs
    Completed_output <- list.files(path=outdir, pattern='coeff_mat', full.names=FALSE)
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
        print(c('reading',nii_name, ss+Z_START,ii,round(ii/N),1))
        NIFTI <- f.read.nifti.slice(nii_name,ss+Z_START,1)
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
          full_mat[count,] <- full_array[ii+Y_START,jj+X_START,]
          print(c('storing',ss,count))
        }
      }
        ### Detrend data
      print('detrending')
      count<-0
      full_mat2<-full_mat
      
	    for (ii in 1:Y_SIZE) {
        for (jj in 1:X_SIZE) {
          count <- count + 1

          full_mat2[count,] <- detrend(full_mat[count,], 'linear')
       
	        snap_mips[ss,ii, jj] <- max(full_mat[count,],na.rm=T)
	        snap_mips_detrend[ss,ii, jj]<-max(full_mat2[count,],na.rm=T)
      	}
	    }
      progress<-ss
      save(progress, file=paste(outdir,'/progress.rdata',sep=''))
      save(snap_mips, file=paste(outdir,'/snap_mips.rdata',sep=''))
      save(snap_mips_detrend, file=paste(outdir,'/snap_mips_detrend.rdata',sep=''))
    }

  }

snap_mips_detrend[is.na(snap_mips_detrend)]<-0
snap_mips[is.na(snap_mips)]<-0
  
f.write.nifti(snap_mips_detrend,file=paste0(outdir,'/time_mips_detrend.nii'), nii=TRUE)
f.write.nifti(snap_mips, file=paste0(outdir,'/time_mips.nii'), nii=TRUE)

print("Done")
