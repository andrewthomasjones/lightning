list.of.packages <- c('parallel', 'lattice', 'AnalyzeFMRI', 'ggplot2', 'reshape2', 'MASS', 'abind', 'fda', 'fields', 'speedglm','pracma', 'tclust', 'signal', 'capushe', 'pryr', 'lowmemtkmeans', 'matrixStats')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

### Load necessary libraries (all CRAN libraries can be acquired using install.packages)
library(lattice)
library(AnalyzeFMRI)
library(ggplot2)
library(reshape2)
library(MASS)
library(abind)
library(fda)
library(fields)
library(speedglm)
library(pracma)
library(tclust)
library(signal)
library(capushe)
library(pryr)
library(lowmemtkmeans)
library(parallel)
library(matrixStats)



### Parse CMD Line Arguments
args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
outdir <- args[2]
FROM <- as.numeric(args[3])
TO <- as.numeric(args[4])
mask_fn <- args[5]

min_activity<-7
activity_range<-0.022
Basis_number <- 100
#Basis_number <- floor(100*(TO-FROM)/2000)
mask_cutoff <-.8
ylimits1<-(c(-0.03,0.03))
f1<-20

X_START <- 0
Y_START <- 0
Z_START <- 0

X_SIZE <- 640#640
Y_SIZE <- 130 #130
Z_SIZE <- 300#300

print(c('1indir', indir))
print(c('2outdir', outdir))
print(c('3FROM', FROM))
print(c('4TO', TO))
print(c('5mask', mask_fn))

header_fn<-paste0(indir, '/R-out-000000.nii')
header<-f.read.header(header_fn)

print(outdir)
if(!dir.exists(outdir)){
  dir.create(outdir, recursive=T)
}

### Declare what files to load in
### Remove time points that do not look stationary
file_list_initial <- list.files(path=indir, pattern='R-out-.*.nii', full.names=FALSE)
file_number_initial <- as.numeric(substr(file_list_initial,7,12))
file_list <- file_list_initial[which(file_number_initial>FROM & file_number_initial<=TO)]
file_number <- as.numeric(substr(file_list,7,12)) - FROM


max_file <- max(file_number)
file_number.old <- file_number
file_number <- (file_number-1)*Z_SIZE + Z_SIZE+Z_START

### Declare number of time slices left and initialize list
N <- length(file_number)
print(paste(N, "timeslices to read."))

if(!file.exists(file=paste(outdir,'/settings.rdata',sep=''))){
  save(header,N,max_file, file_number.old, file_list , X_START,Y_START,Z_START,X_SIZE,Y_SIZE,Z_SIZE,indir,outdir,FROM ,header,f1,ylimits1,TO , mask_cutoff,mask_fn ,min_activity,activity_range,Basis_number,file=paste(outdir,'/settings.rdata',sep=''))
}else{
  load(file=paste(outdir,'/settings.rdata',sep=''))
}

dir.create(paste0(outdir,'/params/'))
dir.create(paste0(outdir,'/mips/'))
dir.create(paste0(outdir,'/coeff_mats/'))


if(!file.exists(file=paste(outdir,'/D_Mask.rdata',sep=''))){
  ### Load Mask File
  print(paste("Loading mask ", mask_fn, sep=""))
  MASK <- f.read.nifti.volume(mask_fn)
  MASK_hdr <- f.read.nifti.header(mask_fn)
  MASK <- MASK[,,,1]
  ### Dummy Mask
  D_Mask <- array(NA,c(X_SIZE,Y_SIZE, Z_SIZE))
  
  for (s in Z_START:(Z_SIZE+Z_START)) {
    for (j in Y_START:(Y_SIZE+Y_START)) {
      for (i in X_START:(X_SIZE+X_START))  {
        D_Mask[i-X_START, j-Y_START, s-Z_START] <- MASK[i, j, s]>= mask_cutoff#.99995
      }
    }
  }
  ssDM <- sum(D_Mask)
  
  save(ssDM, file=paste(outdir,'/ssDM.rdata',sep=''))
  save(MASK_hdr, file=paste(outdir,'/MASK_hdr.rdata',sep=''))
  save(D_Mask, file=paste(outdir,'/D_Mask.rdata',sep=''))
  f.write.nifti(D_Mask,file=paste0(outdir,'/binary_mask.nii'), nii=TRUE, L=header)
  
}else{
  load(file=paste(outdir,'/ssDM.rdata',sep=''))
  load(file=paste(outdir,'/MASK_hdr.rdata',sep=''))
  load(file=paste(outdir,'/D_Mask.rdata',sep=''))
}

# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
#cl <- makeCluster(no_cores)
print(no_cores)

### Get already existing outputs
#Completed_output <- list.files(path=paste0(outdir, '/coeff_mats/'), pattern='coeff_mat', full.names=FALSE)
#Completed_numbers<-(as.numeric(unlist(strsplit(unlist(Completed_output), "[^0-9]+")))
                    #print(Completed_numbers)
                    
                    
                    if(!file.exists(paste(outdir,'/MASK_active.rdata',sep=''))){
                      print("Creating mask files")
                      active_mask1 <- array(NA,c(X_SIZE,Y_SIZE, Z_SIZE))
                      active_mask2 <- array(NA,c(X_SIZE,Y_SIZE, Z_SIZE))
                      active_mask3 <- array(NA,c(X_SIZE,Y_SIZE, Z_SIZE))
                      active_mask4 <- array(NA,c(X_SIZE,Y_SIZE, Z_SIZE))
                      k<-1
                      save(active_mask1,active_mask2,active_mask3,active_mask4,k, file=paste(outdir,'/MASK_active.rdata',sep=''))
                    }else{
                      load(file=paste(outdir,'/MASK_active.rdata',sep=''))
                    }
                    
                    
                    
                    ### Fit Splines to the time series in each of the slice files (s loops over slices)
                    for (s in k:Z_SIZE) {
                      
                      print(mem_used())
                      rm(list = ls()[-c(which(ls()=="outdir"), which(ls()=="s"))])
                      gc()
                      load(file=paste(outdir,'/MASK_hdr.rdata',sep=''))
                      load(file=paste(outdir,'/settings.rdata',sep=''))
                      load(file=paste(outdir,'/ssDM.rdata',sep=''))
                      load(file=paste(outdir,'/D_Mask.rdata',sep=''))
                      load(file=paste(outdir,'/MASK_active.rdata',sep=''))
                      print(mem_used())
                      print(c('Doing ', s))
                      
                      ### Get already existing outputs
                      Completed_output <- list.files(path=paste0(outdir, '/coeff_mats/'), pattern='coeff_mat', full.names=FALSE)
                      #Completed_output <- Completed_output[grep('coeff_mat',Completed_output)]
                      ### Skip things already in Completed_output
 
                        file_number <- (file_number.old-1)*Z_SIZE +s
                        
                        full_array <- list()
                        
                        ### Read in Z slices
                        for (j in 1:N) {
                          file_name <- substring(file_list[j],1,12)
                          nii_name <- paste(indir, '/', file_name,'.nii',sep='')
                          #print(c('reading',nii_name, s+Z_START,j,round(j/N),1))
                          if(j %% 25 == 0){
                            print(paste("Reading Slice", s+Z_START, round(100*(j/N),1), '% Complete'))
                          }
                          NIFTI <- f.read.nifti.slice(nii_name,s+Z_START,1)
                          full_array[[j]] <- NIFTI
                        }
                        
                        
                        ### Combine all of the time slices
                        full_array <- abind(full_array,along=3)
                        
                        print('storing')
                        ### Store data into matrix instead of array
                        count <- 0
                        full_mat <- matrix(NA,prod(Y_SIZE*X_SIZE),N)
                        for (j in 1:Y_SIZE) {
                          for (i in 1:X_SIZE) {
                            count <- count + 1
                            full_mat[count,] <- full_array[i+X_START,j+Y_START,]
                            #print(c('storing',s,count))
                          }
                        }
                        rm(full_array)
                        
                        odds<-seq(1,ncol(full_mat),2)
                        ## take every second point
                        mat_odds<-full_mat[,odds]
                        mat_odds_detr<-array(NA,dim(mat_odds))
                        mat_odds_filt<-array(NA,dim(mat_odds))
                        
                        snap_mips <- array(NA,c(X_SIZE,Y_SIZE))
                        snap_mips_detrend <- array(NA,c(X_SIZE,Y_SIZE))
                        snap_mips_filt <- array(NA,c(X_SIZE,Y_SIZE))
                        
                        ### Detrend data
                        print('detrending')
                        count<-0
                        
                        for (j in 1:Y_SIZE) {
                          for (i in 1:X_SIZE) {
                            count <- count + 1
                            mat_odds_detr[count,] <- detrend(mat_odds[count,], 'linear')
                          }
                        }
                        
                        bf <- butter(5, (1/f1), type="high")
                        
                        print('filtering')
                        count<-0
                        
                        for (j in 1:Y_SIZE) {
                          for (i in 1:X_SIZE) {
                            count <- count + 1
                            mat_odds_filt[count,] <- filter(bf, mat_odds_detr[count,])
                          }
                        }
                        
                        count<-0
                        
                        #temp1<- quantile(mat_odds[count,], probs = seq(0.025, 0.975),na.rm=T)
                        temp1<- rowQuantiles(mat_odds,  probs = c(0.025, 0.975),na.rm=T)
                        print(paste("Looking at activity levels in layer ", s))
                        for (j in 1:Y_SIZE) {
                          for (i in 1:X_SIZE) {
                            count <- count + 1
                            active_mask1[i,j,s] <-  (max(mat_odds_filt[count,],na.rm=T) - min(mat_odds_filt[count,],na.rm=T))
                            active_mask2[i,j,s] <-  (max(mat_odds[count,],na.rm=T) - min(mat_odds[count,],na.rm=T))
                            active_mask3[i,j,s] <-  temp1[i,2]-temp1[i,1]
                            active_mask4[i,j,s] <-  var(mat_odds_filt[count,],na.rm=T)
                            #print(c(active_mask1[i,j,s],active_mask2[i,j,s],active_mask3[i,j,s],active_mask4[i,j,s]))
                          }
                        }
                        print("saving")
                        print(mem_used())
                        k<-s+1
                        save(active_mask1,active_mask2,active_mask3,active_mask4, k, file=paste(outdir,'/MASK_active.rdata',sep=''))
                        
                        
                      }
                      

                      