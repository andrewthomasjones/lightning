
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
library(devtools)
install_github("andrewthomasjones/lowmemtkmeans")
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

#NAME<-'/media/andrew/Port/20150629_6dpf_H2BS_PTZ_LR_F1_1/'
NAME<-'/media/andrew/Port/20150628_5dpf_H2BS_CON_LR_F03/'
indir<-  paste0(NAME,"03-chunk02_old")
outdir <- paste0(NAME,"slice_147c")
FROM <- 60
TO <- 1200 #2200 #1200
mask_fn  <-paste0(indir, "/mask.nii")
min_activity<-0
Basis_number <- 100
mask_cutoff <-.8
ylimits1<-(c(-0.03,0.03))

X_START <- 0
Y_START <- 0
Z_START <- 146


X_SIZE <- 640#640
Y_SIZE <- 130 #130
Z_SIZE <- 1 #300


### These lines are kind of data dependent. I've determined what works for F11, F03, and F04.
## FROM and TO determine the start and end points of the trimming


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


if(!dir.exists(outdir)){
  dir.create(outdir)
}
if(!file.exists(file=paste(outdir,'/settings.rdata',sep=''))){
  save(X_START,Y_START,Z_START,X_SIZE,Y_SIZE,Z_SIZE,NAME,indir,outdir,FROM ,ylimits1,TO ,mask_cutoff,mask_fn ,min_activity,Basis_number,file=paste(outdir,'/settings.rdata',sep=''))
}else{
  load(file=paste(outdir,'/settings.rdata',sep=''))
}


dir.create(paste0(outdir,'/cluster_timeseries_plots/'))
dir.create(paste0(outdir,'/images/'))
dir.create(paste0(outdir,'/clusters/'))
dir.create(paste0(outdir,'/cluster_functions/'))
dir.create(paste0(outdir,'/params/'))
dir.create(paste0(outdir,'/mips/'))
dir.create(paste0(outdir,'/coeff_mats/'))





#save master params list


# new crop sizes BAD post-doc... Swap X and Y later...
#this will be approximately 18.6 GB in RAM with 100 basis function params



# X_SIZE <- 100
# Y_SIZE <- 250
# Z_SIZE <- 175


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
f.write.nifti(D_Mask,file=paste0(outdir,'/binary_mask.nii'), nii=TRUE)


MASK_hdr2<-MASK_hdr
MASK_hdr2$dim<-c(3,X_SIZE,Y_SIZE,Z_SIZE,1,1,0,0)
save(MASK_hdr2, file=paste(outdir,'/MASK_hdr2.rdata',sep=''))

### Declare what files to load in
### Remove time points that do not look stationary
file_list <- list.files(path=indir, pattern='R-out-.*.nii', full.names=FALSE)
#file_list <- file_list[grep('.nii',file_list)]
file_number <- as.numeric(substr(file_list,7,12))

file_list <- file_list[which(file_number>FROM & file_number<=TO)]
file_number <- as.numeric(substr(file_list,7,12)) - FROM

file_number.old <- file_number
max_file <- max(file_number)
file_number <- (file_number-1)*Z_SIZE + s


#probably put a better check here
if(!file.exists(paste0(outdir,'/params/clustering.rdata')) | !file.exists(paste0(outdir,'/params/centers.rdata')) ){
  
  print("Pre-saved clustering doesn't exist. Beginning clustering process..")
  #clear and reload
  
  
  ### Set working directory for computation
  # BAD post-doc, do not do this. (this made things fail later)
  # see changes to list.files below also
  # setwd(indir)
  
  # SPLINE FITTING --------------------------------------
  #snap_500 <- array(NA,c(X_SIZE,Y_SIZE, Z_SIZE))
  
  active_mask <- array(NA,c(X_SIZE,Y_SIZE, Z_SIZE))
  ### Fit Splines to the time series in each of the slice files (s loops over slices)
  for (s in 1:Z_SIZE) {
    rm(list = ls()[-c(which(ls()=="outdir"), which(ls()=="active_mask"), which(ls()=="s"))])
    gc()
    load(file=paste(outdir,'/MASK_hdr.rdata',sep=''))
    load(file=paste(outdir,'/MASK_hdr2.rdata',sep=''))
    load(file=paste(outdir,'/settings.rdata',sep=''))
    load(file=paste(outdir,'/ssDM.rdata',sep=''))
    
    
    print(c('Doing ', s))
    ### Get already existing outputs
    Completed_output <- list.files(path=paste0(outdir, '/coeff_mats/'), pattern='coeff_mat', full.names=FALSE)
    #Completed_output <- Completed_output[grep('coeff_mat',Completed_output)]
    ### Skip things already in Completed_output
    
    if (paste('coeff_mat_',s,'.rdata',sep='')%in%Completed_output){
      print(paste('coeff_mat_', s, '.rdata', ' exists, loading', sep=''))
      
      
      ### Declare what files to load in
      ### Remove time points that do not look stationary
      file_list <- list.files(path=indir, pattern='R-out-.*.nii', full.names=FALSE)
      #file_list <- file_list[grep('.nii',file_list)]
      file_number <- as.numeric(substr(file_list,7,12))
      
      file_list <- file_list[which(file_number>FROM & file_number<=TO)]
      file_number <- as.numeric(substr(file_list,7,12)) - FROM
      
      file_number.old <- file_number
      max_file <- max(file_number)
      file_number <- (file_number-1)*Z_SIZE + s
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
      file_number <- (file_number-1)*Z_SIZE + s
      
      ### Declare number of time slices left and initialize list
      N <- length(file_number)
      print(paste(N, "timeslices to read."))
      full_array <- list()
      
      ### Read in Z slices
      for (j in 1:N) {
        file_name <- substring(file_list[j],1,12)
        nii_name <- paste(indir, '/', file_name,'.nii',sep='')
        print(c('reading',nii_name, s+Z_START,j,round(j/N),1))
        NIFTI <- f.read.nifti.slice(nii_name,s+Z_START,1)
        full_array[[j]] <- NIFTI
      }
      
      
      ### Combine all of the time slices
      full_array <- abind(full_array,along=3)
      
      
      ### Store data into matrix instead of array
      count <- 0
      full_mat <- matrix(NA,prod(Y_SIZE*X_SIZE),N)
      for (j in 1:Y_SIZE) {
        for (i in 1:X_SIZE) {
          count <- count + 1
          full_mat[count,] <- full_array[i+X_START,j+Y_START,]
          print(c('storing',s,count))
        }
      }
      
      odds<-seq(1,ncol(full_mat),2)
      ## take every second point
      mat_odds<-full_mat[,odds]
      mat_odds_detr<-array(NA,dim(mat_odds))
      #mat_evens<-full_mat[,-odds]
      #mat_evens_detr<-array(NA,dim(mat_evens))
      #full_mat_detr<-array(NA,dim(full_mat))
      
      snap_mips <- array(NA,c(X_SIZE,Y_SIZE))
      snap_mips_detrend <- array(NA,c(X_SIZE,Y_SIZE))
      
      ### Detrend data
      print('detrending')
      count<-0
      
      for (j in 1:Y_SIZE) {
        for (i in 1:X_SIZE) {
          count <- count + 1
          
          mat_odds_detr[count,] <- detrend(mat_odds[count,], 'linear')
          #mat_evens_detr[count,] <- detrend(mat_evens[count,], 'linear')
        }
      }
      
      
      full_mat_detr<-mat_odds_detr
      #full_mat_detr[,odds]<-mat_odds_detr #SPLIT
      #full_mat_detr[,-odds]<-mat_evens_detr
      
      count<-0
      
      for (j in 1:Y_SIZE) {
        for (i in 1:X_SIZE) {
          count <- count + 1
          active_mask[i,j,s] <-  (max(full_mat_detr[count,],na.rm=T) - min(full_mat_detr[count,],na.rm=T)) > min_activity*sd(full_mat_detr[count,], na.rm=T)
          snap_mips_detrend[i,j] <- max(full_mat_detr[count,],na.rm=T)
          snap_mips[i,j] <- max(mat_odds[count,],na.rm=T)
        }
      }
      
      #full_mat2[count,] <- full_mat2[count,]-speedlm.fit(full_mat2[count,],cbind(1,file_number))$coefficients[2]*file_number
      # full_mat[j,] <- detrend(full_mat[j,])
      # print(c('detrending',s,j))

      
      #snap_500[i,j,s] <- full_mat2[count,250]
    
    
   
    
    
    
    
    pred_splines <- matrix(NA,X_SIZE*Y_SIZE,length(file_number))
    
    
    
    
    ### Declare number of bases to use
    
    Basis <- create.bspline.basis(c(0,(max_file-1)*Z_SIZE+Z_SIZE),nbasis=Basis_number)
    BS <- eval.basis(file_number,Basis)
    # Each row is a mean time series over the pred_range values
    
    
    
    
    ### Fit B-spline to all time series
    print(c('Spline fitting',s))
    
    which_odds<-(1:length(file_number))%%2 #SPLIT
    #which_odds<-(1:length(file_number))
    
    #print(dim(t(full_mat2)))
    
    FD <- smooth.basisPar(file_number[!!which_odds],t(mat_odds_detr),Basis)
    #FD <- smooth.basisPar(file_number,t(full_mat_detr),Basis)
    
    coeff_mat <- as.matrix(t(FD$fd$coefs))
    
    for (i in 1:(X_SIZE*Y_SIZE)){
      pred_splines[i,] <- apply(BS ,1,function(x) {sum(x*coeff_mat[i,])})
    }
    
    pred_splines2<-data.frame(t(pred_splines))
    pred_splines2$time<-seq(from=0,to=0.2*(TO-FROM), length.out=(length(file_number)))
    
    dir.create(paste0(outdir,'/individual_plots/'))
    #ggplot()+geom_line(data=pred_splines2,aes(y=X5,x=time), color='blue')+ geom_point(aes(x=2*(1:dim(full_mat2)[2]-1), y=full_mat2[5,]))
    random_splines<-sample(X_SIZE*Y_SIZE,20)
    for(i in 1:20){
      p<-paste0("X",random_splines[i])
      plot<-ggplot()+geom_point(aes(x=seq(from=0,to=0.2*(TO-FROM), length.out=(length(file_number)/2)), y=mat_odds_detr[random_splines[i],]))+geom_line(data=pred_splines2,aes(y=get(p),x=time), color='blue') +scale_y_continuous("Intensity", limits=ylimits1)+scale_x_continuous("Time (s)")
      pdf(paste0(outdir,'/individual_plots/point_', random_splines[i], '_slice_', s+Z_START, '.pdf',sep=''),paper='a4r')
      print(plot)
      dev.off()
    }
    
    
    
    ### Save the results in the format coeff_mat_(SLICE NUMBER).rdata
    print(paste('Saving as: ', outdir, 'coeff_mats/coeff_mat_', s, '.rdata', sep=''))
    save(coeff_mat,file=paste(outdir,'/coeff_mats/coeff_mat_',s,'.rdata',sep=''))
    save(snap_mips, file=paste(outdir,'/mips/snap_mips_',s,'.rdata',sep=''))
    save(snap_mips_detrend, file=paste(outdir,'/mips/snap_mips_detrend_',s,'.rdata',sep=''))
    #save(snap_500, file=paste(outdir,'/images/snap_500.rdata',sep=''))
    #save(snap_mips, file=paste(outdir,'/images/snap_mips.rdata',sep=''))
    
  }
    
    
}
  if(!file.exists(file=paste0(outdir,'/mips/time_mips.nii'))){
    mips <- array(NA,c(X_SIZE,Y_SIZE, Z_SIZE))
    mips_detrend<- array(NA,c(X_SIZE,Y_SIZE, Z_SIZE))
    for (s in 1:Z_SIZE) {
      load(file=paste(outdir,'/mips/snap_mips_',s,'.rdata',sep=''))
      load(file=paste(outdir,'/mips/snap_mips_detrend_',s,'.rdata',sep=''))
      mips[,,s] <- snap_mips
      mips_detrend[,,s] <- snap_mips_detrend
    }
    mips[is.na(mips)]<-0
    mips_detrend[is.na(mips_detrend)]<-0
    f.write.nifti(mips,file=paste0(outdir,'/mips/time_mips.nii'), nii=TRUE )
    f.write.nifti(mips_detrend,file=paste0(outdir,'/mips/time_mips_detrend.nii'), nii=TRUE )
    
  }
  
  
  
  
  
  
  
  #f.write.nifti(snap_500,file=paste0(outdir,'/images/snapshot.nii'), nii=TRUE )
  
  ### Make big matrix to store all series voxels that are unmasked
  #Basis_number <- 100
  print(paste("Trying to allocate matrix of size (approx)", round(X_SIZE*Y_SIZE*Z_SIZE*Basis_number*8/(1024^3),2), " GB..."))
  big_mat <- matrix(NA,ssDM,Basis_number)
  Count <- 0
  
  active_mask[ is.na(active_mask) ] <- FALSE
  load(file=paste(outdir,'/D_Mask.rdata',sep=''))
  for (s in 1:Z_SIZE) {
    load(paste(outdir,'/coeff_mats/coeff_mat_',s,'.rdata',sep=''))
    InCount <- 0
    for (j in 1:Y_SIZE) {
      for (i in 1:X_SIZE) {
        InCount <- InCount + 1
        if(D_Mask[i,j,s] & active_mask[i,j,s]) {
          Count <- Count + 1
          big_mat[Count,] <- coeff_mat[InCount,]
        }
      }
    }
    print(paste("Loading slice", c(s), "of", Z_SIZE))
  }
  #& (active_mask[i,j,s] & !is.na(active_mask[i,j,s]))
  ### Scale the coeffient matrix for K-means
  #big_mat <- scale(big_mat)
  ## Remove big_mat for memory saving
  #rm(big_mat)
  
  #scales inplace, returns means and sd for later use
  mean_sd_from_unscaled<-scale_mat_inplace(big_mat)
  save(mean_sd_from_unscaled, file = paste0(outdir,"/params/scalingparams.rdata"))
  print("Matrix successfully loaded and rescaled.")
  
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
  
  if(file.exists(paste0(outdir,'/params/centers.rdata'))){
    load(file=paste0(outdir,'/params/centers.rdata'))
    comp<-dim(clustering)[1]
  }else{
    print("Doing clustering")
    ### Compute BIC over a range of K (here 50--100)
    max_clust<-50 #speed up by making this look at smaller
    BIC_val <- array(0, max_clust)
    TIME_STORE <- array(0, max_clust)
    
    
    
    # load any partial results already saved
    if(file.exists(paste0(outdir,'/params/Time_store.rdata'))){load(file=paste0(outdir,'/params/Time_store.rdata'))}
    if(file.exists(paste0(outdir,'/params/BIC_values.rdata'))){load(file=paste0(outdir,'/params/BIC_values.rdata'))}
    
    
    start_kk = max(which(BIC_val > 0),2)
    
    
    if(start_kk <= max_clust){ #check if havent gone too far
      print(paste("Any previous results saved, begining from number of clusters =", start_kk,". Search continues up to", max_clust))
      
      for (kk in (start_kk:max_clust)) {
        # Conduct timing while computing BIC values
        TIME <- proc.time()
        BIC_Temp<- cluster_BIC(big_mat,tkmeans(big_mat,kk,.9 ))
        #print(BIC_Temp)
        BIC_val[kk] <- BIC_Temp
        #print(BIC_val[kk])
        TIME_DUMMY <- proc.time() - TIME
        #print(TIME_STORE)
        TIME_STORE[kk] <- TIME_DUMMY[1] + TIME_DUMMY[2]
        #print(c(kk,BIC_val[kk]))
        # Save the results
        save(TIME_STORE,file=paste0(outdir,'/params/Time_store.rdata'))
        save(BIC_val,file=paste0(outdir,'/params/BIC_values.rdata'))
      }
    }
    
    
    
    ### Get the optimal K and computer clustering under optimal K
    n <- dim(big_mat)[1]
    m <- Basis_number
    BIC_val<-BIC_val[-1]
    neg_like <- BIC_val - log(n)*(m*(1:length(BIC_val)))
    log_like <- neg_like/(2)
    ave_log_like <- log_like/n
    names_vec <- (1:length(BIC_val))+1
    complexity_h <- shape_h <- (m*(1:length(BIC_val)))
    SHDATA <- cbind(names_vec, shape_h, complexity_h, ave_log_like)
    DD <- DDSE(SHDATA)
    comp <- as.numeric(attributes(DD)$model)
    
    # ### Cluster using the optimal value for K
    # #clustering <- tkmeans(x=big_mat,k=comp,alpha=.9,nstart=5,iter.max=20)
    n_starts  = 5
    # starts_list<-list()
    # BIC_list<-array(0, n_starts)
    #
    # print(paste("Running multiple (", n_starts, ") starts on optimal clustering number,", comp))
    # for(j in 1:n_starts){
    #
    #   starts_list[[j]] <- tkmeans_lowmem(big_mat, comp,.9,1,20) #only one start
    #
    #   BIC_list[j]<-BIC_lowmem(starts_list[[j]] ,big_mat)
    #
    # }
    
    clustering<-tkmeans(big_mat, comp, .9, nstart = n_starts)
    #save centres
    save(clustering,file=paste0(outdir,'/params/centers.rdata'))
    print("Clustering Done")
  }
  
  ### Function for allocating observations to cluster from a tkmeans clustering
  # tmeansClust <- function(fit,data) {
  #   apply(as.matrix(pdist2(data,t(fit$center))),1,which.min)
  # }
  
  print("Doing some more analysis on clustering results...")
  ### Get a clustering using the tkmeans clustering form variable "clustering"
  clustering_cluster <- nearest_cluster(big_mat,clustering)
  ## save cluster allocations
  save(clustering_cluster,file=paste(outdir,'/params/clustering.rdata',sep=''))
}

#reload params
load(file = paste0(outdir,"/params/scalingparams.rdata"))
load(file=paste0(outdir,'/params/centers.rdata'))
load(file=paste(outdir,'/params/clustering.rdata',sep=''))


#load(file=paste(indir,'/image_hold.rdata',sep=''))
#load(file=paste(indir,'/predicted_means.rdata',sep=''))
#load(file=paste(indir,'/correlation_matrix.rdata',sep=''))
#load(file=paste(indir,'/mean_image.rdata',sep=''))


### Make big matrix to store all series voxels that are unmasked
#Basis_number <- 100
print(paste("Trying to allocate matrix of size (approx)", round(X_SIZE*Y_SIZE*Z_SIZE*Basis_number*8/(1024^3),2), " GB..."))
big_mat <- matrix(NA,ssDM,Basis_number)
Count <- 0
for (s in 1:Z_SIZE) {
  load(paste(outdir,'/coeff_mats/coeff_mat_',s,'.rdata',sep=''))
  InCount <- 0
  for (j in 1:Y_SIZE) {
    for (i in 1:X_SIZE) {
      InCount <- InCount + 1
      if(D_Mask[i,j,s]) {
        Count <- Count + 1
        big_mat[Count,] <- coeff_mat[InCount,]
      }
    }
  }
  print(paste("Loading slice", c(s), "of", Z_SIZE))
}

### Scale the coeffient matrix for K-means
#big_mat <- scale(big_mat)
## Remove big_mat for memory saving
#rm(big_mat)

#scales inplace, returns means and sd for later use
mean_sd_from_unscaled<-scale_mat_inplace(big_mat)
save(mean_sd_from_unscaled, file = paste0(outdir,"/params/scalingparams.rdata"))
print("Matrix successfully loaded and rescaled. 2 ")


comp<-dim(clustering)[1]

clustering_cluster <- nearest_cluster(big_mat,clustering)
save(clustering_cluster, file=paste0(outdir,"/params/cluster_allocations.rdata"))

#print(dim(clustering_cluster))
### Produce Volume with cluster labels coordinates are (z,x,y)
image_hold <- array(NA,c(X_SIZE,Y_SIZE, Z_SIZE))
Count <- 0
for (s in 1:Z_SIZE) {
  for (j in 1:Y_SIZE) {
    for (i in 1:X_SIZE) {
      if (D_Mask[i,j,s]) {
        Count <- Count + 1
        image_hold[i,j,s] <- clustering_cluster[Count]
      }
    }
  }
}

save(image_hold,file=paste(outdir,'/clusters/image_hold.rdata',sep=''))
image_hold[is.na(image_hold)]<-0

f.write.nifti(image_hold,file=paste0(outdir,'/clusters/clusters.nii'), nii=TRUE, )

for(g in 1:comp){
  temp_mat<-array(0,dim(image_hold))
  temp_mat[image_hold==g]<-1
  f.write.nifti(temp_mat,file=paste0(outdir,'/clusters/cluster_', g ,'mask.nii'), nii=TRUE )
  
}



#image.plot(1:Y_SIZE,1:X_SIZE,image_hold[Z_SIZE,,])


#
# ### Obtain the cluster mean time series
# # Reload Graphing Parameters (These values are specific to F03 and F04)
file_number <- file_number.old
max_file <- max(file_number)
file_number <- (file_number-1)*Z_SIZE + s
Basis_number <- Basis_number
Basis <- create.bspline.basis(c(0,(max_file-1)*Z_SIZE+Z_SIZE),nbasis=Basis_number)
BS <- eval.basis(file_number,Basis)
# Compute the mean time series
centers <- clustering
#centers <- sweep(centers,2,attributes(big_mat)$'scaled:scale','*')
#centers <- sweep(centers,2,attributes(big_mat)$'scaled:center','+')

#reload scaling if need be
if(!exists("mean_sd_from_unscaled")){
  if(file.exists(paste0(outdir,"/params/scalingparams.rdata"))){
    load(paste0(outdir,"/params/scalingparams.rdata"))
  }
}

#centers <- sweep(centers,2,mean_sd_from_unscaled[1,],'*')
#centers <- sweep(centers,2,mean_sd_from_unscaled[2,],'+')
centers <-sweep(sweep(centers,2,mean_sd_from_unscaled[2,],'*'),2,mean_sd_from_unscaled[1,], '+')

pred_range <- eval.basis(seq(from=1,to=((max_file-1)*Z_SIZE+Z_SIZE),length.out=1000),Basis)

# Each row is a mean time series over the pred_range values
PRED <- matrix(NA,dim(centers)[1],dim(pred_range)[1])
for (j in 1:dim(centers)[1]) {
  PRED[j,] <- apply(pred_range,1,function(x) {sum(x*centers[j,])})
}


# The time average values of each series
MEAN_PRED <- rowMeans(PRED)
save(PRED,file=paste(outdir,'/images/predicted_means.rdata',sep=''))




big_mat<-sweep(sweep(big_mat,2,mean_sd_from_unscaled[2,],'*'),2,mean_sd_from_unscaled[1,], '+')

for (p in 1:dim(centers)[1]) {
  
  whom<-which(clustering_cluster==p)
  
  if(length(whom)>10){
    n_samp<-min(100, length(whom))
    
    pred_to_plot <- array(NA,c(n_samp,dim(pred_range)[1]))
    sample_to_plot<- array(NA,c(n_samp,dim(big_mat)[2]))
    
    sample_to_plot<-big_mat[sample(whom,n_samp),]
    
    for(k in 1:n_samp){
      pred_to_plot[k,] <- apply(pred_range,1,function(x) {sum(x*sample_to_plot[k,])})
    }
    
    pred_to_plot2<-data.frame(t(pred_to_plot))
    pred_to_plot2$time<-seq(from=0,to=0.2*(TO-FROM), length.out=dim(pred_range)[1])
    pred_to_plot3<-melt(pred_to_plot2, id='time')
    pdf(paste0(outdir,'/cluster_timeseries_plots/cluster_', p,'.pdf',sep=''),paper='a4r')
    plot<-ggplot(data=pred_to_plot3,aes(y=value,x=time,group=variable))+geom_line(size=0.1)+scale_y_continuous("Intensity", limits=ylimits1)+ggtitle(paste("Cluster ",p ))+scale_x_continuous("Time(s)")+theme_bw()
    print(plot)
    dev.off()
  }
  
}



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
    #what does this do
    counter_corr = 0
    while( xor(CORR[c1,c2] > 1, CORR[c1,c2] < -1) & counter_corr < 100 ) {
      if (CORR[c1,c2] > 1) {INTEGRAL <- INTEGRAL - PREC*INT_OLD}
      if (CORR[c1,c2] < -1) {INTEGRAL <- INTEGRAL + PREC*INT_OLD}
      CORR[c1,c2] <- INTEGRAL/sqrt(COVAR[c1]*COVAR[c2])
      counter_corr=counter_corr+1
    }
  }
}

save(CORR,file=paste(outdir,'/images/correlation_matrix.rdata',sep=''))

### Create a Mean value image (using the predicted mean values from MEAN_PRED)
which_clust <- 6
mean_image <- array(NA,c(X_SIZE,Y_SIZE, Z_SIZE))
count <- 0
for (s in 1:Z_SIZE)
{
  for (j in 1:Y_SIZE) {
    for (i in 1:X_SIZE) {
      if (D_Mask[i,j,s]) {
        count <- count + 1
        mean_image[i,j,s] <- MEAN_PRED[clustering_cluster[count]]
      }
    }
  }
}

### Save All Intermediate Results

save(mean_image,file=paste(outdir,'/images/mean_image.rdata',sep=''))
mean_image[is.na(mean_image)]<-0
f.write.nifti(mean_image,file=paste0(outdir,'/images/mean_image.nii'), nii=TRUE )



print("Analysis done and saved. Producing and saving plots.")
# Graphing --------------------------------------------

# ### First set of graphs
# # Graph clustering on every 5th slice
# for (slice in seq(5,Z_SIZE,by=5)) {
#   pdf(paste(indir,'/Clustering_Slice_',slice,'.pdf',sep=''),paper='a4r')
#   ### Plot Clustering Image
#   image_hold <- array(NA,c(X_SIZE,Y_SIZE, Z_SIZE))
#   count <- 0
#   for (s in 1:Z_SIZE) {
#     for (j in 1:Y_SIZE) {
#       for (i in 1:X_SIZE) {
#         if (D_Mask[i,j,s]) {
#           count <- count + 1
#           image_hold[i,j,s] <- clustering_cluster[count]
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
#   mean_image <- array(NA,c(X_SIZE,Y_SIZE, Z_SIZE))
#   count <- 0
#   for (s in 1:Z_SIZE)
#   {
#     for (j in 1:Y_SIZE) {
#       for (i in 1:X_SIZE) {
#         if (D_Mask[i,j,s]) {
#           count <- count + 1
#           mean_image[i,j,s] <- MEAN_PRED[clustering_cluster[count]]
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
#       for (s in 1:Z_SIZE) {
#         for (j in 1:Y_SIZE) {
#           for (i in 1:X_SIZE) {
#             if (D_Mask[i,j,s]) {
#               count <- count + 1
#               if (clustering_cluster[count]==cc & s==75) {
#                 points(j,i,pch=15,cex=1,col='green')
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
#       for (s in 1:Z_SIZE) {
#         for (j in 1:Y_SIZE) {
#           for (i in 1:X_SIZE) {
#             if (D_Mask[i,j,s]) {
#               count <- count + 1
#               if (clustering_cluster[count]==cc & s==50) {
#                 points(j,i,pch=15,cex=1,col='green')
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
#       for (s in 1:Z_SIZE) {
#         for (j in 1:Y_SIZE) {
#           for (i in 1:X_SIZE) {
#             if (D_Mask[i,j,s]) {
#               count <- count + 1
#               if (clustering_cluster[count]==cc & s==100) {
#                 points(j,i,pch=15,cex=1,col='green')
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

# Graph the Mean functions
for (cc in 1:comp) {
  pdf(paste(outdir,'/cluster_functions/Cluster_Mean_Function_',cc,'.pdf',sep=''),paper='a4r')
  #ggplot(data=pred.df.flat)+geom_line(aes(y=Signal, x=Time))+facet_wrap(~variable)+theme_bw()
  
  plot(seq(from=1,to=((max_file-1)*Z_SIZE+Z_SIZE),length.out=1000),PRED[cc,],type='l',xlab='Time',ylab='signal',main=cc)
  dev.off()
}

##R version issue here
pred.df<-data.frame(t(PRED))
names(pred.df) <- paste0("Cluster ", 1:comp)
pred.df$Time<-seq(0,dim(PRED)[2]-1,1)
pred.df.flat<-melt(pred.df, id='Time')
pred.df.flat$Signal<-pred.df.flat$value

pdf(paste(outdir,'/cluster_functions/All_Cluster_Mean_Function.pdf',sep=''),paper='a4')
ggplot(data=pred.df.flat)+geom_line(aes(y=Signal, x=Time))+facet_wrap(~variable)+theme_bw()
dev.off()

# Plot the Correlation matrix
pdf(paste(outdir,'/images/Correlation_matrix.pdf',sep=''),paper='a4r')
image.plot(1:comp,1:comp,CORR)
dev.off()

# Hierachical Clustering ------------------------------
# # Make a distance metric
DIST <- as.dist(1-t(CORR))
HCLUST <- hclust(DIST,method='average')
pdf(paste(outdir,'/cluster_functions/Cluster_dendrogram.pdf',sep=''),paper='a4r')
plot(HCLUST)
dev.off()

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
# for (j in 1:max(CUT)) {
#   pdf(paste(indir,'/HCLUST_',j,'.pdf',sep=''),paper='a4r')
#   plot(c(1,((max_file-1)*Z_SIZE+Z_SIZE)),c(min(PRED),max(PRED)),
#        type='n',xlab='Time',ylab='signal',main=j)
#   for (s in 1:comp) {
#     if (CUT[s]==j) {
#       lines(seq(1,((max_file-1)*Z_SIZE+Z_SIZE),1000),PRED[s,],col=tim.colors(comp)[s])
#     }
#   }
#   dev.off()
# }
#
# Plot Frequency Histogram
pdf(paste0(outdir,'/cluster_functions/Frequency_of_clusters.pdf',sep=''),paper='a4r')
plot(table(clustering_cluster),xlab='cluster',ylab='Frequency')
dev.off()

# # Plot K means by HCLUST reference
# pdf(paste(indir,'/Cluster_by_HCLUST.pdf',sep=''),paper='a4r')
# plot(CUT,col=tim.colors(comp),xlim=c(0,comp+1),ylim=c(0,8),xlab='Cluster',ylab='HCLUST')
# grid()
# dev.off()


# # Graph clustering based on cuts on every 5th slice
# for (slice in seq(5,Z_SIZE,by=5)) {
#   pdf(paste(indir,'/CUT_Slice_',slice,'.pdf',sep=''),paper='a4r')
#   ### Plot Clustering Image
#   image_hold <- array(NA,c(X_SIZE,Y_SIZE, Z_SIZE))
#   count <- 0
#   for (s in 1:Z_SIZE) {
#     for (j in 1:Y_SIZE) {
#       for (i in 1:X_SIZE) {
#         if (D_Mask[i,j,s]) {
#           count <- count + 1
#           image_hold[i,j,s] <- CUT[clustering_cluster[count]]
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
print("04-Rstats02.R complete.")
