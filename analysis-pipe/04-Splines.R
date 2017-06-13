
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
Z_START <- 250

X_SIZE <- 640#640
Y_SIZE <- 130 #130
Z_SIZE <- 10#300

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


### Fit Splines to the time series in each of the slice files (s loops over slices)
for (s in 1:Z_SIZE) {
  
  print(mem_used())
  rm(list = ls()[-c(which(ls()=="outdir"), which(ls()=="s"))])
  gc()
  load(file=paste(outdir,'/MASK_hdr.rdata',sep=''))
  load(file=paste(outdir,'/settings.rdata',sep=''))
  load(file=paste(outdir,'/ssDM.rdata',sep=''))
  load(file=paste(outdir,'/D_Mask.rdata',sep=''))
  print(mem_used())
  print(c('Doing ', s))
  
  ### Get already existing outputs
  Completed_output <- list.files(path=paste0(outdir, '/coeff_mats/'), pattern='coeff_mat', full.names=FALSE)
  #Completed_output <- Completed_output[grep('coeff_mat',Completed_output)]
  ### Skip things already in Completed_output
  
  if (paste('coeff_mat_',s,'.rdata',sep='')%in%Completed_output){
    print(paste('coeff_mat_', s, '.rdata', ' exists', sep=''))
  }
  else{

    file_number <- (file_number.old-1)*Z_SIZE +s
    
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
    
    if(!file.exists(paste(outdir,'/MASK_active.rdata',sep=''))){
      active_mask <- array(NA,c(X_SIZE,Y_SIZE, Z_SIZE))
    }else{
      load(file=paste(outdir,'/MASK_active.rdata',sep=''))
    }
    
    for (j in 1:Y_SIZE) {
      for (i in 1:X_SIZE) {
        count <- count + 1
        active_mask[i,j,s] <-  (max(mat_odds_filt[count,],na.rm=T) - min(mat_odds_filt[count,],na.rm=T))
        snap_mips_detrend[i,j] <- max(mat_odds_detr[count,],na.rm=T)
        snap_mips_filt[i,j] <- max(mat_odds_filt[count,],na.rm=T)
        snap_mips[i,j] <- max(mat_odds[count,],na.rm=T)
      }
    }
  
    print(mem_used())
    save(active_mask,file=paste(outdir,'/MASK_active.rdata',sep=''))
    save(snap_mips, file=paste(outdir,'/mips/snap_mips_',s,'.rdata',sep=''))
    save(snap_mips_detrend, file=paste(outdir,'/mips/snap_mips_detrend_',s,'.rdata',sep=''))
    save(snap_mips_filt, file=paste(outdir,'/mips/snap_mips_filt_',s,'.rdata',sep=''))
    pred_splines <- matrix(NA,X_SIZE*Y_SIZE,length(file_number))
    
    ### Declare number of bases to use
    
    Basis <- create.bspline.basis(c(0,(max_file-1)*Z_SIZE+Z_SIZE),nbasis=Basis_number)
    BS <- eval.basis(file_number,Basis)
    # Each row is a mean time series over the pred_range values
    
    ### Fit B-spline to all time series
    print(c('Spline fitting',s))
    print(length(file_number))
    
    which_odds<-(1:length(file_number))%%2 #SPLIT
    
    FD <- smooth.basisPar(file_number[!!which_odds],t(mat_odds_filt),Basis)

    coeff_mat <- as.matrix(t(FD$fd$coefs))
    
    for (i in 1:(X_SIZE*Y_SIZE)){
      pred_splines[i,] <- apply(BS ,1,function(x) {sum(x*coeff_mat[i,])})
    }
    
    pred_splines2<-data.frame(t(pred_splines))
    pred_splines2$time<-seq(from=0,to=0.2*(TO-FROM), length.out=(length(file_number)))
    
    dir.create(paste0(outdir,'/individual_plots/'))
    
    random_splines<-array(0,20)
    i=1
    j=1
    
    while(i<21 & j<1000){
      temp<-sample(X_SIZE*Y_SIZE,1)-1
      if(active_mask[temp%%X_SIZE+1, temp%/%X_SIZE+1,s]){
        random_splines[i]<-temp
        i<-i+1
      }
      j<-j+1
    }
    
    print('Individual plots')


      for(i in 1:20){
        if(random_splines[i]!=0){

        p<-paste0("X",random_splines[i])
        plot<-ggplot()+geom_point(aes(x=seq(from=0,to=0.2*(TO-FROM), length.out=(length(file_number)/2)), y=mat_odds_filt[random_splines[i],]))+geom_line(data=pred_splines2,aes(y=get(p),x=time), color='blue') +scale_y_continuous("Intensity", limits=ylimits1)+scale_x_continuous("Time (s)")
        x_loc <- random_splines[i]%%X_SIZE+1
        y_loc <-random_splines[i]%/%X_SIZE+1
        pdf(paste0(outdir,'/individual_plots/point_', x_loc+X_START, '_', y_loc+Y_START, '_', s+Z_START, '.pdf',sep=''),paper='a4r')
        print(plot)
        dev.off()
        }
      }

    ### Save the results in the format coeff_mat_(SLICE NUMBER).rdata
    print(paste('Saving as: ', outdir, 'coeff_mats/coeff_mat_', s, '.rdata', sep=''))
    save(coeff_mat,file=paste(outdir,'/coeff_mats/coeff_mat_',s,'.rdata',sep=''))
    print(mem_used())
  }
  
  
}

#save mips when finish all slices
if(!file.exists(file=paste0(outdir,'/mips/time_mips.nii'))){
  mips <- array(NA,c(X_SIZE,Y_SIZE, Z_SIZE))
  mips_detrend<- array(NA,c(X_SIZE,Y_SIZE, Z_SIZE))
  mips_filt<- array(NA,c(X_SIZE,Y_SIZE, Z_SIZE))
  for (s in 1:Z_SIZE) {
    load(file=paste(outdir,'/mips/snap_mips_',s,'.rdata',sep=''))
    load(file=paste(outdir,'/mips/snap_mips_detrend_',s,'.rdata',sep=''))
    load(file=paste(outdir,'/mips/snap_mips_filt_',s,'.rdata',sep=''))
    mips[,,s] <- snap_mips
    mips_detrend[,,s] <- snap_mips_detrend
    mips_filt[,,s] <- snap_mips_detrend
  }
  mips[is.na(mips)]<-0
  mips_detrend[is.na(mips_detrend)]<-0
  f.write.nifti(mips,file=paste0(outdir,'/time_mips.nii'), nii=TRUE, L=header)
  f.write.nifti(mips_detrend,file=paste0(outdir,'/time_mips_detrended.nii'), nii=TRUE,L=header )
  f.write.nifti(mips_filt,file=paste0(outdir,'/time_mips_filter_detrended.nii'), nii=TRUE,L=header )
}

print("04-Rstats02.R complete.")
