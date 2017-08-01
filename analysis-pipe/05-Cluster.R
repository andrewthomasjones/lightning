list.of.packages <- c('lattice', 'AnalyzeFMRI', 'ggplot2', 'reshape2', 'MASS', 'abind', 'fda', 'fields', 'speedglm','pracma', 'tclust', 'signal', 'capushe', 'pryr', 'lowmemtkmeans')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

### Load necessary libraries (all CRAN libraries can be acquired using install.packages)
#library(compiler)
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
indir_new <- args[1]
outdir_new <- args[2]
FROM <- as.numeric(args[3])
TO <- as.numeric(args[4])
mask_fn <- args[5]


load(file=paste(indir_new,'/settings.rdata',sep=''))

outdir<-outdir_new
indir<-indir_new

max_clust<-30 #speed up by making this look at smaller
cut_p<-0.04

# X_START <- 0
# Y_START <- 0
# Z_START <- 0
# 
# X_SIZE <- 640#640
# Y_SIZE <- 130 #130
# Z_SIZE <- 300#300


### These lines are kind of data dependent. I've determined what works for F11, F03, and F04.
## FROM and TO determine the start and end points of the trimming


print(c('1indir', indir))
print(c('2outdir', outdir))
print(c('3FROM', FROM))
print(c('4TO', TO))
print(c('5mask', mask_fn))

print(outdir)
if(!dir.exists(outdir)){
  dir.create(outdir, recursive=T)
}

dir.create(paste0(outdir,'/cluster_timeseries_plots/'))
dir.create(paste0(outdir,'/clusters/'))
dir.create(paste0(outdir,'/cluster_functions/'))

load(file=paste(indir,'/ssDM.rdata',sep=''))
load(file=paste(indir,'/MASK_hdr.rdata',sep=''))
load(file=paste(indir,'/D_Mask.rdata',sep=''))
load(file=paste(indir,'/MASK_active.rdata',sep=''))

if(!file.exists(paste0(outdir,'/clustering.rdata'))){
 
  cutoff<-quantile(active_mask[D_Mask], probs=1-cut_p, na.rm=T) 
  active_mask2 <- active_mask > cutoff
  active_mask2[ is.na(active_mask2) ] <- FALSE
  save(active_mask2,file=paste(outdir,'/MASK_active_adjusted.rdata',sep=''))
  ssDM_active <-sum(active_mask2 & D_Mask)
  save(ssDM_active,file=paste(outdir,'/ssDM_active.rdata', sep=''))
  print(paste("Total voxels ", ssDM_active))
  
  print(paste("Trying to allocate matrix of size (approx)", round(ssDM_active*Basis_number*8/(1024^3),2), " GB..."))
  big_mat <- matrix(NA,ssDM_active,Basis_number)
  
 
  Count <- 0
  for (s in 1:Z_SIZE) {
    load(paste(indir,'/coeff_mats/coeff_mat_',s,'.rdata',sep=''))
    InCount <- 0
    for (j in 1:Y_SIZE) {
      for (i in 1:X_SIZE) {
        InCount <- InCount + 1
        if(D_Mask[i,j,s] & active_mask2[i,j,s]) {
          Count <- Count + 1
          big_mat[Count,] <- coeff_mat[InCount,]
        }
      }
    }
    print(paste("Loading slice", c(s), "of", Z_SIZE))
  }

  #scales inplace, returns means and sd for later use
  mean_sd_from_unscaled<-scale_mat_inplace(big_mat)
  save(mean_sd_from_unscaled, file = paste0(outdir,"/scalingparams.rdata"))
  print("Matrix successfully loaded and rescaled.")
  
  # TRIMMED K-MEANS CLUSTERING ------------------------

  if(file.exists(paste0(outdir,'/centers.rdata'))){
    load(file=paste0(outdir,'/centers.rdata'))
    comp<-dim(clustering)[1]
  }else{
    print("Doing clustering")
    ### Compute BIC over a range of K (here 50--100)
    
    BIC_val <- array(0, max_clust)
    TIME_STORE <- array(0, max_clust)
    
    #load(file=paste(outdir,'/settings.rdata',sep=''))
    
    # load any partial results already saved
    if(file.exists(paste0(outdir,'/Time_store.rdata'))){load(file=paste0(outdir,'/Time_store.rdata'))}
    if(file.exists(paste0(outdir,'/BIC_values.rdata'))){load(file=paste0(outdir,'/BIC_values.rdata'))}
    
    
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
        print(c(kk,BIC_val[kk]))
        # Save the results
        save(TIME_STORE,file=paste0(outdir,'/Time_store.rdata'))
        save(BIC_val,file=paste0(outdir,'/BIC_values.rdata'))
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
    
    n_starts  = 5
 
    
    clustering<-tkmeans(big_mat, comp, .9, nstart = n_starts)
    #save centres
    save(clustering,file=paste0(outdir,'/centers.rdata'))
    print("Clustering Done")
  }
  
  print("Doing some more analysis on clustering results...")
  ### Get a clustering using the tkmeans clustering form variable "clustering"
  clustering_cluster <- nearest_cluster(big_mat,clustering)
  ## save cluster allocations
  save(clustering_cluster,file=paste(outdir,'/clustering.rdata',sep=''))
}

#reload params
load(file = paste0(outdir,"/scalingparams.rdata"))
load(file=paste0(outdir,'/centers.rdata'))
load(file=paste(outdir,'/clustering.rdata',sep=''))#????

load(file=paste(indir,'/D_Mask.rdata',sep=''))
load(file=paste(outdir,'/MASK_active_adjusted.rdata',sep=''))
load(file=paste(outdir,'/ssDM_active.rdata', sep=''))

print(paste("Trying to allocate matrix of size (approx)", round(ssDM_active*Basis_number*8/(1024^3),2), " GB..."))
big_mat <- matrix(NA,ssDM_active,Basis_number)
Count <- 0

#load(file=paste(outdir,'/D_Mask.rdata',sep=''))
for (s in 1:Z_SIZE) {
  load(paste(indir,'/coeff_mats/coeff_mat_',s,'.rdata',sep=''))
  InCount <- 0
  for (j in 1:Y_SIZE) {
    for (i in 1:X_SIZE) {
      InCount <- InCount + 1
      if(D_Mask[i,j,s] & active_mask2[i,j,s]) {
        Count <- Count + 1
        big_mat[Count,] <- coeff_mat[InCount,]
      }
    }
  }
  print(paste("Loading slice", c(s), "of", Z_SIZE))
}

#scales inplace, returns means and sd for later use
mean_sd_from_unscaled<-scale_mat_inplace(big_mat)
save(mean_sd_from_unscaled, file = paste0(outdir,"/scalingparams.rdata"))
print("Matrix successfully loaded and rescaled. 2 ")

comp<-dim(clustering)[1]

clustering_cluster <- nearest_cluster(big_mat,clustering)
save(clustering_cluster, file=paste0(outdir,"/cluster_allocations.rdata"))

### Produce Volume with cluster labels coordinates are (z,x,y)
image_hold <- array(NA,c(X_SIZE,Y_SIZE, Z_SIZE))
Count <- 0
for (s in 1:Z_SIZE) {
  for (j in 1:Y_SIZE) {
    for (i in 1:X_SIZE) {
      if (D_Mask[i,j,s] & active_mask2[i,j,s]) {
        Count <- Count + 1
        image_hold[i,j,s] <- clustering_cluster[Count]
      }
    }
  }
}

save(image_hold,file=paste(outdir,'/clusters/image_hold.rdata',sep=''))
image_hold[is.na(image_hold)]<-0

f.write.nifti(image_hold,file=paste0(outdir,'/clusters/clusters.nii'), nii=TRUE, L=header)

active_mask2[is.na(active_mask2)]<-0
f.write.nifti(active_mask2,file=paste0(outdir,'/active_mask_adjusted.nii'), nii=TRUE, L=header )

for(g in 1:comp){
  temp_mat<-array(0,dim(image_hold))
  temp_mat[image_hold==g]<-1
  f.write.nifti(temp_mat,file=paste0(outdir,'/clusters/cluster_', g ,'mask.nii'), nii=TRUE, L=header )
}


# ### Obtain the cluster mean time series
# # Reload Graphing Parameters (These values are specific to F03 and F04)

file_number <- (file_number.old -1)*Z_SIZE + Z_SIZE
Basis <- create.bspline.basis(c(0,(max_file-1)*Z_SIZE+Z_SIZE),nbasis=Basis_number)
BS <- eval.basis(file_number,Basis)

# Compute the mean time series
load(file=paste0(outdir,'/centers.rdata'))
load(file = paste0(outdir,"/scalingparams.rdata"))
centers <- clustering
centers <-sweep(sweep(centers,2,mean_sd_from_unscaled[2,],'*'),2,mean_sd_from_unscaled[1,], '+')
clustering_cluster <- nearest_cluster(big_mat,centers)


#reload scaling if need be
if(!exists("mean_sd_from_unscaled")){
  if(file.exists(paste0(outdir,"/scalingparams.rdata"))){
    load(paste0(outdir,"/scalingparams.rdata"))
  }
}

pred_range <- eval.basis(seq(from=1,to=((max_file-1)*Z_SIZE+Z_SIZE),length.out=1000),Basis)

# Each row is a mean time series over the pred_range values
PRED <- matrix(NA,dim(centers)[1],dim(pred_range)[1])
for (j in 1:dim(centers)[1]) {
  PRED[j,] <- apply(pred_range,1,function(x) {sum(x*centers[j,])})
}


# The time average values of each series
MEAN_PRED <- rowMeans(PRED)


for (p in 1:dim(centers)[1]) {
  
  whom<-which(clustering_cluster==p)
  
  if(length(whom)>10){
    n_samp<-min(50, length(whom))
    
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
for (cc in 1:dim(centers)[1]) {
  S1 <- function(x) Spline_function(x,cc)
  S2 <- function(x) Spline_function(x,cc)
  S_Prod <- function(x) S1(x)*S2(x)
  INTEGRAL <- quadv(S_Prod,1,(max_file-1)*Z_SIZE+Z_SIZE)$Q
  COVAR[cc] <- INTEGRAL
}
comp<-dim(centers)[1]
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


#also graph here  with gggplot2
##R version issue here
pred.df<-data.frame(t(PRED))
names(pred.df) <- paste0("Cluster ", 1:comp)
pred.df$Time<-seq(from=0,to=(TO-FROM)/5,length.out=1000)
pred.df.flat<-melt(pred.df, id='Time')
pred.df.flat$Signal<-pred.df.flat$value

# Graph the Mean functions
for (cc in 1:comp) {
  pdf(paste(outdir,'/cluster_functions/Cluster_Mean_Function_',cc,'.pdf',sep=''),paper='a4r')
  p1<-ggplot(data=subset(pred.df.flat, variable==paste0("Cluster ",cc)))+geom_line(aes(y=Signal, x=Time))+theme_bw()+scale_y_continuous("Intensity", limits=ylimits1)+scale_x_continuous("Time (s)")
  print(p1)
  #plot(seq(from=1,to=(TO-FROM)/5,length.out=1000),PRED[cc,],type='l',xlab='Time',ylab='signal',main=cc)
  
  #plot(seq(from=1,to=((max_file-1)*Z_SIZE+Z_SIZE),length.out=1000),PRED[cc,],type='l',xlab='Time',ylab='signal',main=cc)
  dev.off()
}


pdf(paste(outdir,'/cluster_functions/All_Cluster_Mean_Function.pdf',sep=''),paper='a4')
ggplot(data=pred.df.flat)+geom_line(aes(y=Signal, x=Time))+facet_wrap(~variable)+theme_bw()+scale_y_continuous("Intensity", limits=ylimits1)+scale_x_continuous("Time (s)")
dev.off()

# Plot the Correlation matrix
pdf(paste(outdir,'/Correlation_matrix.pdf',sep=''),paper='a4r')
image.plot(1:comp,1:comp,CORR)
dev.off()

# Plot Frequency Histogram
pdf(paste0(outdir,'/cluster_functions/Frequency_of_clusters.pdf',sep=''),paper='a4r')
plot(table(clustering_cluster),xlab='cluster',ylab='Frequency')
dev.off()

# Hierachical Clustering ------------------------------
# # Make a distance metric
DIST <- as.dist(1-t(CORR))
HCLUST <- hclust(DIST,method='average')

if(comp>2){
  pdf(paste(outdir,'/cluster_functions/Cluster_dendrogram.pdf',sep=''),paper='a4r')
  plot(HCLUST)
  dev.off()
}

print("05-Cluster.R complete.")
