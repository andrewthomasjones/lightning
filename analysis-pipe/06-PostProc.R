list.of.packages <- c('matrixStats', 'lattice', 'AnalyzeFMRI', 'ggplot2', 'reshape2', 'MASS', 'abind', 'fda', 'fields', 'speedglm','pracma', 'tclust', 'signal', 'capushe', 'pryr', 'lowmemtkmeans')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos="https://cloud.r-project.org")

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
library(matrixStats)
### Parse CMD Line Arguments
args <- commandArgs(trailingOnly = TRUE)
indir_new <- args[1] # $input = "04-Splines"; 
outdir_new <- args[2]
FROM <- as.numeric(args[3])
TO <- as.numeric(args[4])
mask_fn <- args[5]
indir_new_2 <- args[6] # $input2 = "05-Clusters";
indir_new_3 <- args[7] # $input3 = "03-chunk02";

load(file=paste(indir_new,'/settings.rdata',sep=''))

outdir<-outdir_new
indir<-indir_new
indir2<-indir_new_2
indir3<-indir_new_3

max_clust<-30 #speed up by making this look at smaller
cut_p<-0.04
merge_height<-0.45  

print(c('1indir', indir))
print(c('2outdir', outdir))
print(c('3FROM', FROM))
print(c('4TO', TO))
print(c('5mask', mask_fn))
print(c('6indir2', indir_new_2))
print(c('7indir3', indir_new_3))

print(outdir)
if(!dir.exists(outdir)){
  dir.create(outdir, recursive=T)
}

#create subfolders
dir.create(paste0(outdir,'/clusters/'))
dir.create(paste0(outdir,'/time_series/'))

#reload misc
load(file=paste(indir,'/ssDM.rdata',sep=''))
load(file=paste(indir,'/MASK_hdr.rdata',sep=''))
load(file=paste(indir,'/D_Mask.rdata',sep=''))
load(file=paste(indir,'/MASK_active.rdata',sep=''))

load(file = paste0(indir2,"/scalingparams.rdata"))
load(file=paste0(indir2,'/centers.rdata'))
load(file=paste(indir2,'/clustering.rdata',sep=''))#????

load(file=paste(indir,'/D_Mask.rdata',sep=''))
load(file=paste(indir2,'/MASK_active_adjusted.rdata',sep=''))
load(file=paste(indir2,'/ssDM_active.rdata', sep=''))
load(file=paste(indir2,'/clusters/image_hold.rdata', sep=''))
centers <- clustering
centers <-sweep(sweep(centers,2,mean_sd_from_unscaled[2,],'*'),2,mean_sd_from_unscaled[1,], '+')


###########################################
if(!file.exists(paste(outdir,'/clusters/image_hold_merge.rdata',sep=''))){

  
  # ### CLUSTER MERGING
  
  
  file_number <- (file_number.old -1)*Z_SIZE + Z_SIZE
  Basis <- create.bspline.basis(c(0,(max_file-1)*Z_SIZE+Z_SIZE),nbasis=Basis_number)
  BS <- eval.basis(file_number,Basis)
  
  pred_range <- eval.basis(seq(from=1,to=((max_file-1)*Z_SIZE+Z_SIZE),length.out=1000),Basis)
  
  # Each row is a mean time series over the pred_range values
  PRED <- matrix(NA,dim(centers)[1],dim(pred_range)[1])
  for (j in 1:dim(centers)[1]) {
    PRED[j,] <- apply(pred_range,1,function(x) {sum(x*centers[j,])})
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
  
  
    # Hierachical Clustering ------------------------------
    # # Make a distance metric
    DIST <- as.dist(1-t(CORR))
    HCLUST <- hclust(DIST,method='average')
  
  
  print(HCLUST$merge)
  print(HCLUST$height)
  save(HCLUST, file=paste(outdir,'/clusters/tree.rdata',sep=''))
  
  
  merge_list<-array(0, comp-1)
  for (i in 1:(comp-1)){
    
    if (HCLUST$height[i] < merge_height){
    print(paste("Merge", i,": "))  
      if (HCLUST$merge[i,1]<0 & HCLUST$merge[i,2]<0){
        image_hold[image_hold== -HCLUST$merge[i,2]]<- -HCLUST$merge[i,1]
        merge_list[i]<- -HCLUST$merge[i,1]
        print(paste(-HCLUST$merge[i,2], "into",-HCLUST$merge[i,1]))
      }
      if (HCLUST$merge[i,1]<0 & HCLUST$merge[i,2]>0){
        image_hold[image_hold== -HCLUST$merge[i,1]]<- merge_list[HCLUST$merge[i,2]]
        merge_list[i]<-merge_list[HCLUST$merge[i,2]]
        print(paste(-HCLUST$merge[i,1],  "into", merge_list[HCLUST$merge[i,2]]))
      }
      if (HCLUST$merge[i,1]>0 & HCLUST$merge[i,2]>0){
        image_hold[image_hold== merge_list[HCLUST$merge[i,1]]]<- merge_list[HCLUST$merge[i,2]]
        merge_list[i]<-merge_list[HCLUST$merge[i,2]]
        print(paste(merge_list[HCLUST$merge[i,1]],  "into", merge_list[HCLUST$merge[i,2]])) 
      }
      
      
    }
    
  }
  
  tab<-table(image_hold)
  cl<-length(tab)
  
  
  
  print("Saving merged clusters")  
  image_hold[is.na(image_hold)]<-0
  f.write.nifti(image_hold,file=paste0(outdir,'/clusters/clusters_merge_orignal_numbering.nii'), nii=TRUE, L=header)
  
  image_hold2<-image_hold
  print(cl)
  
  for(g in 1:cl){
    image_hold2[image_hold2==as.numeric(names(tab))[g]]<-(g)
  }
  
  tab<-table(image_hold2)
  cl<-length(tab)
  print(tab)
  
  for(g in 1:cl){ #need to skip background
    temp_mat<-array(0,dim(image_hold2))
    temp_mat[image_hold2==as.numeric(names(tab))[g]]<-1
    if((g-1)>0){
      f.write.nifti(temp_mat,file=paste0(outdir,'/clusters/cluster_', g-1 ,'mask.nii'), nii=TRUE, L=header )
    }
  }
  
  # Plot Frequency Histogram
  pdf(paste0(outdir,'/clusters/Frequency_of_clusters.pdf',sep=''),paper='a4r')
  plot(table(image_hold2),xlab='cluster',ylab='Frequency')
  dev.off()

  
  f.write.nifti(image_hold2,file=paste0(outdir,'/clusters/clusters_merge.nii'), nii=TRUE, L=header)
  tab<-table(image_hold2)
  cl<-length(tab)
  print(tab)
  save(image_hold2,file=paste(outdir,'/clusters/image_hold_merge.rdata',sep=''))
  
}else{
  load(file=paste(outdir,'/clusters/image_hold_merge.rdata',sep=''))
  load(file=paste(outdir,'/clusters/tree.rdata',sep=''))
  
}

###########################################
if(!file.exists(paste(outdir,'/time_series/bg_data.rdata',sep=''))){
  if(!file.exists(paste(outdir,'/time_series/temp_data.rdata',sep=''))){
    print("Reloading raw time series data") 
    tab<-table(image_hold2)
    cl<-length(tab)
    print(tab)
    active_vox<-sum(image_hold2>0)
    ### Store data into matrix instead of array
    count <- 0
    count_bg <- 0
    print(paste(active_vox, N))
    print(paste("Trying to allocate matrix of size (approx)", round(active_vox*(N+4)*8/(1024^3),2), " GB..."))
    full_mat <- matrix(NA,active_vox,N)
    bg_mat <- matrix(0,2,N)
    meta_mat <- matrix(NA,active_vox,4)
    k=1
    
  }else{
    print("Reloading raw time series data") 
    print(paste(outdir,'/time_series/temp_data.rdata',sep='')) 
    load(paste(outdir,'/time_series/temp_data.rdata',sep=''))
  }
  
  for (s in k:Z_SIZE) {
    
    sink(file=paste(outdir,'/time_series/prog.txt',sep=''))
    cat(s)
    cat("\n")
    sink()
    print(paste("Loading slice",s,"of", Z_SIZE)) 
    file_number <- (file_number.old-1)*Z_SIZE +s+Z_START
      
      full_array <- list()
      
      ### Read in Z slices
      for (j in 1:N) {
        file_name <- substring(file_list[j],1,12)
        nii_name <- paste(indir3, '/', file_name,'.nii',sep='')
        #print(c('reading',nii_name, s+Z_START,j,round(j/N),1))
        NIFTI <- f.read.nifti.slice(nii_name,s+Z_START,1)
        full_array[[j]] <- NIFTI
      }
      
      ### Combine all of the time slices
      full_array <- abind(full_array,along=3)
      
      s<-s+Z_START
      for (j in 1:Y_SIZE) {
        for (i in 1:X_SIZE) {
          
          if(image_hold2[i,j,s]>0){
            count <- count + 1
            full_mat[count,] <- full_array[i+X_START,j+Y_START,]
            #x,y,z, cluster, 
            meta_mat[count,] <- c(i,j,s,image_hold2[i,j,s])
          }
          
          if(D_Mask[i,j,s]==T){
            count_bg <- count_bg + 1
            temp1<-full_array[i+X_START,j+Y_START,]
            temp1[is.nan(temp1)]<-0
            
            bg_mat[1,]<-bg_mat[1,]+temp1
            bg_mat[2,]<-bg_mat[2,]+(temp1)^2
          }
          #print(c('storing',s,count))
        }
      }
      
      rm(full_array)
      k<-s+1
      save(active_vox, count, count_bg, full_mat, bg_mat, meta_mat, k, file= paste(outdir,'/time_series/temp_data.rdata',sep=''))
  
  }
  
  print(paste(count_bg,  all(!is.na(bg_mat[1,]))))
  print(head(bg_mat[1,]))
  bg_mat[2,]<-(bg_mat[2,]-(bg_mat[1,]^2))/count_bg
  bg_mat[1,]<-bg_mat[1,]/count_bg
  
  print(head(bg_mat[1,]))
  
  
  print("Saving within cluster time series data") 
  save(full_mat,file=paste(outdir,'/time_series/main_data.rdata',sep=''))
  save(meta_mat,file=paste(outdir,'/time_series/meta_data.rdata',sep=''))
  save(bg_mat,file=paste(outdir,'/time_series/bg_data.rdata',sep=''))
}else{
  load(file=paste(outdir,'/time_series/main_data.rdata',sep=''))
  load(file=paste(outdir,'/time_series/meta_data.rdata',sep=''))
  load(file=paste(outdir,'/time_series/bg_data.rdata',sep=''))
}


###########################################
#delta F
#clusters
cl<-length(table(meta_mat[,4]))
time<-dim(full_mat)[2]
means<-array(0, c(time,cl))
sds<-array(0, c(time,cl))
q1<-array(0, c(2, time,cl))


for(i in 1:cl){
  
  means[,i]<-colMeans(full_mat[meta_mat[,4]==i,])
  
  sds[,i]<-sqrt(colMeans((full_mat[meta_mat[,4]==i,])^2)-(means[,i])^2)
  
  q1[,,i]<-t(colQuantiles(full_mat[meta_mat[,4]==i,], probs=c(.25,.75)))
  
}

evens<-as.vector((rep(c(F,T), length.out=time)))
grand_mean<-mean(full_mat)

if(cl>1){
  CORR<-cor(means[evens==F,])
  save(CORR,file=paste(outdir,'/time_series/correlation_data.rdata',sep=''))
  pdf(paste(outdir,'/Correlation_matrix.pdf',sep=''),paper='a4r')
  image.plot(1:cl,1:cl,CORR)
  dev.off()
}

    
means<-cbind(means, (bg_mat[1,]))
  
means2<-melt(means[evens==F,])
names(means2)<-c("Index", "Cluster", "Value")
means2$F0<-means2$Value/grand_mean
freq<-5 #5hz

means2$Seconds<-means2$Index/freq
means2$Cluster<-factor(means2$Cluster)
levels(means2$Cluster)[length(levels(means2$Cluster))]<-"Background"
# Plots
plot<-ggplot(data=means2, aes(y=F0,x=Seconds, colour=Cluster))+geom_line()+theme_bw()



pdf(paste(outdir,'/Mean_time_series_by_cluster.pdf',sep=''),paper='a4r')
print(plot)
dev.off()

plot<-ggplot(data=means2, aes(y=F0,x=Seconds))+geom_line()+theme_bw()+facet_wrap(~Cluster)
pdf(paste(outdir,'/Mean_time_series_by_cluster_facet.pdf',sep=''),paper='a4r')
print(plot)
dev.off()

clust_n<-length(levels(means2$Cluster))

for(i in 1:clust_n){
  temp<-subset(means2, as.numeric(means2$Cluster)==i)
  plot<-ggplot(data=temp, aes(y=F0,x=Seconds))+geom_line()+theme_bw()
  pdf(paste(outdir,'/Mean_time_cluster', levels(means2$Cluster)[i] ,'.pdf',sep=''),paper='a4r')
  print(plot)
  dev.off()
}


for(i in 1:clust_n){
  # #whole brain
  temp<-subset(means2, as.numeric(means2$Cluster)==i)
  
  if(i !=clust_n){
    temp$min<-q1[1,evens==F,i]/grand_mean
    temp$max<-q1[2,evens==F,i]/grand_mean
    plot<-ggplot(data=temp, aes(y=F0,x=Seconds))+geom_line()+theme_bw()+geom_ribbon(aes(ymin=min, ymax=max), alpha=0.3)+ggtitle(paste0("Cluster ", levels(means2$Cluster)[i]), "\n with centre 50% empirical distribution")
    pdf(paste(outdir,'/Mean_time_cluster', levels(means2$Cluster)[i] ,'_with_interval.pdf',sep=''),paper='a4r')
    print(plot)
    dev.off()
  }else{
    temp$min<-temp$F0-.77*sqrt(bg_mat[2,evens==F])/grand_mean
    temp$max<-temp$F0+.77*sqrt(bg_mat[2,evens==F])/grand_mean
    plot<-ggplot(data=temp, aes(y=F0,x=Seconds))+geom_line()+theme_bw()+geom_ribbon(aes(ymin=min, ymax=max), alpha=0.3)+ggtitle(paste0("Cluster ", levels(means2$Cluster)[i]), "\n with centre 50% normal CI")
    pdf(paste(outdir,'/Mean_time_cluster', levels(means2$Cluster)[i] ,'_with_interval.pdf',sep=''),paper='a4r')
    print(plot)
    dev.off()
  }
  # 

}
###########################################
#correlation analysis

###########################################
print("06-PostProc.R complete.")
