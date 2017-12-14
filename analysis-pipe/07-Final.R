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

a<-"/five"
b<-"/eleven"
c<-"/reduce"

who<-strsplit(indir, "\\/")[[1]][5]
#need to genralise this

if(who=="20150628_5dpf_H2BS_CON_LR_F01_0"){
   fold<-c
   clust_ID<-15  
}

if(who=="20150629_6dpf_H2BS_PTZ_LR_F1_1"){
    fold<-a
    clust_ID<-3
}




outdir2 <-paste0(outdir,fold)

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


load(file=paste(outdir2,'/clusters/image_hold_merge.rdata',sep=''))
load(file=paste(outdir2,'/clusters/tree.rdata',sep=''))

if(file.exists(file= paste(outdir2,'/time_series/temp_data2.rdata',sep=''))){  
  load(file= paste(outdir2,'/time_series/temp_data2.rdata',sep=''))
}else{
    k=1
    clust_members <-sum(image_hold2==clust_ID)
    full_mat<-matrix(NA, c(clust_members, N))
    meta_mat<-matrix(NA, c(clust_members ,4))
    count<-0
}
    if(k<Z_SIZE){
      for (s in k:Z_SIZE){

        sink(file=paste(outdir2,'/time_series/prog2.txt',sep=''))
        cat(s)
        cat("\n")
        sink()
        print(paste("Loading slice",s,"of", Z_SIZE))
        file_number <- (file_number.old-1)*Z_SIZE +s+Z_START
        print(mem_used())
        full_array <- list()

          ### Read in Z slices
          for (j in 1:N) {
            file_name <- substring(file_list[j],1,12)
            nii_name <- paste(indir3, '/', file_name,'.nii',sep='')
            print(c('reading',nii_name, s+Z_START,j,round(j/N),1))
            NIFTI <- f.read.nifti.slice(nii_name,s+Z_START,1)
            full_array[[j]] <- NIFTI
          }

          ### Combine all of the time slices
          full_array <- abind(full_array,along=3)

          s<-s+Z_START
          temp1<-image_hold2[,,s]
          
          for (j in 1:Y_SIZE) {
            for (i in 1:X_SIZE) {
      
              if(temp1[i,j]==clust_ID){
                  temp2<-full_array[i+X_START,j+Y_START,]
                  temp2[is.nan(temp2)]<-0
                  full_mat[count,] <- temp2
                  meta_mat[count,] <- c(i,j,s,image_hold2[i,j,s])
                  count<-count+1
              }
              #print(c('storing',s,count))
            }
          }

          rm(full_array)
          k<-s-Z_START+1
          
          save(full_mat, meta_mat, k, count,  file= paste(outdir2,'/time_series/temp_data2.rdata',sep=''))

      }
    }

  save(full_mat, meta_mat, k, file= paste(outdir2,'/time_series/for_correlation.rdata',sep=''))
