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
indir_new <- args[1] # $input = "04-Splines"; 
outdir_new <- args[2]
FROM <- as.numeric(args[3])
TO <- as.numeric(args[4])
mask_fn <- args[5]
indir_new_2 <- args[6] # $input2 = "05-Clusters";

load(file=paste(indir_new,'/settings.rdata',sep=''))

outdir<-outdir_new
indir<-indir_new
indir2<-indir_new_2

max_clust<-30 #speed up by making this look at smaller
cut_p<-0.04

print(c('1indir', indir))
print(c('2outdir', outdir))
print(c('3FROM', FROM))
print(c('4TO', TO))
print(c('5mask', mask_fn))
print(c('6indir2', indir_new_2))

print(outdir)
if(!dir.exists(outdir)){
  dir.create(outdir, recursive=T)
}

#create subfolders
dir.create(paste0(outdir,'/clusters/'))

#reload misc
load(file=paste(indir,'/ssDM.rdata',sep=''))
load(file=paste(indir,'/MASK_hdr.rdata',sep=''))
load(file=paste(indir,'/D_Mask.rdata',sep=''))
load(file=paste(indir,'/MASK_active.rdata',sep=''))

#reload params
load(file = paste0(indir2,"/scalingparams.rdata"))
load(file=paste0(indir2,'/centers.rdata'))
load(file=paste(indir2,'/clustering.rdata',sep=''))#????

load(file=paste(indir,'/D_Mask.rdata',sep=''))
load(file=paste(indir2,'/MASK_active_adjusted.rdata',sep=''))
load(file=paste(indir2,'/ssDM_active.rdata', sep=''))
load(file=paste(indir,'/clusters/image_hold.rdata', sep=''))













# # ### Obtain the cluster mean time series
# # # Reload Graphing Parameters (These values are specific to F03 and F04)
# 
# file_number <- (file_number.old -1)*Z_SIZE + Z_SIZE
# Basis <- create.bspline.basis(c(0,(max_file-1)*Z_SIZE+Z_SIZE),nbasis=Basis_number)
# BS <- eval.basis(file_number,Basis)
# 
# # Compute the mean time series
# load(file=paste0(indir2,'/centers.rdata'))
# load(file = paste0(indir2,"/scalingparams.rdata"))
# 
# pred_range <- eval.basis(seq(from=1,to=((max_file-1)*Z_SIZE+Z_SIZE),length.out=1000),Basis)
# 
# # Each row is a mean time series over the pred_range values
# PRED <- matrix(NA,dim(centers)[1],dim(pred_range)[1])
# for (j in 1:dim(centers)[1]) {
#   PRED[j,] <- apply(pred_range,1,function(x) {sum(x*centers[j,])})
# }
# 
# # The time average values of each series
# MEAN_PRED <- rowMeans(PRED)
# 
# ### Make a set of functions for evaluating splines and convolutions of splines
# Spline_function <- function(x,cc) {sum(eval.basis(x,Basis)*centers[cc,])}
# S1 <- function(x) Spline_function(x,1)
# S2 <- function(x) Spline_function(x,2)
# S_Prod <- function(x) S1(x)*S2(x)
# # INTEGRAL <- integrate(Vectorize(S_Prod),1,(max_file-1)*Z_SIZE+Z_SIZE)$value
# 
# ### Compute the Covariance of each mean function
# COVAR <- c()
# for (cc in 1:dim(centers)[1]) {
#   S1 <- function(x) Spline_function(x,cc)
#   S2 <- function(x) Spline_function(x,cc)
#   S_Prod <- function(x) S1(x)*S2(x)
#   INTEGRAL <- quadv(S_Prod,1,(max_file-1)*Z_SIZE+Z_SIZE)$Q
#   COVAR[cc] <- INTEGRAL
# }
# 
# comp<-dim(centers)[1]
# ### Compute the Correlation between pairs of mean functions
# CORR <- matrix(NA,comp,comp)
# for (c1 in 1:comp) {
#   for (c2 in c1:comp) {
#     S1 <- function(x) Spline_function(x,c1)
#     S2 <- function(x) Spline_function(x,c2)
#     S_Prod <- function(x) S1(x)*S2(x)
#     QUAD <- quadv(S_Prod,1,(max_file-1)*Z_SIZE+Z_SIZE)
#     INTEGRAL <- QUAD$Q
#     INT_OLD <- INTEGRAL
#     PREC <- QUAD$estim.prec
#     CORR[c1,c2] <- INTEGRAL/sqrt(COVAR[c1]*COVAR[c2])
#     #what does this do
#     counter_corr = 0
#     while( xor(CORR[c1,c2] > 1, CORR[c1,c2] < -1) & counter_corr < 100 ) {
#       if (CORR[c1,c2] > 1) {INTEGRAL <- INTEGRAL - PREC*INT_OLD}
#       if (CORR[c1,c2] < -1) {INTEGRAL <- INTEGRAL + PREC*INT_OLD}
#       CORR[c1,c2] <- INTEGRAL/sqrt(COVAR[c1]*COVAR[c2])
#       counter_corr=counter_corr+1
#     }
#   }
# }
# 
# if(comp>2){
#   # Hierachical Clustering ------------------------------
#   # # Make a distance metric
#   DIST <- as.dist(1-t(CORR))
#   HCLUST <- hclust(DIST,method='average')
# 
# }

print("06-PostProc.R complete.")
