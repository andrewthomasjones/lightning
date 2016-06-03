### Load necessary libraries (all CRAN libraries can be acquired using install.packages(<pkgname>))

library(AnalyzeFMRI)
library(abind)


args <- commandArgs(trailingOnly = TRUE)
indir <- args[1]
outfile <- args[2]

print(c('INDIR: ', indir))
print(c('OUTFILE: ', outfile))

### Read in 75th Z slice
for (ss in 75) {

  ### Declare what files to load in
  ### Remove time points that do not look stationary
  file_list <- list.files(path=indir)
  file_list <- file_list[grep('.nii', file_list)]

  ### Declare number of time slices left and initialize list
  N <- length(file_list)
  full_array <- list()

  for (ii in 1:N) {
    NIFTI <- f.read.nifti.slice(paste(indir, file_list[ii], sep='/'), ss, 1)
    full_array[[ii]] <- NIFTI
    print(paste('reading slice', ss, ii, ii/N, sep=" - "))
  }

  ### Combine all of the time slices
  full_array <- abind(full_array, along=3)

  ### Store data into matrix instead of array
  count <- 0
  full_mat <- matrix(NA,prod(170*70),N)
  for (ii in 1:170) {
    for (jj in 1:70) {
      count <- count + 1
      full_mat[count,] <- full_array[ii,jj,]
      print(paste('storing', ss, count, sep=" - "))
    }
  }
}

# Plot stuff
pdf(outfile, width=50, height=10)
par(mfrow=c(2,5))
for (ii in 1:10) {
  plot(full_mat[sample(10000,1),],type='l')
}
dev.off()
