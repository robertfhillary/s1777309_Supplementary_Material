## Normal-exponential out-of-band (noob) normalisation was performed on LBC1936 methylation data 
## in order to avoid missing CpG values for the calculation of DNAm GrimAge as recommended by Horvath 

## Installing Requisite Packages 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("minfi")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("IlluminaHumanMethylation450kmanifest")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")


library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)


cat("Packages installed and loaded...")

## Read and load green and red channel methylation data 

cat("About to read IDATs \n")
RGset=minfi::read.metharray.exp("IDATs/", verbose=T)
save(RGset, file="RGset.Robject")

cat("About to load RGset file \n")
load(file="RGset.Robject")
cat("RGset file loaded \n")

pvals=minfi::detectionP(RGset, type = "m+u")
cat("Pval Detection Complete \n")
save(pvals, file="DetectionPvalues.Robject")
rm(pvals)

## Noob preprocessing

cat("About to Noob preprocess \n")
MSet.noob=minfi::preprocessNoob(RGset, dyeCorr = TRUE, dyeMethod="single")
cat("MSet.noob calculated \n")
save(MSet.noob, file="MSet_preprocessNoob.Robject")
cat("MSet.noob saved \n")

## Obtain noob-normalised beta values and save them

betas=minfi::getBeta(MSet.noob, betaThreshold = 0.001)
cat("noob betas calculated \n")
save(betas, file="Beta_preprocessNoob_BetaThreshold0pt001.Robject")
cat("noob betas saved \n")
rm(MSet.noob)
rm(betas)







