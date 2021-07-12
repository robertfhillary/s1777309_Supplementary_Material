##### PREPARATION OF PHENOTYPE FILE #####

phenos = read.csv("phenos.csv") 

phenos$g[phenos$g %in% NA] <- mean(phenos$g, na.rm = T)

write.table(x = t(as.matrix(as.numeric(scale(phenos$g)))),file = "pheno.csvphen" ,quote = F, sep = ",", row.names = F, col.names = F)



##### PREPARATION OF METHYLATION (or Single OMICS) FILE #####

## read in GS methylation file

meth = readRDS("mvals.rds") 

## Read in phenotypes file 

phenos = read.csv("phenos.csv") 

## Match order of ids in both files 

ids = phenos$meth_id 
meth1 = meth[,match(ids,colnames(meth))] 
table(colnames(meth1) == phenos$meth_id) 

## Mean impute NAs - do this for every CpG site (transpose df so CpGs are columns) 

tmeth1 = t(meth1) 
tmeth1 = as.data.frame(tmeth1) 
library(imputeTS) 
tmeth1 = na_mean(tmeth1) 

## Scale every row - should be CpGs as rows (transpose df) 

tmeth1 = t(tmeth1) 
tmeth1 = as.data.frame(tmeth1) 

tmeth1 = apply(tmeth1, 1, scale)

## Final file should be 4450 individuals as columns and ~800k probes as rows 

write.table(x = t(tmeth1), "Single_Omics.csv", sep = ",", row.names = F, col.names = F, quote = F) 

## Keep CpG Identifiers as row.names are deleted from file for BayesR+ 

names1 <- colnames(tmeth1) 
saveRDS(names1, "CpG_Names_BayesR.rds") 

###### To scale up to two omics type, prior to processing, rbind/cbind file of CpGs and SNPs so that for every individual their value for each CpG and SNP is shown. Apply same processing principles. 

## PREPARATION OF COVARIATES FILE 

cov <- read.csv("Some_Covariates_File.csv")

## Scale all of the columns of interest and make sure no NA are present
## Final structure should be columns as covariates and rows as individuals matched in order to that of other files  

write.table(x = as.matrix(cov),file = "Some_Scaled_Covariates_File.csv" ,quote = F, sep = ",", row.names = F, col.names = F)



## PREPARATION OF GROUPS FILE

## This is a .txt file 
## First column shows sites - in my example, CpGs and SNPs
## Second column shows group - either 0 for CpGs or 1 for SNPs
## Order of CpGs and SNPs matches the order of the sites in the row.names of Omics.csv file  
