## Epigenome Wide Analyses of Protein Levels 


## Linear Regression Method 


## Installing Requisite Packages

if(!require(limma)){
  install.packages("limma")
}

if(!require(dplyr)){
  install.packages("dplyr")
}


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")

library(limma)
library(dplyr)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)



pData <- read.csv("pData.csv")

meth <- load("methylation.RObject")

a = which(colnames(dat) %in% pData$Basename)
dat1 = dat[,a]
ids = pData$Basename
dat1 = dat1[,match(ids, colnames(dat1))]

rownames(pData) <- pData$Basename


all.equal(rownames(pData), colnames(dat1))


saveRDS(dat1, "Methylation_Inflam.rds")
saveRDS(pData, "Phenotypic_Inflam.rds")


## Models  
IL8 <- model.matrix(~IL8 + agedays_w1 + as.factor(sex.x) + Mono + Gran + CD8T + NK + Bcell + date + plate + set + pos + array, data=pData)

VEGFA <- model.matrix(~VEGFA + agedays_w1 + as.factor(sex.x) + Mono + Gran + CD8T + NK + Bcell + date + plate + set + pos + array, data=pData)

MCP3 <- model.matrix(~MCP.3 + agedays_w1 + as.factor(sex.x) + Mono + Gran + CD8T + NK + Bcell + date + plate + set + pos + array, data=pData)



fits <- list()
models <- c("IL8", "VEGFA", "MCP3")

for(i in models) {
print(i)
timestamp()
fits[[i]] <- lmFit(dat1[,rownames(get(i))], get(i))
}

saveRDS(fits, file="filename.rds")

fits2 <- list()
for(i in models) {
print(i)
timestamp()
fits2[[i]] <- eBayes(fits[[i]])
}
saveRDS(fits2, file="filename.rds")

anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
genes <- data.frame(ID=anno$Name, geneSymbol=anno$UCSC_RefGene_Name, CHR=anno$chr, MAPINFO=anno$pos, FEATURE=anno$UCSC_RefGene_Group, CpGISLAND=anno$Relation_to_Island)
rownames(genes) <- genes$ID

tt <- list()
for(i in models){
tt[[i]] <- topTable(fits2[[i]], coef=2, adjust.method="none", number=nrow(fits2[[i]]),
                  genelist=genes[rownames(fits2[[i]]),c("geneSymbol", "CHR", "MAPINFO", "FEATURE", "CpGISLAND")])
 }
saveRDS(tt, file="filename.rds")


## Mixed Model Method 

 ## Preparation for OSCA run in Terminal 
 
 ## (run in R)
 osca_dat <- data.frame(IID=colnames(dat), FID=colnames(dat))
 osca_dat <- cbind(osca_dat, t(dat))
 write.table(osca_dat, file="myprofile.txt", row.names=F, sep=' ')


## #Youll have to manually make the .opi annotation file:
## Manual Annotation 
 opi <- anno[colnames(osca_dat)[3:ncol(osca_dat)],c("chr","Name", "pos","UCSC_RefGene_Name", "strand")]
 opi$chr <- gsub("chr", "", opi$chr)
 opi$chr <- gsub("X", "23", opi$chr)
 opi$chr <- gsub("Y", "24", opi$chr)
 opi$chr <- as.numeric(opi$chr)
 opi$UCSC_RefGene_Name	 <- as.factor(opi$UCSC_RefGene_Name	)
 opi$strand <- as.factor(opi$strand)
 opi[which(opi$UCSC_RefGene_Name	==""), "UCSC_RefGene_Name"] <- NA
 write.table(opi, file="output_filenam.opi", 
 	             col.names=F, 
 	             row.names=F, 
 	             quote=F, sep='\t')


## Make phenotype file
 ## 3 columns: IID, FID, Phenotype
 pheno <- data.frame(FID = phenos$Basename,
 	                IID = phenos$Basename,
 	                simd = scale(phenos$ pheno))

 ## Make quantitative covariates file (e.g. age, pack years, cellcounts etc) 
 ## First two columns are IID and FID as before

 quant_cov <- data.frame(FID = phenos$Barcode, 
 	                    IID = phenos$Barcode,
 	                    CD8T = phenos$CD8T, 
 	                    CD4T = phenos$CD4T,
 	                    Mono = phenos$Mono,
 	                    Bcell = phenos$Bcell,
 	                    Gran = phenos$Gran,
 	                    Age = phenos$Age,
 	                    pack_years = phenos$pack_years)
 write.table(quant_cov, file="quant.cov", row.names=F, sep=' ')


# Make factor covariates file (e.g. smoking status, sex, batch etc) 
# First two columns are IID and FID as before
 fact_cov <- data.frame(FID =phenos$Barcode,
 	                   IID = phenos$Barcode,
 	                   ever_smoke = phenos$ever_smoke,
 	                   batch = phenos$plate_processing_batch)
 write.table(fact_cov, file="factors.cov", row.names=F, sep=' ')


## Run OSCA in terminal ##
# Create binary methylation data --methylation-beta if beta values, --methylation-m if M-values

system("osca_Linux --efile myprofile.txt --methylation-beta --make-bod --out output_filename")

osca_Linux \
	   --moment \
 	   --befile ../output_filename \
 	   --pheno ../simd.phen \
 	   --qcovar ../quant.cov \
 	   --covar ../factors.cov \
 	   --fast-linear \
 	   --out osca_simd_ewas \
 	   --methylation-m
 	   
 	   
 	   
## Automation of Genome Wide Association Study 

## GWAS Preparation Code 

 ## Here columns 2:92 represent levels of 91 different proteins - run in R 
 
R 
for(i in 2:92){
out = prot[,c(1,1,103,103,101,i)]
var = names(prot)[i]

var = gsub("\\.", "_", var)

dir.create(paste("Inputs_For_GWAS/", var, sep=""))
write.table(out, file=paste("Inputs_For_GWAS/", var, "/", var, "_res.ped", sep=""), quote=F, col.names=F, row.names=F)
dat = t(c("T", var))
write.table(dat, file=paste("Inputs_For_GWAS/", var, "/", var, "_res.dat", sep=""), quote=F, col.names=F, row.names=F, sep="\t")
}


## Automation of GWAS - run in terminal: mach2qtl for imputation 

cd /Inputs_For_GWAS/

FILES=(*) 
## First 15 proteins
for i in ${FILES[@]:0:15}
do 
for x in {1..22} 
do 
data=/Inputs_For_GWAS/$i/$i_res.dat
pedi=/Inputs_For_GWAS/$i/$i_res.ped
info=/Data/CHR${x}_results.info; dose=/Data/CHR${x}_results.dose 
mach2qtl-ext -d $data -p $pedi -i $info --dosefile $dose --samplesize > /Inputs_For_GWAS/$i/$i_output_chr${x}.out
done
done




