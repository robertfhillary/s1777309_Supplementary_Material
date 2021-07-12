setwd("./Epigenetic Clocks WHO/Incidence")

## Preparing of Phenotype Files 

DNAm_w1 <- read.csv("./STRADL_w1_DNAmAge_combined_agemonths.csv")
DNAm_w3 <- read.csv("./STRADL_w3_DNAmAge_combined.csv")

demog_w1 <- read.csv("./stradl-samples-5101.csv")
demog_w3 <-  read.csv("./samplesheet.final.csv")

demog_w1$Sample_Sentrix_ID <- paste(demog_w1$Sentrix_ID, demog_w1$Sentrix_Position, sep="_")
demog <- rbind(demog_w1[,c("Sample_Name", "Sample_Sentrix_ID")], demog_w3[,c("Sample_Name", "Sample_Sentrix_ID")])
DNAm_w1_1 <- DNAm_w1[,c(58,60,61,2,57,62:74,105:111,96:103)]
DNAm_w3_1 <- DNAm_w3[,c(58,60,61,2,57,62:74,105:111,96:103)]

DNAm <- rbind(DNAm_w1_1, DNAm_w3_1)

DNAm1 <- merge(DNAm, demog, by.x="id", by.y="Sample_Sentrix_ID")

pck_yrs <- read.csv("updated_smoking_Jan2019/pack_years.csv")
body <- read.csv("clinical/body.csv")
simd <- read.csv("Clinical/SIMD.csv")
ea <- read.csv("PCQ/education.csv")
alc <- read.csv("PCQ/alcohol.csv")
alc$units_usual <- ifelse(alc$drink_status==1 & alc$usual==2, alc$units, NA)

pheno <- merge(pck_yrs, body, by.x="Sample_Name", by.y="id", all=T)
pheno <- merge(pheno, ea, by.x="Sample_Name", by.y="id", all=T)
pheno <- merge(pheno, simd, by.x="Sample_Name", by.y="id", all=T)
pheno <- merge(pheno, alc, by.x="Sample_Name", by.y="id", all=T)
DNAm2 <- merge(DNAm1, pheno, by.x = "Sample_Name", by.y = "Sample_Name")

index1 = read.csv("../Discovery_Phenotype_1.csv")
index2 = read.csv("../Replication_Phenotype_1.csv")

age_alive <- read.csv("../age_may19.csv")
age_alive = age_alive[age_alive$dead %in% 0,]
age_dead = read.csv("../age_at_death.csv")
names(age_alive)[7] <- "aged"
age = rbind(age_alive, age_dead)
DNAm2 = merge(DNAm2, age, by.x = "Sample_Name", by.y = "id")


## Read in Incidence Files 

AD <- read.csv("ALZ.csv")
COPD <- read.csv("COPD.csv")
Depression <- read.csv("Dep.csv")
Heart.Disease <- read.csv("Heart.csv")
Stroke <- read.csv("Stroke.csv")
Lung.Cancer <- read.csv("Lung.csv")
Breast.Cancer <- read.csv("Breast.csv")
Bowel.Cancer <- read.csv("Bowel.csv")
Diabetes <- read.csv("T2D.csv")
Pain <- read.csv("Back.csv")


clock <- c("AgeAccelGrim","AgeAccelPheno","EEAA","IEAA","DNAmTLAdjAge") 
k <- list("AgeAccelGrim"=0, "AgeAccelPheno"=7, "EEAA"=14, "IEAA"=21, "DNAmTLAdjAge"=28)

## BASIC MODEL 

## Run Incidence of Alzheimer's Disease Analysis - Restricted to Participants > 60 

index1 = read.csv("../Discovery_Phenotype_1.csv")
index2 = read.csv("../Replication_Phenotype_1.csv")
index1$AD <- 0 
index2$AD <- 0
index1$AD <- ifelse(index1$Alzheimers.Disease...Maternal %in% 1 | index1$Alzheimers.Disease...Paternal %in% 1,1,0)
index2$AD <- ifelse(index2$Alzheimers.Disease.Maternal %in% 1 | index2$Alzheimers.Disease.Paternal %in% 1,1,0)

mat_ad <- matrix(nrow=1,ncol=7*5)
output_AD<- as.data.frame(mat_ad)

ind1 = index1[index1[,"AD"] %in% "1", "Sample_Name"]
ind2 = index2[index2[,"AD"] %in% "1", "Sample_Name"]
list= c(ind1, ind2)

for(j in clock) { 
  
  DNAm2$AD<- 0
  DNAm2[which(DNAm2$Sample_Name %in% AD$id),"AD"] <- 1
  DNAm2_AD<-DNAm2[-which(DNAm2$Sample_Name %in% list),]
  DNAm2_AD <- DNAm2_AD[DNAm2_AD$Age >= 60,]
  
output_AD[1,1+k[[j]]] <- "Alzheimer's Disease"
mod = glm(DNAm2_AD[,"AD"] ~ scale(DNAm2_AD[,j]) + scale(DNAm2_AD$Age) + factor(DNAm2_AD$Female), family = "binomial")
output_AD[1, 2:4+k[[j]]] <- exp(cbind(coef(mod), confint(mod)))[2,1:3]
output_AD[1,5+k[[j]]] <- summary(mod)$coefficients[2,4]
output_AD[1,6+k[[j]]] <- length(which(DNAm2_AD$AD %in% 1))
output_AD[1,7+k[[j]]] <- length(which(DNAm2_AD$AD %in% 0))
} 
## Create List of Remaining Dataframes - not Cancer  

my.list = list(COPD,Depression,Heart.Disease,Stroke,Diabetes,Pain)
names = list("COPD","Depression","Heart.Disease","Stroke","Diabetes","Pain")
names(my.list) <- names 

mat <- matrix(nrow=length(my.list),ncol=7*5)
output<- as.data.frame(mat)

l=lapply(my.list, "[", 1)

index1$Pain <- 0 
index2$Pain <- 0
index1$Pain <- ifelse(index1$Back.Pain %in% 1 | index1$Neck.Pain %in% 1,1,0)
index2$Pain <- ifelse(index2$Back.Pain %in% 1 | index2$Neck.Pain %in% 1,1,0)


## Loop of Logistic Regression Models - Cross-Sectional Associations


for(j in clock) { 
for(i in 1:length(l)){ 
  
  tmp <- l[[i]]
  
  
  ind1 = index1[index1[,names[[i]]] %in% "1", "Sample_Name"]
  ind2 = index2[index2[,names[[i]]] %in% "1", "Sample_Name"]
  list= c(ind1, ind2)
  DNAm2[,names[[i]]] <- 0
  DNAm2[which(DNAm2$Sample_Name %in% tmp$id),names[[i]]] <- 1
  tmp2=DNAm2[-which(DNAm2$Sample_Name %in% list),]
  
  output[i,1+k[[j]]] <- names[[i]]
  mod = glm(tmp2[,names[[i]]] ~ scale(tmp2[,j]) + scale(tmp2$Age) + factor(tmp2$Female), family = "binomial")
  output[i, 2:4+k[[j]]] <- exp(cbind(coef(mod), confint(mod)))[2,1:3]
  output[i,5+k[[j]]] <- summary(mod)$coefficients[2,4]
  output[i,6+k[[j]]] <- length(which(tmp2[,names[[i]]] %in% 1))
  output[i,7+k[[j]]] <- length(which(tmp2[,names[[i]]] %in% 0))

  } 
} 

## Create List of Remaining Dataframes - Cancer  

mat_cancer <- matrix(nrow=length(my.list.cancer),ncol=7*5)
output_cancer<- as.data.frame(mat_cancer)

my.list.cancer = list(Bowel.Cancer,Bowel.Cancer,Lung.Cancer)
names = list("Breast.Cancer","Bowel.Cancer","Lung.Cancer") 
names(my.list.cancer) <- names 

l.cancer=lapply(my.list.cancer, "[", 1)
smr=lapply(my.list.cancer, "[", c(1,4))
smr1=lapply(smr,subset, smr==1)


## Loop of Logistic Regression Models - Cross-Sectional Associations

for(j in clock){ 
for(i in 1:length(l.cancer)){ 
  
  tmp <- l.cancer[[i]]
  smr_tmp <-  smr1[[i]]

  ind1 = index1[index1[,names[[i]]] %in% "1", "Sample_Name"]
  ind2 = index2[index2[,names[[i]]] %in% "1", "Sample_Name"]
  list= c(ind1, ind2)
  DNAm2[,names[[i]]] <- 0
  DNAm2[which(DNAm2$Sample_Name %in% tmp$id),names[[i]]] <- 1
  tmp2=DNAm2[-which(DNAm2$Sample_Name %in% list),]

## Exclude Indiviudals Reported on Cancer Inpatient Day Visit List - might not have had cancer 
   tmp2=tmp2[-which(tmp2$Sample_Name %in% smr_tmp$id),]
  
  
  output_cancer[i,1+k[[j]]] <- names[[i]]
  mod = glm(tmp2[,names[[i]]] ~ scale(tmp2[,j]) + scale(tmp2$Age) + factor(tmp2$Female), family = "binomial")
  output_cancer[i, 2:4+k[[j]]] <- exp(cbind(coef(mod), confint(mod)))[2,1:3]
  output_cancer[i,5+k[[j]]] <- summary(mod)$coefficients[2,4]
  output_cancer[i,6+k[[j]]] <- length(which(tmp2[,names[[i]]] %in% 1))
  output_cancer[i,7+k[[j]]] <- length(which(tmp2[,names[[i]]] %in% 0))
} 

} 
## Combine Outputs into One Dataframe 

comb = rbind(output_AD, output)
comb = rbind(comb, output_cancer)
comb <- comb[,c(1,6,7,2:5,(c(1,6,7,2:5)+7),(c(1,6,7,2:5)+14),(c(1,6,7,2:5)+21),(c(1,6,7,2:5)+28))]
names(comb) <- rep(c("Trait", "no. of cases", "no. of controls", "Hazard Ratio", "LCI", "UCI", "P"),5)

## FULLY-ADJUSTED MODEL 

## Run Incidence of Alzheimer's Disease Analysis - Restricted to Participants > 60 

index1 = read.csv("../Discovery_Phenotype_1.csv")
index2 = read.csv("../Replication_Phenotype_1.csv")
index1$AD <- 0 
index2$AD <- 0
index1$AD <- ifelse(index1$Alzheimers.Disease...Maternal %in% 1 | index1$Alzheimers.Disease...Paternal %in% 1,1,0)
index2$AD <- ifelse(index2$Alzheimers.Disease.Maternal %in% 1 | index2$Alzheimers.Disease.Paternal %in% 1,1,0)

mat_ad <- matrix(nrow=1,ncol=7*5)
output_AD<- as.data.frame(mat_ad)

ind1 = index1[index1[,"AD"] %in% "1", "Sample_Name"]
ind2 = index2[index2[,"AD"] %in% "1", "Sample_Name"]
list= c(ind1, ind2)

DNAm2$AD<- 0
DNAm2[which(DNAm2$Sample_Name %in% AD$id),"AD"] <- 1
DNAm2_AD<-DNAm2[-which(DNAm2$Sample_Name %in% list),]
DNAm2_AD <- DNAm2_AD[DNAm2_AD$Age >= 60,]
for(j in clock){ 
output_AD[1,1+k[[j]]] <- "Alzheimer's Disease"
mod = glm(DNAm2_AD[,"AD"] ~ scale(DNAm2_AD[,j]) + scale(DNAm2_AD$Age) + factor(DNAm2_AD$Female) + scale(DNAm2_AD$pack_years) + scale(DNAm2_AD$units_usual) + scale(DNAm2_AD$rank) + scale(DNAm2_AD$years) + scale(DNAm2_AD$bmi), family = "binomial")
output_AD[1, 2:4+k[[j]]] <- exp(cbind(coef(mod), confint(mod)))[2,1:3]
output_AD[1,5+k[[j]]] <- summary(mod)$coefficients[2,4]
output_AD[1,6+k[[j]]] <- length(which(DNAm2_AD$AD %in% 1))
output_AD[1,7+k[[j]]] <- length(which(DNAm2_AD$AD %in% 0))
} 
## Create List of Remaining Dataframes - not Cancer  

my.list = list(COPD,Depression,Heart.Disease,Stroke,Diabetes,Pain)
names = list("COPD","Depression","Heart.Disease","Stroke","Diabetes","Pain")
names(my.list) <- names 

mat <- matrix(nrow=length(my.list),ncol=7*5)
output<- as.data.frame(mat)

l=lapply(my.list, "[", 1)

index1$Pain <- 0 
index2$Pain <- 0
index1$Pain <- ifelse(index1$Back.Pain %in% 1 | index1$Neck.Pain %in% 1,1,0)
index2$Pain <- ifelse(index2$Back.Pain %in% 1 | index2$Neck.Pain %in% 1,1,0)

## Loop of Logistic Regression Models - Cross-Sectional Associations
for(j in clock){
for(i in 1:length(l)){ 
  
  tmp <- l[[i]]
  
  
  ind1 = index1[index1[,names[[i]]] %in% "1", "Sample_Name"]
  ind2 = index2[index2[,names[[i]]] %in% "1", "Sample_Name"]
  list= c(ind1, ind2)
  DNAm2[,names[[i]]] <- 0
  DNAm2[which(DNAm2$Sample_Name %in% tmp$id),names[[i]]] <- 1
  tmp2=DNAm2[-which(DNAm2$Sample_Name %in% list),]
  
  output[i,1+k[[j]]] <- names[[i]]
  mod = glm(tmp2[,names[[i]]] ~ scale(tmp2[,j]) + scale(tmp2$Age) + factor(tmp2$Female) + scale(tmp2$pack_years) + scale(tmp2$units_usual) + scale(tmp2$years) + scale(tmp2$rank) + scale(tmp2$bmi), family = "binomial")
  output[i, 2:4+k[[j]]] <- exp(cbind(coef(mod), confint(mod)))[2,1:3]
  output[i,5+k[[j]]] <- summary(mod)$coefficients[2,4]
  output[i,6+k[[j]]] <- length(which(tmp2[,names[[i]]] %in% 1))
  output[i,7+k[[j]]] <- length(which(tmp2[,names[[i]]] %in% 0))
  
} 
} 
## Create List of Remaining Dataframes - Cancer  

mat_cancer <- matrix(nrow=length(my.list.cancer),ncol=7*5)
output_cancer<- as.data.frame(mat_cancer)

my.list.cancer = list(Bowel.Cancer,Bowel.Cancer,Lung.Cancer)
names = list("Breast.Cancer","Bowel.Cancer","Lung.Cancer") 
names(my.list.cancer) <- names 

l.cancer=lapply(my.list.cancer, "[", 1)
smr=lapply(my.list.cancer, "[", c(1,4))
smr1=lapply(smr,subset, smr==1)


## Loop of Logistic Regression Models - Cross-Sectional Associations
for(j in clock){ 
for(i in 1:length(l.cancer)){ 
  
  tmp <- l.cancer[[i]]
  smr_tmp <-  smr1[[i]]
  
  ind1 = index1[index1[,names[[i]]] %in% "1", "Sample_Name"]
  ind2 = index2[index2[,names[[i]]] %in% "1", "Sample_Name"]
  list= c(ind1, ind2)
  DNAm2[,names[[i]]] <- 0
  DNAm2[which(DNAm2$Sample_Name %in% tmp$id),names[[i]]] <- 1
  tmp2=DNAm2[-which(DNAm2$Sample_Name %in% list),]
 
## Exclude Indiviudals Reported on Cancer Inpatient Day Visit List - might not have had cancer 
  
   tmp2=tmp2[-which(tmp2$Sample_Name %in% smr_tmp$id),]
  
  
  output_cancer[i,1+k[[j]]] <- names[[i]]
  mod = glm(tmp2[,names[[i]]] ~ scale(tmp2[,j]) + scale(tmp2$Age) + factor(tmp2$Female) + scale(tmp2$pack_years) + scale(tmp2$units_usual) + scale(tmp2$years) + scale(tmp2$rank) + scale(tmp2$bmi), family = "binomial")
  output_cancer[i,2:4+k[[j]]] <- exp(cbind(coef(mod), confint(mod)))[2,1:3]
  output_cancer[i,5+k[[j]]] <- summary(mod)$coefficients[2,4]
  output_cancer[i,6+k[[j]]] <- length(which(tmp2[,names[[i]]] %in% 1))
  output_cancer[i,7+k[[j]]] <- length(which(tmp2[,names[[i]]] %in% 0))
} 
} 
## Combine Outputs into One Dataframe 

comb_full = rbind(output_AD, output)
comb_full = rbind(comb_full, output_cancer)
comb_full <- comb_full[,c(1,6,7,2:5,(c(1,6,7,2:5)+7),(c(1,6,7,2:5)+14),(c(1,6,7,2:5)+21),(c(1,6,7,2:5)+28))]
names(comb_full) <- rep(c("Trait", "no. of cases", "no. of controls", "Hazard Ratio", "LCI", "UCI", "P"),5)

## Final Files 

View(comb)
View(comb_full)
