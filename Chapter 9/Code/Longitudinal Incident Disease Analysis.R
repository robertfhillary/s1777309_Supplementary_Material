setwd("./Epigenetic Clocks WHO//Incidence")

## Installing Requisite Packages 

if(!require(survival)){
  install.packages("survival")
}

library(survival)

## Preparing of Phenotype Files

d1 =readRDS("PheWAS_phenotypes_GS_10k.rds")
pain =read.csv("chronic_painv5.csv")

age_alive <- read.csv("../age_may19.csv")
age_alive = age_alive[age_alive$dead %in% 0,]
age_dead = read.csv("../age_at_death.csv")
names(age_alive)[7] <- "aged"
age = rbind(age_alive, age_dead)

d1 <- merge(d1,age,by.x = "Sample_Name", by.y="id")
d1$Pain <- 0
d1[d1$Sample_Name %in% pain[(pain$back_pain %in% 1 | pain$neck_pain %in% 1),"ID"],"Pain"] <- 1

d1$AD <- 0 
d1[d1$Sample_Name %in% d1[(d1$alzheimers_M %in% 1 | d1$alzheimers_F %in% 1),"Sample_Name"],"AD"] <- 1

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



clock <- c("AgeAccelGrim", "Dunedin_PoAm", "AgeAccelPheno","EEAA","IEAA","DNAmTLAdjAge") 
k <- list("AgeAccelGrim"=0, "Dunedin_PoAm" = 7, "AgeAccelPheno"=14, "EEAA"=21, "IEAA"=28, "DNAmTLAdjAge"=35)


## HAZARD MODELS - PREDICTING TIME-TO-ONSET OF DISEASES FROM STUDY BASELINE (2006)

## BASIC MODEL


## Run Incidence of Alzheimer's Disease Analysis - Restricted to Participants > 60 

d1_AD <- d1[d1$Age >= 60,]

mat_hazard_ad <- matrix(nrow=1,ncol=7*6)
output_hazard_AD<- as.data.frame(mat_hazard_ad)
for(j in clock){ 
  dat1= d1_AD[-which(d1_AD[,"AD"] %in% 1),]
  tmp1 = AD[which(AD$id %in% dat1$Sample_Name),]
  
  ## Obtain Age of Onset 
  affected = dat1[which(dat1$Sample_Name %in% tmp1$id),] 
  age_onset = AD[,c("first", "id")]
  affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
  affected$Event = 1
  affected$yoe = substring(affected$first, 1, 4)
  affected$moe = substring(affected$first, 5,6)
  affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
  affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
  affected$age_event = affected$age_event1 + affected$month_event1
  affected$first = NULL
  affected$yoe = NULL 
  affected$moe = NULL
  affected$month_event1 = NULL 
  affected$age_event1 = NULL
  
  healthy = dat1[-which(dat1$Sample_Name %in% AD$id),]
  healthy$Event = 0
  healthy$age_event = 0 
  affected$id.y <- NULL
  healthy$id <- NULL
  names(affected)[names(affected)=="id"] <- "Sample_Name"
  cox = rbind(affected, healthy)
  
  cox$age_death = 0
  cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
  cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
  cox$tte = cox$age_at_event - cox$Age
  cox$tte = as.numeric(cox$tte)
  cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
  cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
  cox$Event = as.numeric(cox$Event)
  cox$tte<-as.numeric(cox$tte)
  
  mod = coxph(Surv(event = cox$Event, cox$tte) ~ scale(cox[,j]) + cox$Age + factor(cox$Female))
  output_hazard_AD[1,1+k[[j]]] <- as.character("Alzheimer's Disease")
  output_hazard_AD[1, 2:4+k[[j]]]<-exp(cbind(coef(mod), confint(mod)))[1,1:3]
  output_hazard_AD[1,5+k[[j]]] <- summary(mod)$coefficients[1,5]
  output_hazard_AD[1,6+k[[j]]] <- length(which(cox$Event %in% 1))
  output_hazard_AD[1,7+k[[j]]] <- length(which(cox$Event %in% 0))
} 

## Create List of Remaining Dataframes - not Cancer  

my.list = list(COPD,Depression,Heart.Disease,Stroke,Diabetes,Pain)
names = list("COPD","Depression","Heart.Disease","Stroke","Diabetes","Pain")
names(my.list) <- names 

l=lapply(my.list, "[", c(1:5))

names(d1)[names(d1) == "COPD_Y"] <- "COPD"
names(d1)[names(d1) == "depression_Y"] <- "Depression"
names(d1)[names(d1) == "heart_disease_Y"] <- "Heart.Disease"
names(d1)[names(d1) == "stroke_Y"] <- "Stroke"
names(d1)[names(d1) == "diabetes_Y"] <- "Diabetes"

mat_hazard <- matrix(nrow=length(my.list),ncol=7*6)
output_hazard <- as.data.frame(mat_hazard)

## Loop of Logistic Regression Models - Cross-Sectional Associations
for(j in clock){
  for(i in 1:length(l)){ 
    
    tmp <- l[[i]]
    
    ## Exclude Indiviudals who Reported Disease at Study Baseline  
    dat1= d1[-which(d1[,names[[i]]] %in% 1),]
    tmp1 = tmp[which(tmp$id %in% dat1$Sample_Name),]
    
    ## Obtain Age of Onset 
    affected = dat1[which(dat1$Sample_Name %in% tmp1$id),] 
    age_onset = tmp[,c("first", "id")]
    affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
    affected$Event = 1
    affected$yoe = substring(affected$first, 1, 4)
    affected$moe = substring(affected$first, 5,6)
    affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
    affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
    affected$age_event = affected$age_event1 + affected$month_event1
    affected$first = NULL
    affected$yoe = NULL 
    affected$moe = NULL
    affected$month_event1 = NULL 
    affected$age_event1 = NULL
    
    healthy = dat1[-which(dat1$Sample_Name %in% tmp$id),]
    healthy$Event = 0
    healthy$age_event = 0 
    affected$id.y <- NULL
    healthy$id <- NULL
    names(affected)[names(affected)=="id"] <- "Sample_Name"
    cox = rbind(affected, healthy)
    
    cox$age_death = 0
    cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
    cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
    cox$tte = cox$age_at_event - cox$Age
    cox$tte = as.numeric(cox$tte)
    cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
    cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
    cox$Event = as.numeric(cox$Event)
    cox$tte<-as.numeric(cox$tte)
    
    mod = coxph(Surv(event = cox$Event, cox$tte) ~ scale(cox[,j]) + cox$Age + factor(cox$Female))
    output_hazard[i,1+k[[j]]] <- as.character(names[[i]])
    output_hazard[i, 2:4+k[[j]]]<-exp(cbind(coef(mod), confint(mod)))[1,1:3]
    output_hazard[i,5+k[[j]]] <- summary(mod)$coefficients[1,5]
    output_hazard[i,6+k[[j]]] <- length(which(cox$Event %in% 1))
    output_hazard[i,7+k[[j]]] <- length(which(cox$Event %in% 0))
    
  } 
} 
## Create List of Remaining Dataframes - Cancer  

names(d1)[names(d1) == "breast_cancer_Y"] <- "Breast.Cancer"
names(d1)[names(d1) == "bowel_cancer_Y"] <- "Bowel.Cancer"
names(d1)[names(d1) == "lung_cancer_Y"] <- "Lung.Cancer"

mat_hazard_cancer <- matrix(nrow=length(my.list.cancer),ncol=7*6)
output_hazard_cancer<- as.data.frame(mat_hazard_cancer)

my.list.cancer = list(Bowel.Cancer,Bowel.Cancer,Lung.Cancer)
names = list("Breast.Cancer","Bowel.Cancer","Lung.Cancer") 
names(my.list.cancer) <- names 

l.cancer=lapply(my.list.cancer, "[", 1:6)
smr=lapply(my.list.cancer, "[", c(1,4))
smr1=lapply(smr,subset, smr==1)

for(j in clock){ 
  for(i in 1:length(l.cancer)){ 
    
    tmp <- l.cancer[[i]]
    smr_tmp <-  smr1[[i]]
    
    ## Exclude Indiviudals who Reported Disease at Study Baseline 
    ## Exclude Indiviudals Reported on Cancer Inpatient Day Visit List - might not have had cancer 
    
    dat1= d1[-which(d1[,names[[i]]] %in% 1),]
    dat1=dat1[-which(dat1$Sample_Name %in% smr_tmp$id),]
    tmp1 = tmp[which(tmp$id %in% dat1$Sample_Name),]
    
    ## Obtain Age of Onset 
    affected = dat1[which(dat1$Sample_Name %in% tmp1$id),] 
    age_onset = tmp[,c("dt", "id")]
    affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
    affected$Event = 1
    affected$yoe = substring(affected$dt, 1, 4)
    affected$moe = substring(affected$dt, 5,6)
    affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
    affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
    affected$age_event = affected$age_event1 + affected$month_event1
    affected$dt = NULL
    affected$yoe = NULL 
    affected$moe = NULL
    affected$month_event1 = NULL 
    affected$age_event1 = NULL
    
    healthy = dat1[-which(dat1$Sample_Name %in% tmp$id),]
    healthy$Event = 0
    healthy$age_event = 0 
    affected$id.y <- NULL
    healthy$id <- NULL
    names(affected)[names(affected)=="id"] <- "Sample_Name"
    cox = rbind(affected, healthy)
    
    cox$age_death = 0
    cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
    cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
    cox$tte = cox$age_at_event - cox$Age
    cox$tte = as.numeric(cox$tte)
    cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
    cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
    cox$Event = as.numeric(cox$Event)
    cox$tte<-as.numeric(cox$tte)
    
    mod = coxph(Surv(event = cox$Event, cox$tte) ~ scale(cox[,j]) + cox$Age + factor(cox$Female))
    output_hazard_cancer[i,1+k[[j]]] <- as.character(names[[i]])
    output_hazard_cancer[i, 2:4+k[[j]]]<-exp(cbind(coef(mod), confint(mod)))[1,1:3]
    output_hazard_cancer[i,5+k[[j]]] <- summary(mod)$coefficients[1,5]
    output_hazard_cancer[i,6+k[[j]]] <- length(which(cox$Event %in% 1))
    output_hazard_cancer[i,7+k[[j]]] <- length(which(cox$Event %in% 0))
    
  } 
}

comb = rbind(output_hazard_AD, output_hazard)
comb = rbind(comb, output_hazard_cancer)


## FULLY-ADJUSTED MODEL 

## Run Incidence of Alzheimer's Disease Analysis - Restricted to Participants > 60 

d1_AD <- d1[d1$Age >= 60,]

mat_hazard_ad <- matrix(nrow=1,ncol=7*6)
output_hazard_AD<- as.data.frame(mat_hazard_ad)

for(j in clock){ 
  dat1= d1_AD[-which(d1_AD[,"AD"] %in% 1),]
  tmp1 = AD[which(AD$id %in% dat1$Sample_Name),]
  
  ## Obtain Age of Onset 
  affected = dat1[which(dat1$Sample_Name %in% tmp1$id),] 
  age_onset = AD[,c("first", "id")]
  affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
  affected$Event = 1
  affected$yoe = substring(affected$first, 1, 4)
  affected$moe = substring(affected$first, 5,6)
  affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
  affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
  affected$age_event = affected$age_event1 + affected$month_event1
  affected$first = NULL
  affected$yoe = NULL 
  affected$moe = NULL
  affected$month_event1 = NULL 
  affected$age_event1 = NULL
  
  healthy = dat1[-which(dat1$Sample_Name %in% AD$id),]
  healthy$Event = 0
  healthy$age_event = 0 
  affected$id.y <- NULL
  healthy$id <- NULL
  names(affected)[names(affected)=="id"] <- "Sample_Name"
  cox = rbind(affected, healthy)
  
  cox$age_death = 0
  cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
  cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
  cox$tte = cox$age_at_event - cox$Age
  cox$tte = as.numeric(cox$tte)
  cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
  cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
  cox$Event = as.numeric(cox$Event)
  cox$tte<-as.numeric(cox$tte)
  
  mod = coxph(Surv(event = cox$Event, cox$tte) ~ scale(cox[,j]) + cox$Age + factor(cox$Female) + cox$units_usual + cox$pack_years + cox$simd + cox$EA + cox$bmi)
  output_hazard_AD[1,1+k[[j]]] <- as.character("Alzheimer's Disease")
  output_hazard_AD[1, 2:4+k[[j]]]<-exp(cbind(coef(mod), confint(mod)))[1,1:3]
  output_hazard_AD[1,5+k[[j]]] <- summary(mod)$coefficients[1,5]
  output_hazard_AD[1,6+k[[j]]] <- length(which(cox$Event %in% 1))
  output_hazard_AD[1,7+k[[j]]] <- length(which(cox$Event %in% 0))
}

## Create List of Remaining Dataframes - not Cancer  

my.list = list(COPD,Depression,Heart.Disease,Stroke,Diabetes,Pain)
names = list("COPD","Depression","Heart.Disease","Stroke","Diabetes","Pain")
names(my.list) <- names 

l=lapply(my.list, "[", c(1:5))

names(d1)[names(d1) == "COPD_Y"] <- "COPD"
names(d1)[names(d1) == "depression_Y"] <- "Depression"
names(d1)[names(d1) == "heart_disease_Y"] <- "Heart.Disease"
names(d1)[names(d1) == "stroke_Y"] <- "Stroke"
names(d1)[names(d1) == "diabetes_Y"] <- "Diabetes"

mat_hazard <- matrix(nrow=length(my.list),ncol=7*6)
output_hazard <- as.data.frame(mat_hazard)

## Loop of Logistic Regression Models - Cross-Sectional Associations

for(j in clock){
  for(i in 1:length(l)){ 
    
    tmp <- l[[i]]
    
    ## Exclude Indiviudals who Reported Disease at Study Baseline  
    dat1= d1[-which(d1[,names[[i]]] %in% 1),]
    tmp1 = tmp[which(tmp$id %in% dat1$Sample_Name),]
    
    ## Obtain Age of Onset 
    affected = dat1[which(dat1$Sample_Name %in% tmp1$id),] 
    age_onset = tmp[,c("first", "id")]
    affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
    affected$Event = 1
    affected$yoe = substring(affected$first, 1, 4)
    affected$moe = substring(affected$first, 5,6)
    affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
    affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
    affected$age_event = affected$age_event1 + affected$month_event1
    affected$first = NULL
    affected$yoe = NULL 
    affected$moe = NULL
    affected$month_event1 = NULL 
    affected$age_event1 = NULL
    
    healthy = dat1[-which(dat1$Sample_Name %in% tmp$id),]
    healthy$Event = 0
    healthy$age_event = 0 
    affected$id.y <- NULL
    healthy$id <- NULL
    names(affected)[names(affected)=="id"] <- "Sample_Name"
    cox = rbind(affected, healthy)
    
    cox$age_death = 0
    cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
    cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
    cox$tte = cox$age_at_event - cox$Age
    cox$tte = as.numeric(cox$tte)
    cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
    cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
    cox$Event = as.numeric(cox$Event)
    cox$tte<-as.numeric(cox$tte)
    
    mod = coxph(Surv(event = cox$Event, cox$tte) ~ scale(cox[,j]) + cox$Age + factor(cox$Female) + cox$units_usual + cox$EA + cox$simd + cox$pack_years + cox$bmi)
    output_hazard[i,1+k[[j]]] <- as.character(names[[i]])
    output_hazard[i, 2:4+k[[j]]]<-exp(cbind(coef(mod), confint(mod)))[1,1:3]
    output_hazard[i,5+k[[j]]] <- summary(mod)$coefficients[1,5]
    output_hazard[i,6+k[[j]]] <- length(which(cox$Event %in% 1))
    output_hazard[i,7+k[[j]]] <- length(which(cox$Event %in% 0))
    
  } 
}
## Create List of Remaining Dataframes - Cancer  

names(d1)[names(d1) == "breast_cancer_Y"] <- "Breast.Cancer"
names(d1)[names(d1) == "bowel_cancer_Y"] <- "Bowel.Cancer"
names(d1)[names(d1) == "lung_cancer_Y"] <- "Lung.Cancer"

mat_hazard_cancer <- matrix(nrow=length(my.list.cancer),ncol=7*6)
output_hazard_cancer<- as.data.frame(mat_hazard_cancer)

my.list.cancer = list(Bowel.Cancer,Bowel.Cancer,Lung.Cancer)
names = list("Breast.Cancer","Bowel.Cancer","Lung.Cancer") 
names(my.list.cancer) <- names 

l.cancer=lapply(my.list.cancer, "[", 1:6)
smr=lapply(my.list.cancer, "[", c(1,4))
smr1=lapply(smr,subset, smr==1)

for(j in clock){
  for(i in 1:length(l.cancer)){ 
    
    tmp <- l.cancer[[i]]
    smr_tmp <-  smr1[[i]]
    
    ## Exclude Indiviudals who Reported Disease at Study Baseline 
    ## Exclude Indiviudals Reported on Cancer Inpatient Day Visit List - might not have had cancer 
    
    dat1= d1[-which(d1[,names[[i]]] %in% 1),]
    dat1=dat1[-which(dat1$Sample_Name %in% smr_tmp$id),]
    tmp1 = tmp[which(tmp$id %in% dat1$Sample_Name),]
    
    ## Obtain Age of Onset 
    affected = dat1[which(dat1$Sample_Name %in% tmp1$id),] 
    age_onset = tmp[,c("dt", "id")]
    affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
    affected$Event = 1
    affected$yoe = substring(affected$dt, 1, 4)
    affected$moe = substring(affected$dt, 5,6)
    affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
    affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
    affected$age_event = affected$age_event1 + affected$month_event1
    affected$dt = NULL
    affected$yoe = NULL 
    affected$moe = NULL
    affected$month_event1 = NULL 
    affected$age_event1 = NULL
    
    healthy = dat1[-which(dat1$Sample_Name %in% tmp$id),]
    healthy$Event = 0
    healthy$age_event = 0 
    affected$id.y <- NULL
    healthy$id <- NULL
    names(affected)[names(affected)=="id"] <- "Sample_Name"
    cox = rbind(affected, healthy)
    
    cox$age_death = 0
    cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
    cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
    cox$tte = cox$age_at_event - cox$Age
    cox$tte = as.numeric(cox$tte)
    cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
    cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
    cox$Event = as.numeric(cox$Event)
    cox$tte<-as.numeric(cox$tte)
    
    mod = coxph(Surv(event = cox$Event, cox$tte) ~ scale(cox[,j]) + cox$Age + factor(cox$Female) + cox$units_usual + cox$EA + cox$simd + cox$pack_years + cox$bmi)
    output_hazard_cancer[i,1+k[[j]]] <- as.character(names[[i]])
    output_hazard_cancer[i, 2:4+k[[j]]]<-exp(cbind(coef(mod), confint(mod)))[1,1:3]
    output_hazard_cancer[i,5+k[[j]]] <- summary(mod)$coefficients[1,5]
    output_hazard_cancer[i,6+k[[j]]] <- length(which(cox$Event %in% 1))
    output_hazard_cancer[i,7+k[[j]]] <- length(which(cox$Event %in% 0))
    
  } 
} 

comb_full = rbind(output_hazard_AD, output_hazard)
comb_full = rbind(comb_full, output_hazard_cancer)

comb <- comb[,c(1,6,7,2:5,(c(1,6,7,2:5)+7),(c(1,6,7,2:5)+14),(c(1,6,7,2:5)+21),(c(1,6,7,2:5)+28), (c(1,6,7,2:5)+35))]
names(comb) <- rep(c("Trait", "no. of cases", "no. of controls", "Hazard Ratio", "LCI", "UCI", "P"),6)

comb_full <- comb_full[,c(1,6,7,2:5,(c(1,6,7,2:5)+7),(c(1,6,7,2:5)+14),(c(1,6,7,2:5)+21),(c(1,6,7,2:5)+28), (c(1,6,7,2:5)+35))]
names(comb_full) <- rep(c("Trait", "no. of cases", "no. of controls", "Hazard Ratio", "LCI", "UCI", "P"),6)


## Final Files 

View(comb)
View(comb_full)


write.csv(comb, "C:/Users/s1777309/Desktop/PhD/Year 2/Epigenetic Clocks WHO/New_Summary/Incidence_Basic.csv",row.names = F)
write.csv(comb_full, "C:/Users/s1777309/Desktop/PhD/Year 2/Epigenetic Clocks WHO/New_Summary/Incidence_FullyAdjusted.csv",row.names = F)
