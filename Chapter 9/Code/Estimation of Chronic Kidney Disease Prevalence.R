setwd("./Epigenetic Clocks WHO") 

## Discovery Dataset 

d1 <- read.csv("Discovery_Phenotype_1.csv", check.names = F)
creat <- read.csv("Corr_Creatinine.csv")
creat = creat[,c(1,5)]
d1 = merge(d1, creat, by.x = "Sample_Name", by.y = "id", all.x = T)
names(d1)[92] <- "Creat"

females = d1[which(d1$sex %in% 1),]
males = d1[which(d1$sex %in% 0),] 
na = d1[which(d1$sex %in% NA),]

females_low <- females[which(females$Creat <= 62),] 
females_high <- females[which(females$Creat > 62),] 
females_na <- females[which(females$Creat %in% NA),] 

males_low <- males[which(males$Creat <= 80),] 
males_high <- males[which(males$Creat > 80),] 
males_na <- males[which(males$Creat %in% NA),] 

females_low$eGFR <- 141*((females_low$Creat/61.9)^-0.329)*(0.993^females_low$age)*1.018
females_high$eGFR <- 141*((females_high$Creat/61.9)^-1.209)*(0.993^females_high$age)*1.018
males_low$eGFR <- 141*((males_low$Creat/79.6)^-0.411)*(0.993^males_low$age)
males_high$eGFR <-141*((males_high$Creat/79.6)^-1.209)*(0.993^males_high$age)
females_na$eGFR <- NA
males_na$eGFR <- NA 
na$eGFR <- NA 
females_1 = rbind(females_low,females_high)
females_1 = rbind(females_1, females_na)
males_1 = rbind(males_low,males_high)
males_1 = rbind(males_1, males_na)
all = rbind(males_1, females_1)
all = rbind(all, na)
d1<-all
d1$CKD = ifelse(d1$eGFR < 60, 1, 0)


names(d1)[names(d1) == "CKD"] <- "Chronic Kidney Disease"
names(d1)[names(d1) == "eGFR"] <- "Estimated Glomerular Filtration Rate"

write.csv(d1, filename = "Discovery_CKD.csv")

## Replication Dataset 

d1 <- read.csv("Replication_Phenotype_1.csv", check.names = F)
creat <- read.csv("Corr_Creatinine.csv")
creat = creat[,c(1,5)]
d1 = merge(d1, creat, by.x = "Sample_Name", by.y = "id", all.x = T)
names(d1)[92] <- "Creat"

females = d1[which(d1$sex %in% 1),]
males = d1[which(d1$sex %in% 0),] 
na = d1[which(d1$sex %in% NA),]

females_low <- females[which(females$Creat <= 62),] 
females_high <- females[which(females$Creat > 62),] 
females_na <- females[which(females$Creat %in% NA),] 

males_low <- males[which(males$Creat <= 80),] 
males_high <- males[which(males$Creat > 80),] 
males_na <- males[which(males$Creat %in% NA),] 

females_low$eGFR <- 141*((females_low$Creat/61.9)^-0.329)*(0.993^females_low$age)*1.018
females_high$eGFR <- 141*((females_high$Creat/61.9)^-1.209)*(0.993^females_high$age)*1.018
males_low$eGFR <- 141*((males_low$Creat/79.6)^-0.411)*(0.993^males_low$age)
males_high$eGFR <-141*((males_high$Creat/79.6)^-1.209)*(0.993^males_high$age)
females_na$eGFR <- NA
males_na$eGFR <- NA 
na$eGFR <- NA 
females_1 = rbind(females_low,females_high)
females_1 = rbind(females_1, females_na)
males_1 = rbind(males_low,males_high)
males_1 = rbind(males_1, males_na)
all = rbind(males_1, females_1)
all = rbind(all, na)
d1<-all
d1$CKD = ifelse(d1$eGFR < 60, 1, 0)


names(d1)[names(d1) == "CKD"] <- "Chronic Kidney Disease"
names(d1)[names(d1) == "eGFR"] <- "Estimated Glomerular Filtration Rate"

write.csv(d1, filename = "Replication_Dataset.csv")
