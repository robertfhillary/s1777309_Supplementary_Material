setwd("./Epigenetic Clocks WHO") 

## Installing Requisite Packages 

if(!require(data.table)){
  install.packages("data.table")
}

library(data.table)


## Summary - Discovery Dataset 

d1 <- read.csv("Discovery_Phenotype_1.csv", check.names = F)
colnames(d1) = gsub(".", " ", colnames(d1), fixed = T)
names(d1)[names(d1) == "Alzheimers Disease   Maternal"] <- "Alzheimers Disease Maternal" 
names(d1)[names(d1) == "Alzheimers Disease   Paternal"] <- "Alzheimers Disease Paternal" 
names(d1)[names(d1) == "Creat"] <- "Creatinine"

cont <- c("age", "Creatinine", "SIMD", "Neuroticism", "Systolic Pressure", "Forced Expiratory Flow", "Forced Expiratory Volume", "Forced Vital Capacity", 
          "Body Mass Index", "Waist:Hip Ratio", "Diastolic Pressure", "g", "gf", "Average Heart Rate", 
          "QT Interval (corrected for heart rate)", "Total Cholesterol", "HDL Cholesterol", "Glucose", "Pack Years")   

categ <- c("sex", "Chronic Kidney Disease", "Back Pain", "Neck Pain", "Stroke", "Heart Disease",         
           "Depression", "Diabetes", "COPD", "Lung Cancer", "Bowel Cancer",   
           "Breast Cancer", "Alzheimers Disease Paternal", "Alzheimers Disease Maternal","SCID Depression")  

cat_sum <- data.frame(Variable=character(), n=double(), n_events = double(), stringsAsFactors=FALSE)

for(i in 1:length(categ)){
  
  var = categ[i]
  n = length(which(!is.na(d1[,categ[i]])))  
  n_events = length(which(d1[,categ[i]] %in% 1))

  cat_sum[i,] <- c(var, n, n_events)
  
}

cont_sum <- data.frame(Variable=character(), n=double(), mean=double(), sd=double(), stringsAsFactors=FALSE)

  for(i in 1:length(cont)){
    
    var = cont[i]
    n = length(which(!is.na(d1[,cont[i]])))  
    mean = mean(d1[,cont[i]], na.rm = T)
    sd = sd(d1[,cont[i]], na.rm = T)
    
    cont_sum[i,] <- c(var, n, signif(mean,3), signif(sd, 3))
    
  }

ea <- data.frame(Variable = "Educational Attainment", 
                 Median = median(d1$EA, na.rm = T),
                 IQR = IQR(d1$EA, na.rm =T))

dis = rbind(setDT(cont_sum), setDT(ea), fill=TRUE)


dis = rbind(setDT(dis), setDT(cat_sum), fill=TRUE)


## Summary - Replication Dataset 

d1 <- read.csv("Replication_Phenotype_1.csv", check.names = F)
names(d1)[names(d1) == "Creat"] <- "Creatinine"
names(d1)[names(d1) == "Waist Hip Ratio"] <- "Waist:Hip Ratio"

cont <- c("age", "Creatinine", "SIMD", "Neuroticism", "Systolic Pressure", "Forced Expiratory Flow", "Forced Expiratory Volume", "Forced Vital Capacity", 
          "Body Mass Index", "Waist:Hip Ratio", "Diastolic Pressure", "g", "gf", "Average Heart Rate", 
          "QT Interval (corrected for heart rate)", "Total Cholesterol", "HDL Cholesterol", "Glucose", "Pack Years")   

categ <- c("sex", "Chronic Kidney Disease", "Back Pain", "Neck Pain", "Stroke", "Heart Disease",         
           "Depression", "Diabetes", "COPD", "Lung Cancer", "Bowel Cancer",   
           "Breast Cancer", "Alzheimers Disease Paternal", "Alzheimers Disease Maternal","SCID Depression")  

cat_sum <- data.frame(Variable=character(), n=double(), n_events = double(), stringsAsFactors=FALSE)


for(i in 1:length(categ)){
  
  var = categ[i]
  n = length(which(!is.na(d1[,categ[i]])))  
  n_events = length(which(d1[,categ[i]] %in% 1))
  
  cat_sum[i,] <- c(var, n, n_events)
  
}

cont_sum <- data.frame(Variable=character(), n=double(), Mean=double(), SD=double(), stringsAsFactors=FALSE)

for(i in 1:length(cont)){
  
  var = cont[i]
  n = length(which(!is.na(d1[,cont[i]])))  
  mean = mean(d1[,cont[i]], na.rm = T)
  sd = sd(d1[,cont[i]], na.rm = T)
  
  cont_sum[i,] <- c(var, n, signif(mean,3), signif(sd, 3))
  
}

ea <- data.frame(Variable = "Educational Attainment", 
                 Median = median(d1$EA, na.rm = T),
                 IQR = IQR(d1$EA, na.rm =T))

rep = rbind(setDT(cont_sum), setDT(ea), fill=TRUE)


rep = rbind(setDT(rep), setDT(cat_sum), fill=TRUE)


