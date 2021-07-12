setwd(",/Epigenetic Clocks WHO")

## Installing Requisite Packages 

if(!require(survival)){
  install.packages("survival")
}

if(!require(ggplot2)){
  install.packages("ggplot2")
}


library(survival)
library(ggplot2)


## DISCOVERY STAGE OF ANALYSIS 

## Clean-up of data 

d1 <- read.csv("Replication_Phenotype_1.csv", check.names = F)
colnames(d1) = gsub(".", " ", colnames(d1), fixed = T)
names(d1)[names(d1) == "Alzheimers Disease   Maternal"] <- "Alzheimers Disease Maternal" 
names(d1)[names(d1) == "Alzheimers Disease   Paternal"] <- "Alzheimers Disease Paternal" 
names(d1)[names(d1) == "Creat"] <- "Creatinine"
names(d1)[names(d1) == "QT Interval  corrected for heart rate "] <- "QT Interval (corrected for heart rate)"
names(d1)[names(d1) == "Waist Hip Ratio"] <- "Waist:Hip Ratio"
clock <- c("IEAA","EEAA","AgeAccelGrim","AgeAccelPheno","DNAmTLAdjAge", "Dunedin_PoAm") 

categ <- c("Chronic Kidney Disease", "Back Pain", "Neck Pain", "Stroke", "Heart Disease",         
           "Depression", "Diabetes", "COPD", "Lung Cancer", "Bowel Cancer",   
           "Breast Cancer", "Alzheimers Disease Paternal", "Alzheimers Disease Maternal","SCID Depression")  


## All-Cause Mortality

out_surv <- data.frame(Clock=character(), n=double(), n_event=double(), HR=double(), LCI=double(), UCI=double(), P=double(), stringsAsFactors=FALSE)

for(j in 1:length(clock)){
  clk = clock[j]
  
  model <- summary(coxph(Surv(tte, event) ~ scale(d1[,clock[j]])+ age + sex + SIMD + EA + d1$units_usual + 
                           d1$`Pack Years` + d1$`Body Mass Index`, data=d1))
  
  n = model$n
  n_event = model$nevent
  HR = model$coefficients[1,2]
  LCI = exp(model$coefficients[1,1] - 1.96*model$coefficients[1,3])
  UCI = exp(model$coefficients[1,1] + 1.96*model$coefficients[1,3])
  P = model$coefficients[1,5]
  
  out_surv[j,] <- c(clk, n, n_event, signif(HR, 3), signif(LCI, 3), signif(UCI, 3), signif(P, 2))
  
}

test <- out_surv
test$Clock <- factor(test$Clock, levels=rev(c("AgeAccelGrim", "Dunedin_PoAm", "AgeAccelPheno", "EEAA", "IEAA", "DNAmTLAdjAge")))
test$upper <- as.numeric(test$UCI)
test$lower <- as.numeric(test$LCI)

p = ggplot(test,aes(y=as.numeric(test$HR), x=test$Clock)) + 
  geom_point(position=position_dodge(width=0.5), size = 3, colour = "gold2") +
  geom_errorbar(aes(ymin = test$lower, ymax = test$upper),
                position = position_dodge(0.5), width = 0.05,
                colour = "black") +
  ylab("Hazard Ratio")+ xlab ("") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), legend.position = "none") +
  ggtitle("Survival - Discovery") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  coord_flip()

print(p)


## Logistic Regression Analysis of Continuous Variables 

out_cat <- data.frame(Clock=character(), Variable=character(), n=double(), n_events=double(), OR=double(), LCI=double(), UCI=double(),P=double(), stringsAsFactors=FALSE)

for(j in 1:length(clock)){
  for(i in 1:length(categ)){
    
    clk <- clock[j]
    var = categ[i]
    p = summary(glm(d1[,categ[i]] ~ scale(d1[,clock[j]]) + age + sex + SIMD + EA + d1$units_usual + d1$`Pack Years` 
                    + d1$`Body Mass Index`, data=d1, family="binomial"))$coefficients[2,4]
    n = nobs(glm(d1[,categ[i]] ~ scale(d1[,clock[j]]) + age + sex + SIMD + EA + d1$units_usual + 
                   d1$`Pack Years` + d1$`Body Mass Index`, data=d1, family="binomial"))
    n_events = table(d1[,categ[i]])[2]
    logodds <- summary(glm(d1[,categ[i]] ~ scale(d1[,clock[j]])  + age + sex + SIMD + EA + d1$units_usual +
                             d1$`Pack Years` + d1$`Body Mass Index`, data=d1, family="binomial"))$coefficients[2,1]
    se <- summary(glm(d1[,categ[i]] ~ scale(d1[,clock[j]])  + age + sex + SIMD + EA + d1$units_usual
                      + d1$`Pack Years` + d1$`Body Mass Index`, data=d1, family="binomial"))$coefficients[2,2]
    
    or <- logodds
    lower <- logodds - 1.96*se
    upper <- logodds + 1.96*se 
    
    out_cat[i + (length(categ)*(j-1)),] <- c(clk, var, n, n_events, signif(or,3), signif(lower, 3), signif(upper, 3), signif(p, 2))
    
  }
}


## Linear Regression Analysis of Continuous Variables 

cont <- c("Neuroticism", "Systolic Pressure", "Forced Expiratory Flow", "Forced Expiratory Volume", "Forced Vital Capacity", 
          "Waist:Hip Ratio", "Diastolic Pressure", "g", "gf", "Average Heart Rate", 
          "QT Interval (corrected for heart rate)", "Total Cholesterol", "HDL Cholesterol", "Glucose")   

out_cont <- data.frame(Clock=character(), Variable=character(),n=double(),Beta=double(),SE=double(),P=double(), stringsAsFactors=FALSE)


for(j in 1:length(clock)){
  for(i in 1:length(cont)){
    clk = clock[j]
    var = cont[i]
    coeft = summary(lm(scale(d1[,cont[i]]) ~ scale(d1[,clock[j]]) + age + sex + SIMD + EA + d1$units_usual + d1$`Pack Years` + d1$`Body Mass Index`, data=d1))$coefficients[2,c(1,2,4)]
    n = nobs(lm(scale(d1[,cont[i]]) ~ scale(d1[,clock[j]]) + age + sex + SIMD + EA + d1$units_usual + d1$`Pack Years` + d1$`Body Mass Index`, data=d1))
    out_cont[i + (length(cont)*(j-1)),] <- c(clk, var, n, signif(coeft, 2))
  }
}

creat <- data.frame(Clock=character(), Variable=character(),n=double(),Beta=double(),SE=double(),P=double(), stringsAsFactors=FALSE)

cont = "Creatinine"
for(j in 1:length(clock)){
  for(i in 1:length(cont)){
    clk = clock[j]
    var = cont[i]
    coeft = summary(lm(scale(d1[,cont[i]]) ~ scale(d1[,clock[j]]) + age + SIMD + EA + d1$units_usual + d1$`Pack Years` + d1$`Body Mass Index`, data=d1))$coefficients[2,c(1,2,4)]
    n = nobs(lm(scale(d1[,cont[i]]) ~ scale(d1[,clock[j]]) + age + SIMD + EA + d1$units_usual + d1$`Pack Years` + d1$`Body Mass Index`, data=d1))
    creat[i + (length(cont)*(j-1)),] <- c(clk, var, n, signif(coeft, 2))
  }
}

SIMD <- data.frame(Clock=character(), Variable=character(),n=double(),Beta=double(),SE=double(),P=double(), stringsAsFactors=FALSE)

cont = "SIMD"
for(j in 1:length(clock)){
  for(i in 1:length(cont)){
    clk = clock[j]
    var = cont[i]
    coeft = summary(lm(scale(d1[,cont[i]]) ~ scale(d1[,clock[j]]) + age  + EA + d1$units_usual + d1$`Pack Years` + d1$`Body Mass Index`, data=d1))$coefficients[2,c(1,2,4)]
    n = nobs(lm(scale(d1[,cont[i]]) ~ scale(d1[,clock[j]]) + age + EA + d1$units_usual + d1$`Pack Years` + d1$`Body Mass Index`, data=d1))
    SIMD[i + (length(cont)*(j-1)),] <- c(clk, var, n, signif(coeft, 2))
  }
}

PackYRS <- data.frame(Clock=character(), Variable=character(),n=double(),Beta=double(),SE=double(),P=double(), stringsAsFactors=FALSE)

cont = "Pack Years"
for(j in 1:length(clock)){
  for(i in 1:length(cont)){
    clk = clock[j]
    var = cont[i]
    coeft = summary(lm(scale(d1[,cont[i]]) ~ scale(d1[,clock[j]]) + age + sex + SIMD + EA + d1$units_usual+ d1$`Body Mass Index`, data=d1))$coefficients[2,c(1,2,4)]
    n = nobs(lm(scale(d1[,cont[i]]) ~ scale(d1[,clock[j]]) + age + sex + SIMD + EA + d1$units_usual + d1$`Body Mass Index`, data=d1))
    PackYRS[i + (length(cont)*(j-1)),] <- c(clk, var, n, signif(coeft, 2))
  }
}

BMI <- data.frame(Clock=character(), Variable=character(),n=double(),Beta=double(),SE=double(),P=double(), stringsAsFactors=FALSE)

cont = "Body Mass Index"
for(j in 1:length(clock)){
  for(i in 1:length(cont)){
    clk = clock[j]
    var = cont[i]
    coeft = summary(lm(scale(d1[,cont[i]]) ~ scale(d1[,clock[j]]) + age + sex + SIMD + EA + d1$units_usual+ d1$`Pack Years`, data=d1))$coefficients[2,c(1,2,4)]
    n = nobs(lm(scale(d1[,cont[i]]) ~ scale(d1[,clock[j]]) + age + sex + SIMD + EA + d1$units_usual + d1$`Pack Years`, data=d1))
    BMI[i + (length(cont)*(j-1)),] <- c(clk, var, n, signif(coeft, 2))
  }
}



out_cont = rbind(out_cont, creat)
out_cont = rbind(out_cont, SIMD)
out_cont = rbind(out_cont, PackYRS)
out_cont = rbind(out_cont, BMI)


out_cont$LCI <- as.numeric(out_cont$Beta) - 1.96*as.numeric(out_cont$SE)
out_cont$UCI <- as.numeric(out_cont$Beta) + 1.96*as.numeric(out_cont$SE)
out_cont  =out_cont[,c(1:5,7,8,6)]


######### Plots for Individual Disease/Variable Pairings
pos <- position_dodge(0.4)
basic_con <- out_cont
basic_cat <- out_cat
basic_cat$n_events <- NULL
basic_con$SE <- NULL
names(basic_con)[4] <- "test"
names(basic_cat)[4] <- "test"
basic_con$test <- as.numeric(basic_con$test)
basic_cat$test <- as.numeric(basic_cat$test)

## Cardiovascular Diseases

##Stroke
cv_con = basic_con[which(basic_con$Variable %in% c("Systolic Pressure", "Diastolic Pressure", "HDL Cholesterol", "Total Cholesterol")),] 
cv_cat = basic_cat[which(basic_cat$Variable %in% "Stroke"),]
cv_con$Type = "Phenotype"
cv_cat$Type = "Disease"
stroke = rbind(cv_cat, cv_con)


stroke$Clock <- factor(stroke$Clock, levels = rev(c("AgeAccelGrim", "Dunedin_PoAm", "AgeAccelPheno", "EEAA", "IEAA", "DNAmTLAdjAge"))) 
stroke$Variable = factor(stroke$Variable, levels=c("Systolic Pressure", "Diastolic Pressure", "HDL Cholesterol","Total Cholesterol","Stroke"))
stroke_graph <- ggplot(stroke,aes(y=as.numeric(stroke$test), 
                                  x=Clock, group = Variable)) +   
  geom_point(data=stroke, size = 3, 
             aes(colour = Variable, shape = factor(Type)), show.legend =T, 
             position = pos) +
  geom_errorbar(aes(ymin = as.numeric(stroke$LCI), 
                    ymax= as.numeric(stroke$UCI)), position = pos,
                colour  ="black", width = 0.05) +
   xlab ("")+ ylab("Standardised Beta/Log Odds \n [95% Confidence Interval]") + 
  geom_hline(yintercept = 0, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(size = 12))+
   guides(colour = guide_legend(reverse = T)) +
  labs(shape = "Type")



stroke_graph + ggtitle("Stroke and Associated Phenotypes") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title.align = 0.5) +  theme(legend.title.align = 0.5) +coord_flip()



##Heart Disease
cv_con = basic_con[which(basic_con$Variable %in% c("Average Heart Rate", "QT Interval (corrected for heart rate)", "Body Mass Index", "Waist:Hip Ratio")),] 
cv_cat = basic_cat[which(basic_cat$Variable %in% "Heart Disease"),]
cv_con$Type = "Phenotype"
cv_cat$Type = "Disease"
cv = rbind(cv_cat, cv_con)

cv$Clock <- factor(cv$Clock, levels = rev(c("AgeAccelGrim", "Dunedin_PoAm", "AgeAccelPheno", "EEAA", "IEAA", "DNAmTLAdjAge"))) 
cv$Variable = factor(cv$Variable, 
                     levels=c("Average Heart Rate", 
                              "QT Interval (corrected for heart rate)", 
                              "Body Mass Index", "Waist:Hip Ratio",
                              "Heart Disease"))

cv_graph <- ggplot(cv,aes(y=as.numeric(cv$test), 
                          x=Clock, group = Variable)) +   
  geom_point(data=cv, size = 3, 
             aes(colour = Variable, shape = factor(Type)), show.legend =T, 
             position = pos) +
  geom_errorbar(aes(ymin = as.numeric(cv$LCI), 
                    ymax= as.numeric(cv$UCI)), position = pos,
                colour  ="black", width = 0.05) + coord_flip() +
  ylab("Standardised Beta/Log Odds \n [95% Confidence Interval]")+ xlab ("")+
  geom_hline(yintercept = 0, linetype = "dotted")+  
  theme(legend.title.align = 0.5) + 
  theme(axis.text.x = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(size = 12))+
  labs(shape = "Type") + guides(colour = guide_legend(reverse = T))  


cv_graph + ggtitle("Heart Disease and Associated Phenotypes") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.title.align = 0.5) 


## Neurological Diseases

##Alzheimer's Disease 
ad_con = basic_con[which(basic_con$Variable %in% c("g", "gf")),] 
ad_cat = basic_cat[which(basic_cat$Variable %in% c("Alzheimers Disease Maternal", "Alzheimers Disease Paternal")),]
ad_con$Type = "Phenotype"
ad_cat$Type = "Disease"
ad = rbind(ad_con, ad_cat)


ad$Clock <- factor(ad$Clock, levels = rev(c("AgeAccelGrim", "Dunedin_PoAm", "AgeAccelPheno", "EEAA", "IEAA", "DNAmTLAdjAge"))) 
ad$Variable = factor(ad$Variable, levels=c("g", "gf", "Alzheimers Disease Maternal", "Alzheimers Disease Paternal"))
ad_graph <- ggplot(ad,aes(y=as.numeric(ad$test), 
                          x=Clock, group = Variable)) +   
  geom_point(data=ad, size = 3, 
             aes(colour = Variable, shape = factor(Type)), show.legend =T, 
             position = pos) +
  geom_errorbar(aes(ymin = as.numeric(ad$LCI), 
                    ymax= as.numeric(ad$UCI)), position = pos,
                colour  ="black", width = 0.05) +
  ylab("Standardised Beta/Log Odds \n [95% Confidence Interval]")+ xlab ("")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(size = 12))+
  guides(colour = guide_legend(reverse = T)) +
  labs(shape = "Type")



ad_graph + ggtitle("Alzheimer's Disease and Associated Phenotype") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title.align = 0.5) + coord_flip()


##Depression 
dep_con = basic_con[which(basic_con$Variable %in% c("SIMD", "Neuroticism")),] 
dep_cat = basic_cat[which(basic_cat$Variable %in% c("SCID Depression", "Depression")),]
dep_con$Type = "Phenotype"
dep_cat$Type = "Disease"
dep = rbind(dep_con, dep_cat)


dep$Clock <- factor(dep$Clock, levels = rev(c("AgeAccelGrim", "Dunedin_PoAm", "AgeAccelPheno", "EEAA", "IEAA", "DNAmTLAdjAge"))) 
dep$Variable = factor(dep$Variable, levels=c("SIMD", "Neuroticism", "SCID Depression", "Depression"))
dep_graph <- ggplot(dep,aes(y=as.numeric(dep$test), 
                            x=Clock, group = Variable)) +   
  geom_point(data=dep, size = 3, 
             aes(colour = Variable, shape = factor(Type)), show.legend =T, 
             position = pos) +
  geom_errorbar(aes(ymin = as.numeric(dep$LCI), 
                    ymax= as.numeric(dep$UCI)), position = pos,
                colour  ="black", width = 0.05) +
  ylab("Standardised Beta/Log Odds \n [95% Confidence Interval]")+ xlab ("")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(size = 12))+
  guides(colour = guide_legend(reverse = T)) +
  labs(shape = "Type")



dep_graph + ggtitle("Depression and Associated Phenotypes") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title.align = 0.5) + coord_flip()

## Pulmonary Diseases 

##COPD 
copd_con = basic_con[which(basic_con$Variable %in% c("Forced Expiratory Flow", "Forced Expiratory Volume", "Forced Vital Capacity")),] 
copd_cat = basic_cat[which(basic_cat$Variable %in% "COPD"),]
copd_con$Type = "Phenotype"
copd_cat$Type = "Disease"
copd = rbind(copd_con, copd_cat)


copd$Clock <- factor(copd$Clock, levels = rev(c("AgeAccelGrim", "Dunedin_PoAm", "AgeAccelPheno", "EEAA", "IEAA", "DNAmTLAdjAge"))) 
copd$Variable = factor(copd$Variable, levels=c("Forced Vital Capacity", "Forced Expiratory Volume", "Forced Expiratory Flow", "COPD"))
copd_graph <- ggplot(copd,aes(y=as.numeric(copd$test), 
                              x=Clock, group = Variable)) +   
  geom_point(data=copd, size = 3, 
             aes(colour = Variable, shape = factor(Type)), show.legend =T, 
             position = pos) +
  geom_errorbar(aes(ymin = as.numeric(copd$LCI), 
                    ymax= as.numeric(copd$UCI)), position = pos,
                colour  ="black", width = 0.05) +
  ylab("Standardised Beta/Log Odds \n [95% Confidence Interval]")+ xlab ("")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(size = 12))+
   guides(colour = guide_legend(reverse = T)) +
  labs(shape = "Type")



copd_graph + ggtitle("COPD and Associated Phenotypes") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title.align = 0.5) + coord_flip()


##Lung Cancer
lc_con = basic_con[which(basic_con$Variable %in% "Pack Years"),] 
lc_cat = basic_cat[which(basic_cat$Variable %in% "Lung Cancer"),]
lc_con$Type = "Phenotype"
lc_cat$Type = "Disease"
lc = rbind(lc_con, lc_cat)


lc$Clock <- factor(lc$Clock, levels = rev(c("AgeAccelGrim", "Dunedin_PoAm", "AgeAccelPheno", "EEAA", "IEAA", "DNAmTLAdjAge"))) 
lc$Variable = factor(lc$Variable, levels=c("Pack Years", "Lung Cancer"))
lc_graph <- ggplot(lc,aes(y=as.numeric(lc$test), 
                          x=Clock, group = Variable)) +   
  geom_point(data=lc, size = 3, 
             aes(colour = Variable, shape = factor(Type)), show.legend =T, 
             position = pos) +
  geom_errorbar(aes(ymin = as.numeric(lc$LCI), 
                    ymax= as.numeric(lc$UCI)), position = pos,
                colour  ="black", width = 0.05) +
  ylab("Standardised Beta/Log Odds \n [95% Confidence Interval]")+ xlab ("")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(size = 12))+
   guides(colour = guide_legend(reverse = T)) + 
  labs(shape = "Type")

lc_graph + ggtitle("Lung Cancer and Associated Phenotype") +
  guides(shape = guide_legend(order=1),
         color = guide_legend(order=2)) +  
  guides(colour = guide_legend(reverse = T)) +
  theme(plot.title = element_text(hjust = 0.5)) +theme(legend.title.align = 0.5) + coord_flip()


## Other Cancers

##Breast Cancer 
bc = basic_cat[which(basic_cat$Variable %in% "Breast Cancer"),]
bc$Type = "Disease"


bc$Clock <- factor(bc$Clock, levels = rev(c("AgeAccelGrim", "Dunedin_PoAm", "AgeAccelPheno", "EEAA", "IEAA", "DNAmTLAdjAge"))) 
bc_graph <- ggplot(bc,aes(y=as.numeric(bc$test), 
                          x=Clock, group = Variable)) +   
  geom_point(data=bc, size = 3, 
             aes(colour = Variable), show.legend =T, 
             position = pos) +
  geom_errorbar(aes(ymin = as.numeric(bc$LCI), 
                    ymax= as.numeric(bc$UCI)), position = pos,
                colour  ="black", width = 0.05) +
  ylab("Standardised Beta/Log Odds \n [95% Confidence Interval]")+ xlab ("")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(size = 12))+
  guides(colour = guide_legend(reverse = T)) 



bc_graph + ggtitle("Breast Cancer") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title.align = 0.5) + coord_flip()


##Bowel Cancer 
bow_c = basic_cat[which(basic_cat$Variable %in% "Bowel Cancer"),]
bow_c$Type = "Disease"


bow_c$Clock <- factor(bow_c$Clock, levels = rev(c("AgeAccelGrim", "Dunedin_PoAm", "AgeAccelPheno", "EEAA", "IEAA", "DNAmTLAdjAge"))) 
bow_c_graph <- ggplot(bow_c,aes(y=as.numeric(bow_c$test), 
                                x=Clock, group = Variable)) +   
  geom_point(data=bow_c, size = 3, 
             aes(colour = Variable), show.legend =T, 
             position = pos) +
  geom_errorbar(aes(ymin = as.numeric(bow_c$LCI), 
                    ymax= as.numeric(bow_c$UCI)), position = pos,
                colour  ="black", width = 0.05) +
  ylab("Standardised Beta/Log Odds \n [95% Confidence Interval]")+ xlab ("")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(size = 12))+
  guides(colour = guide_legend(reverse = T)) 



bow_c_graph + ggtitle("Bowel Cancer") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title.align = 0.5) + coord_flip()

## Renal and Endocrinology Diseases

##Diabetes
diab_con = basic_con[which(basic_con$Variable %in% "Glucose"),] 
diab_cat = basic_cat[which(basic_cat$Variable %in% "Diabetes"),]
diab_con$Type = "Phenotype"
diab_cat$Type = "Disease"
diab = rbind(diab_con, diab_cat)


diab$Clock <- factor(diab$Clock, levels = rev(c("AgeAccelGrim", "Dunedin_PoAm", "AgeAccelPheno", "EEAA", "IEAA", "DNAmTLAdjAge"))) 
diab$Variable = factor(diab$Variable, levels=c("Glucose", "Diabetes"))
diab_graph <- ggplot(diab,aes(y=as.numeric(diab$test), 
                              x=Clock, group = Variable)) +   
  geom_point(data=diab, size = 3, 
             aes(colour = Variable, shape = factor(Type)), show.legend =T, 
             position = pos) +
  geom_errorbar(aes(ymin = as.numeric(diab$LCI), 
                    ymax= as.numeric(diab$UCI)), position = pos,
                colour  ="black", width = 0.05) +
  ylab("Standardised Beta/Log Odds \n [95% Confidence Interval]")+ xlab ("")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(size = 12))+
  guides(colour = guide_legend(reverse = T)) +
  labs(shape = "Type")



diab_graph + ggtitle("Diabetes and Associated Phenotype") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title.align = 0.5) + coord_flip()


##Chronic Kidney Disease

ckd_con = basic_con[which(basic_con$Variable %in% "Creatinine"),] 
ckd_cat = basic_cat[which(basic_cat$Variable %in% "Chronic Kidney Disease"),]
ckd_con$Type = "Phenotype"
ckd_cat$Type = "Disease"
ckd = rbind(ckd_con, ckd_cat)

ckd$Clock <- factor(ckd$Clock, levels = rev(c("AgeAccelGrim", "Dunedin_PoAm", "AgeAccelPheno", "EEAA", "IEAA", "DNAmTLAdjAge"))) 
ckd$Variable = factor(ckd$Variable, levels=c("Creatinine", "Chronic Kidney Disease"))
ckd_graph <- ggplot(ckd,aes(y=as.numeric(ckd$test), 
                            x=Clock, group = Variable)) +   
  geom_point(data=ckd, size = 3, 
             aes(colour = Variable, shape = factor(Type)), show.legend =T, 
             position = pos) +
  geom_errorbar(aes(ymin = as.numeric(ckd$LCI), 
                    ymax= as.numeric(ckd$UCI)), position = pos,
                colour  ="black", width = 0.05) +
  ylab("Standardised Beta/Log Odds \n [95% Confidence Interval]")+ xlab ("")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(size = 12))+
  guides(colour = guide_legend(reverse = T)) +
  labs(shape = "Type")



ckd_graph + ggtitle("Chronic Kidney Disease and Associated Phenotype") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title.align = 0.5) + coord_flip()

## Pain Disorders 

##Back Pain
back = basic_cat[which(basic_cat$Variable %in% "Back Pain"),]
back$Type = "Disease"


back$Clock <- factor(back$Clock, levels = rev(c("AgeAccelGrim", "Dunedin_PoAm", "AgeAccelPheno", "EEAA", "IEAA", "DNAmTLAdjAge"))) 
back_graph <- ggplot(back,aes(y=as.numeric(back$test), 
                              x=Clock, group = Variable)) +   
  geom_point(data=back, size = 3, 
             aes(colour = Variable), show.legend =T, 
             position = pos) +
  geom_errorbar(aes(ymin = as.numeric(back$LCI), 
                    ymax= as.numeric(back$UCI)), position = pos,
                colour  ="black", width = 0.05) +
  ylab("Standardised Beta/Log Odds \n [95% Confidence Interval]")+ xlab ("")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(size = 12))+
  guides(colour = guide_legend(reverse = T)) 


back_graph + ggtitle("Back Pain") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title.align = 0.5) + coord_flip()



##Neck Pain 
neck = basic_cat[which(basic_cat$Variable %in% "Neck Pain"),]
neck$Type = "Disease"


neck$Clock <- factor(neck$Clock, levels = rev(c("AgeAccelGrim", "Dunedin_PoAm", "AgeAccelPheno", "EEAA", "IEAA", "DNAmTLAdjAge"))) 
neck_graph <- ggplot(neck,aes(y=as.numeric(neck$test), 
                              x=Clock, group = Variable)) +   
  geom_point(data=neck, size = 3, 
             aes(colour = Variable), show.legend =T, 
             position = pos) +
  geom_errorbar(aes(ymin = as.numeric(neck$LCI), 
                    ymax= as.numeric(neck$UCI)), position = pos,
                colour  ="black", width = 0.05) +
  ylab("Standardised Beta/Log Odds \n [95% Confidence Interval]")+ xlab ("")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(size = 12))+
  guides(colour = guide_legend(reverse = T)) 



neck_graph + ggtitle("Neck Pain") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title.align = 0.5) + coord_flip()

