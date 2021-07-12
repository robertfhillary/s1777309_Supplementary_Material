setwd("./GrimAge Cognitive/")

## Installing Requisite Packages

if(!require(corrplot)){
  install.packages("corrplot")
}


if(!require(ggplot2)){
  install.packages("ggplot2")
}


if(!require(foreign)){
  install.packages("foreign")
}

if(!require(psych)){
  install.packages("psych")
}

library(foreign)
library(psych)
library(corrplot)
library(ggplot2)


## Read in Phenotype File 

dataset = read.spss("LBC1936_GrimAge_RH_13MAY2019.sav", to.data.frame = TRUE)
head(dataset)
names(dataset)[1] <- "ID"
dataset$sex = NULL

## Read in GrimAge file

grim <- read.csv("Combined_Grim_LBC_Info.csv")
nrow(grim)

phenotype = merge(grim, dataset, by = "ID")
nrow(phenotype)

## COGNITIVE TRAITS 

## ApoE genotype - non significant

levels(as.factor(phenotype$APOEe4))
phenotype$APOEe4 <- as.character(phenotype$APOEe4)
phenotype$APOEe4[phenotype$APOEe4 == "No e4 allele"] <- 0
phenotype$APOEe4[phenotype$APOEe4 == "e4 allele"] <- 1
phenotype$APOEe4 <- as.factor(phenotype$APOEe4)

summary(glm(scale(phenotype$DNAmGrimAge) ~ phenotype$APOEe4))



## ApoE - e2/e3 vs e3/e3 

a = which(phenotype$APOEgenotype %in% c("e3/e3", "e2/e3"))
apoe = phenotype[a,]
apoe$APOEgenotype = as.character(apoe$APOEgenotype)
apoe$APOEgenotype[apoe$APOEgenotype %in% "e3/e3"] <- 0
apoe$APOEgenotype[apoe$APOEgenotype %in% "e2/e3"] <- 1
apoe$APOEgenotype = as.factor(apoe$APOEgenotype)
summary(glm(scale(apoe$AgeAccelGrim) ~ apoe$APOEgenotype))
plot(apoe$APOEgenotype, apoe$AgeAccelGrim)

## ApoE - e2/e2 vs e2/e3 

a = which(phenotype$APOEgenotype %in% c("e2/e2", "e2/e3"))
apoe = phenotype[a,]
apoe$APOEgenotype = as.character(apoe$APOEgenotype)
apoe$APOEgenotype[apoe$APOEgenotype %in% "e2/e2"] <- 0
apoe$APOEgenotype[apoe$APOEgenotype %in% "e2/e3"] <- 1
apoe$APOEgenotype = as.factor(apoe$APOEgenotype)
summary(glm(scale(apoe$AgeAccelGrim) ~ apoe$APOEgenotype))
plot(apoe$APOEgenotype, apoe$AgeAccelGrim)

## ApoE - e3/e3 vs e3/e4 

a = which(phenotype$APOEgenotype %in% c("e3/e3", "e3/e4"))
apoe = phenotype[a,]
apoe$APOEgenotype = as.character(apoe$APOEgenotype)
apoe$APOEgenotype[apoe$APOEgenotype %in% "e3/e3"] <- 0
apoe$APOEgenotype[apoe$APOEgenotype %in% "e3/e4"] <- 1
apoe$APOEgenotype = as.factor(apoe$APOEgenotype)
summary(glm(scale(apoe$AgeAccelGrim) ~ apoe$APOEgenotype))
plot(apoe$APOEgenotype, apoe$AgeAccelGrim)

## ApoE - e3/e4 vs e4/e4 

a = which(phenotype$APOEgenotype %in% c("e3/e4", "e4/e4"))
apoe = phenotype[a,]
apoe$APOEgenotype = as.character(apoe$APOEgenotype)
apoe$APOEgenotype[apoe$APOEgenotype %in% "e3/e4"] <- 0
apoe$APOEgenotype[apoe$APOEgenotype %in% "e4/e4"] <- 1
apoe$APOEgenotype = as.factor(apoe$APOEgenotype)
summary(glm(scale(apoe$AgeAccelGrim) ~ apoe$APOEgenotype))
plot(apoe$APOEgenotype, apoe$AgeAccelGrim)

## ApoE vs cognitive scores

apoe4_cognitive <- list()
for(i in colnames(phenotype)[42:73]){ 
  apoe4_cognitive[[i]] <- summary(glm(scale(phenotype[i]) ~ phenotype$APOEe4 + scale(phenotype$age) + as.factor(phenotype$sex)))
  }
coefs<-lapply(apoe4_cognitive,function(x)coef(x)[2,c(1:4)])
coefficient_df <- as.data.frame(coefs) 
coefficient_df <- t(coefficient_df)
coefficient_df <- as.data.frame(coefficient_df)
coefficient_df <- coefficient_df[order(coefficient_df$`Pr(>|t|)`),]


## General Factor of Mean Diffusivity 

phenotype$g_mean_diffusivity = as.numeric(principal(phenotype[,c(147, 149, 151, 153, 155, 157, 159, 161, 163, 165, 167, 169, 171, 173)], rotate = "none", residuals = FALSE, nfactors = 1)$scores)
summary(lm(scale(phenotype$g_mean_diffusivity) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$age) + factor(phenotype$sex))) 
summary(lm(scale(phenotype$g_mean_diffusivity) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$age) + factor(phenotype$sex) + scale(phenotype$age11IQ))) 

## General Factor of Fractional Anistropy 


phenotype$g_fractional_anisotropy = as.numeric(principal(phenotype[,c(146, 148, 150, 152, 154, 156, 158, 160, 162, 164, 166, 168, 170, 172)], rotate = "none", residuals = FALSE, nfactors = 1)$scores)
summary(lm(scale(phenotype$g_fractional_anisotropy) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$age) + factor(phenotype$sex))) 
summary(lm(scale(phenotype$g_fractional_anisotropy) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$age) + factor(phenotype$sex) + scale(phenotype$age11IQ))) 


## Brain_ICV 

summary(lm(scale(phenotype$brainIcv_ratio_w2) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$age) + factor(phenotype$sex)))
summary(lm(scale(phenotype$brainIcv_ratio_w2) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$age) + factor(phenotype$sex) + scale(phenotype$age11IQ)))



## WMH_ICV 

summary(lm(scale(phenotype$wmhIcv_ratio_w2) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$age) + factor(phenotype$sex)))
summary(lm(scale(phenotype$wmhIcv_ratio_w2) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$age) + factor(phenotype$sex) + scale(phenotype$age11IQ)))



## GM_ICV 

summary(lm(scale(phenotype$gmIcv_ratio_w2) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$age) + factor(phenotype$sex)))
summary(lm(scale(phenotype$gmIcv_ratio_w2) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$age) + factor(phenotype$sex) + scale(phenotype$age11IQ)))


## Nawm_ICV 

summary(lm(scale(phenotype$wmIcv_ratio_w2) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$age) + factor(phenotype$sex)))
summary(lm(scale(phenotype$wmIcv_ratio_w2) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$age) + factor(phenotype$sex) + scale(phenotype$age11IQ)))


## Derivation of general cognitive factor 


phenotype$g_factor = as.numeric(principal(phenotype[,c("digback_w2", "lnseq_w2", "matreas_w2", "blkdes_w2", "digsym_w2", "symsear_w2")], rotate = "none", residuals = FALSE, nfactors = 1)$scores)
summary(lm(scale(phenotype$g_factor) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$age) + factor(phenotype$sex))) 
summary(lm(scale(phenotype$g_factor) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$age) + factor(phenotype$sex) + scale(phenotype$age11IQ))) 


## DNAm GrimAge versus Cognitive Traits - without Age 11 IQ

grim_cognitive <- list()

for(i in colnames(phenotype)[c(42:55,176)]){ 
    grim_cognitive[[i]] <- summary(lm(scale(phenotype[i]) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$age) + as.factor(phenotype$sex)))
}
coefs<-lapply(grim_cognitive,function(x)coef(x)[2,c(1:4)])
coefficient_df <- as.data.frame(coefs) 
coefficient_df <- t(coefficient_df)
coefficient_df <- as.data.frame(coefficient_df)
coefficient_df <- coefficient_df[order(coefficient_df$`Pr(>|t|)`),]

grim_v_cognitive <- coefficient_df
names(grim_v_cognitive)
grim_v_cognitive$'97.5% CI' <- grim_v_cognitive$Estimate + 1.96*grim_v_cognitive$`Std. Error`
grim_v_cognitive$'2.5% CI' <- grim_v_cognitive$Estimate - 1.96*grim_v_cognitive$`Std. Error`
names(grim_v_cognitive)[1] <- "Standarised Beta"
names(grim_v_cognitive)[2] <- "Standard Error"
names(grim_v_cognitive)[3] <- "t Statistic"
names(grim_v_cognitive)[4] <- "P Value"
grim_v_cognitive$Trait <- rownames(grim_v_cognitive)
grim_v_cognitive <- grim_v_cognitive[,c(7,1,6,5, 2:4)]
grim_v_cognitive$Trait <- gsub("_w2", "", grim_v_cognitive$Trait)

View(grim_v_cognitive)
write.csv(grim_v_cognitive, "Regressions/DNAmGrimAge_Cognitive.csv", row.names = F)

## DNAm GrimAge versus Cognitive Traits - with Age 11 IQ

grim_cognitive_11 <- list()
for(i in colnames(phenotype)[c(42:55,176)]){ 
  grim_cognitive_11[[i]] <- summary(lm(scale(phenotype[i]) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$age) + as.factor(phenotype$sex) + scale(phenotype$age11IQ)))
}
coefs<-lapply(grim_cognitive_11,function(x)coef(x)[2,c(1:4)])
coefficient_df <- as.data.frame(coefs) 
coefficient_df <- t(coefficient_df)
coefficient_df <- as.data.frame(coefficient_df)
coefficient_df <- coefficient_df[order(coefficient_df$`Pr(>|t|)`),]

grim_v_cognitive_11 <- coefficient_df
names(grim_v_cognitive_11)
grim_v_cognitive_11$'97.5% CI' <- grim_v_cognitive_11$Estimate + 1.96*grim_v_cognitive_11$`Std. Error`
grim_v_cognitive_11$'2.5% CI' <- grim_v_cognitive_11$Estimate - 1.96*grim_v_cognitive_11$`Std. Error`
names(grim_v_cognitive_11)[1] <- "Standarised Beta"
names(grim_v_cognitive_11)[2] <- "Standard Error"
names(grim_v_cognitive_11)[3] <- "t Statistic"
names(grim_v_cognitive_11)[4] <- "P Value"
grim_v_cognitive_11$Trait <- rownames(grim_v_cognitive_11)
grim_v_cognitive_11 <- grim_v_cognitive_11[,c(7,1,6,5, 2:4)]
grim_v_cognitive_11$Trait <- gsub("_w2", "", grim_v_cognitive_11$Trait)

View(grim_v_cognitive_11)
write.csv(grim_v_cognitive_11, "Regressions/GrimAge_Cognitive_Age11IQ.csv", row.names = F)


## BLOOD TRAITS 


## Folate

phenotype$bld_folate_w2 <- as.character(phenotype$bld_folate_w2)
phenotype$bld_folate_w2[phenotype$bld_folate_w2 %in% ">20 ug/l"] <- NA
phenotype$bld_folate_w2 <- as.numeric(phenotype$bld_folate_w2)


## B12

phenotype$bld_b12_w2 <- as.character(phenotype$bld_b12_w2)
phenotype$bld_b12_w2[phenotype$bld_b12_w2 %in% ">1000"] <- NA
phenotype$bld_b12_w2 <- as.numeric(phenotype$bld_b12_w2)



## CRProt

phenotype$bld_crprot_w2 <- as.character(phenotype$bld_crprot_w2)
phenotype$bld_crprot_w2[phenotype$bld_crprot_w2 %in% "<3 mg/l"] <- 0
phenotype$bld_crprot_w2 <- as.numeric(phenotype$bld_crprot_w2)



## TSH

phenotype$bld_tsh_w2 <- as.character(phenotype$bld_tsh_w2)
phenotype$bld_tsh_w2[phenotype$bld_tsh_w2 %in% "<0.001"] <- 0
phenotype$bld_tsh_w2 <- as.numeric(phenotype$bld_tsh_w2)


## S1000

phenotype$bld_S100_w2 <- as.character(phenotype$bld_S100_w2)
phenotype$bld_S100_w2[phenotype$bld_S100_w2 %in% "<0.02"] <- 0
phenotype$bld_S100_w2 <- as.numeric(phenotype$bld_S100_w2)

##hsCRP

phenotype$bld_hsCRP_w2 <- as.character(phenotype$bld_hsCRP_w2)
phenotype$bld_hsCRP_w2[phenotype$bld_hsCRP_w2 %in% "<0.154"] <- 0
phenotype$bld_hsCRP_w2 <- as.numeric(phenotype$bld_hsCRP_w2)



## DNAm GrimAge versus Blood Traits - without Age 11 IQ

grim_v_blood <- list()

for(i in colnames(phenotype)[c(81:103, 105:126)]){ 
  grim_v_blood[[i]] <- summary(lm(scale(phenotype[i]) ~ scale(phenotype$AgeAccelGrim) + + scale(phenotype$age) + as.factor(phenotype$sex)))
}
coefs<-lapply(grim_v_blood,function(x)coef(x)[2,c(1:4)])
coefficient_df <- as.data.frame(coefs) 
coefficient_df <- t(coefficient_df)
coefficient_df <- as.data.frame(coefficient_df)
coefficient_df <- coefficient_df[order(coefficient_df$`Pr(>|t|)`),]
coefficient_df$FDR <- p.adjust(coefficient_df$`Pr(>|t|)`, method = "BH")

grim_v_blood <- coefficient_df
names(grim_v_blood)
grim_v_blood$'97.5% CI' <- grim_v_blood$Estimate + 1.96*grim_v_blood$`Std. Error`
grim_v_blood$'2.5% CI' <- grim_v_blood$Estimate - 1.96*grim_v_blood$`Std. Error`
names(grim_v_blood)[1] <- "Standarised Beta"
names(grim_v_blood)[2] <- "Standard Error"
names(grim_v_blood)[3] <- "t Statistic"
names(grim_v_blood)[4] <- "P Value"
names(grim_v_blood)[5] <- "FDR-adjusted P Value"
grim_v_blood$Trait <- rownames(grim_v_blood)
grim_v_blood <- grim_v_blood[,c(8,1,7,6, 2:5)]
grim_v_blood$Trait <- gsub("_w2", "", grim_v_blood$Trait)
grim_v_blood$Trait <- gsub("bld_", "", grim_v_blood$Trait)


View(grim_v_blood)
write.csv(grim_v_blood, "Regressions/GrimAge_Blood.csv", row.names = F)



## DNAm GrimAge versus Blood Traits - with Age 11 IQ

grim_blood_11 <- list()

for(i in colnames(phenotype)[c(81:103, 105:126)]){ 
  grim_blood_11[[i]] <- summary(lm(scale(phenotype[i]) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$age) + as.factor(phenotype$sex) + scale(phenotype$age11IQ)))
}
coefs<-lapply(grim_blood_11,function(x)coef(x)[2,c(1:4)])
coefficient_df <- as.data.frame(coefs) 
coefficient_df <- t(coefficient_df)
coefficient_df <- as.data.frame(coefficient_df)
coefficient_df <- coefficient_df[order(coefficient_df$`Pr(>|t|)`),]
coefficient_df$FDR <- p.adjust(coefficient_df$`Pr(>|t|)`, method = "BH")

grim_v_blood_11 <- coefficient_df
names(grim_v_blood_11)
grim_v_blood_11$'97.5% CI' <- grim_v_blood_11$Estimate + 1.96*grim_v_blood_11$`Std. Error`
grim_v_blood_11$'2.5% CI' <- grim_v_blood_11$Estimate - 1.96*grim_v_blood_11$`Std. Error`
names(grim_v_blood_11)[1] <- "Standarised Beta"
names(grim_v_blood_11)[2] <- "Standard Error"
names(grim_v_blood_11)[3] <- "t Statistic"
names(grim_v_blood_11)[4] <- "P Value"
names(grim_v_blood_11)[5] <- "FDR-adjusted P Value"
grim_v_blood_11$Trait <- rownames(grim_v_blood_11)
grim_v_blood_11 <- grim_v_blood_11[,c(8,1,7,6, 2:5)]
grim_v_blood_11$Trait <- gsub("_w2", "", grim_v_blood_11$Trait)
grim_v_blood_11$Trait <- gsub("bld_", "", grim_v_blood_11$Trait)


View(grim_v_blood_11)


## PHYSICAL TRAITS 


phenotype$bld_eGFR_w2 <- as.character(phenotype$bld_eGFR_w2)
phenotype$bld_eGFR_w2[phenotype$bld_eGFR_w2 %in% ">60 ml/min"] <- NA
phenotype$bld_eGFR_w2 <- as.numeric(phenotype$bld_eGFR_w2)


## DNAm GrimAge versus Physical Traits - without Age 11 IQ


physical_continuous <- list()
for(i in colnames(phenotype[c(78, 79, 107, 117, 128:132)])){
  physical_continuous[[i]] <- summary(lm(scale(phenotype[,i]) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$agedays_w2) + as.factor(phenotype$sex)))
}
coefs<-lapply(physical_continuous,function(x)coef(x)[2,c(1:4)])
coefficient_df <- as.data.frame(coefs) 
coefficient_df <- t(coefficient_df)
coefficient_df <- as.data.frame(coefficient_df)
coefficient_df <- coefficient_df[order(coefficient_df$`Pr(>|t|)`),]
coefficient_df$Trait <- rownames(coefficient_df)

fev = summary(lm(scale(phenotype$fev_w2) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$agedays_w2) + as.factor(phenotype$sex) + scale(phenotype$height_w2) + scale(phenotype$DNAmPACKYRS)))
fev = as.data.frame(t(fev$coefficients[2,c(1:4)]))
fev$Trait <- "FEV"

fer = summary(lm(scale(phenotype$fer_w2) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$agedays_w2) + as.factor(phenotype$sex) + scale(phenotype$height_w2) + scale(phenotype$DNAmPACKYRS)))
fer = as.data.frame(t(fer$coefficients[2,c(1:4)]))
fer$Trait <- "FER"

fvc = summary(lm(scale(phenotype$fvc_w2) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$agedays_w2) + as.factor(phenotype$sex) + scale(phenotype$height_w2) + scale(phenotype$DNAmPACKYRS)))
fvc = as.data.frame(t(fvc$coefficients[2,c(1:4)]))
fvc$Trait <- "FVC"

pef = summary(lm(scale(phenotype$pef_w2) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$agedays_w2) + as.factor(phenotype$sex) + scale(phenotype$height_w2) + scale(phenotype$DNAmPACKYRS)))
pef = as.data.frame(t(pef$coefficients[2,c(1:4)]))
pef$Trait <- "PEF"

physical_coefficients <- rbind(fev, fer)
physical_coefficients <- rbind(physical_coefficients, fvc)
physical_coefficients <- rbind(physical_coefficients, pef)
names(physical_coefficients)[4] <- "Pr(>|t|)"

physical_traits = rbind(coefficient_df, physical_coefficients)
physical_traits$FDR <- p.adjust(physical_traits$`Pr(>|t|)`, method = "BH")
physical_traits <- physical_traits[order(physical_traits$FDR),] 

physical_traits$'97.5% CI' <- physical_traits$Estimate + 1.96*physical_traits$`Std. Error`
physical_traits$'2.5% CI' <- physical_traits$Estimate - 1.96*physical_traits$`Std. Error`
names(physical_traits)[1] <- "Standarised Beta"
names(physical_traits)[2] <- "Standard Error"
names(physical_traits)[3] <- "t Statistic"
names(physical_traits)[4] <- "P Value"
names(physical_traits)[5] <- "Trait"
names(physical_traits)[6] <- "FDR-adjusted P Value"
physical_traits <- physical_traits[,c(5,1,8,7, 2:4, 6)]
physical_traits$Trait <- gsub("_w2", "", physical_traits$Trait)
physical_traits$Trait <- gsub("_w1", "", physical_traits$Trait)

View(physical_traits)
write.csv(physical_traits, "Regressions/DNAmGrimAge_Physical.csv", row.names = F)


## DNAm GrimAge versus Physical Traits - with Age 11 IQ


physical_continuous <- list()
for(i in colnames(phenotype[c(78, 79, 104, 128:130, 132)])){
  physical_continuous[[i]] <- summary(lm(scale(phenotype[,i]) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$agedays_w2) + as.factor(phenotype$sex) + scale(phenotype$age11IQ)))
}
coefs<-lapply(physical_continuous,function(x)coef(x)[2,c(1:4)])
coefficient_df <- as.data.frame(coefs) 
coefficient_df <- t(coefficient_df)
coefficient_df <- as.data.frame(coefficient_df)
coefficient_df <- coefficient_df[order(coefficient_df$`Pr(>|t|)`),]
coefficient_df$Trait <- rownames(coefficient_df)

fev = summary(lm(scale(phenotype$fev_w2) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$agedays_w2) + as.factor(phenotype$sex) + scale(phenotype$height_w2)  + scale(phenotype$age11IQ)))
fev = as.data.frame(t(fev$coefficients[2,c(1:4)]))
fev$Trait <- "FEV"

fer = summary(lm(scale(phenotype$fer_w2) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$agedays_w2) + as.factor(phenotype$sex) + scale(phenotype$height_w2) + scale(phenotype$age11IQ)))
fer = as.data.frame(t(fer$coefficients[2,c(1:4)]))
fer$Trait <- "FER"

fvc = summary(lm(scale(phenotype$fvc_w2) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$agedays_w2) + as.factor(phenotype$sex) + scale(phenotype$height_w2)  + scale(phenotype$age11IQ)))
fvc = as.data.frame(t(fvc$coefficients[2,c(1:4)]))
fvc$Trait <- "FVC"

pef = summary(lm(scale(phenotype$pef_w2) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$agedays_w2) + as.factor(phenotype$sex) + scale(phenotype$height_w2) + scale(phenotype$age11IQ)))
pef = as.data.frame(t(pef$coefficients[2,c(1:4)]))
pef$Trait <- "PEF"

physical_coefficients <- rbind(fev, fer)
physical_coefficients <- rbind(physical_coefficients, fvc)
physical_coefficients <- rbind(physical_coefficients, pef)
names(physical_coefficients)[4] <- "Pr(>|t|)"

physical_traits = rbind(coefficient_df, physical_coefficients)
physical_traits$FDR <- p.adjust(physical_traits$`Pr(>|t|)`, method = "BH")
physical_traits <- physical_traits[order(physical_traits$FDR),] 

physical_traits$'97.5% CI' <- physical_traits$Estimate + 1.96*physical_traits$`Std. Error`
physical_traits$'2.5% CI' <- physical_traits$Estimate - 1.96*physical_traits$`Std. Error`
names(physical_traits)[1] <- "Standarised Beta"
names(physical_traits)[2] <- "Standard Error"
names(physical_traits)[3] <- "t Statistic"
names(physical_traits)[4] <- "P Value"
names(physical_traits)[5] <- "Trait"
names(physical_traits)[6] <- "FDR-adjusted P Value"
physical_traits <- physical_traits[,c(5,1,8,7, 2:4, 6)]
physical_traits$Trait <- gsub("_w2", "", physical_traits$Trait)
physical_traits$Trait <- gsub("_w1", "", physical_traits$Trait)

physical_traits_11 <- physical_traits
View(physical_traits_11)
write.csv(physical_traits_11, "Regressions/GrimAge_Physical_age11.csv", row.names = F)


## PLOTS 

## Plot All Significant Results 

plot <- read.csv("all_significant_1.csv")

## Order traits by magnitude of coefficient - need to specify unique(trait) as traits are entered twice 
plot$Trait <- factor(plot$Trait, levels=unique(as.character(plot$Trait)[rev(order(plot$Standarised_Beta))]))

## ggplot step 
all_significant_graph <- ggplot(plot,aes(y=plot$Standarised_Beta, x=plot$Trait, colour = plot$Type, shape = plot$Age11IQ)) +
  geom_point(position=position_dodge(width=0.5), size = 3) +
  geom_errorbar(aes(ymin = plot$LCI, ymax = plot$UCI),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("Standardised Beta")+ xlab ("")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 10, vjust = 0.5),
        plot.title = element_text(size = 12))+
  scale_y_continuous(limit = c(-0.402, 0.6))+ guides(colour = F) + labs(shape = "Age 11 IQ") +
  coord_flip()

# split by group
all_significant_graph + facet_wrap(vars(Type), scales = "free_y")

## Plot of Grim Age versus Cognitive Traits - not adjusted for age 11 IQ 

grim_cognitive <- list()

for(i in colnames(phenotype)[42:73]){ 
  grim_cognitive[[i]] <- summary(lm(scale(phenotype[i]) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$age) + as.factor(phenotype$sex)))
}
coefs<-lapply(grim_cognitive,function(x)coef(x)[2,c(1:4)])
coefficient_df <- as.data.frame(coefs) 
coefficient_df <- t(coefficient_df)
coefficient_df <- as.data.frame(coefficient_df)
coefficient_df <- coefficient_df[order(coefficient_df$`Pr(>|t|)`),]
coefficient_df$FDR <- p.adjust(coefficient_df$`Pr(>|t|)`, method = "BH")

grim_v_cognitive <- coefficient_df
names(grim_v_cognitive)

a = which(grim_v_cognitive$FDR < 0.05)
cognitive_plot <- grim_v_cognitive[a,]
cognitive_plot$variable <- rownames(cognitive_plot)
cognitive_plot$upper_95 <- cognitive_plot$Estimate + 1.96*cognitive_plot$`Std. Error`
cognitive_plot$lower_95 <- cognitive_plot$Estimate - 1.96*cognitive_plot$`Std. Error`

cognitive_plot$variable <- factor(cognitive_plot$variable, levels=cognitive_plot$variable[order(cognitive_plot$Estimate, decreasing=F)])

ggplot(cognitive_plot,aes(y=cognitive_plot$Estimate, x=cognitive_plot$variable)) + 
  geom_point(position=position_dodge(width=0.5), size = 3, colour = "#7CAE00")+
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("standardised beta")+ xlab ("")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), legend.position = "none")+
  coord_flip()


## DNAm GrimAge versus Cognitive Traits - adjsuted for age 11 IQ

grim_cognitive_11 <- list()
for(i in colnames(phenotype)[42:73]){ 
  grim_cognitive_11[[i]] <- summary(lm(scale(phenotype[i]) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$age) + as.factor(phenotype$sex) + scale(phenotype$age11IQ)))
}
coefs<-lapply(grim_cognitive_11,function(x)coef(x)[2,c(1:4)])
coefficient_df <- as.data.frame(coefs) 
coefficient_df <- t(coefficient_df)
coefficient_df <- as.data.frame(coefficient_df)
coefficient_df <- coefficient_df[order(coefficient_df$`Pr(>|t|)`),]
coefficient_df$FDR <- p.adjust(coefficient_df$`Pr(>|t|)`, method = "BH")

grim_v_cognitive_11 <- coefficient_df
names(grim_v_cognitive_11)

a = which(grim_v_cognitive_11$FDR < 0.05)
cognitive_plot_11 <- grim_v_cognitive_11[a,]
cognitive_plot_11$variable <- rownames(cognitive_plot_11)
cognitive_plot_11$upper_95 <- cognitive_plot_11$Estimate + 1.96*cognitive_plot_11$`Std. Error`
cognitive_plot_11$lower_95 <- cognitive_plot_11$Estimate - 1.96*cognitive_plot_11$`Std. Error`

cognitive_plot_11$variable <- factor(cognitive_plot_11$variable, levels=cognitive_plot_11$variable[order(cognitive_plot_11$Estimate, decreasing=F)])

ggplot(cognitive_plot_11,aes(y=cognitive_plot_11$Estimate, x=cognitive_plot_11$variable)) + 
  geom_point(position=position_dodge(width=0.5), size = 3, colour = "#7CAE00")+
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("standardised beta")+ xlab ("")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), legend.position = "none")+
  coord_flip()

## Plot of Grim Age versus Blood Traits - not adjusted for age 11 IQ 


grim_v_blood <- list()

for(i in colnames(phenotype)[c(81:103, 105:126)]){ 
  grim_v_blood[[i]] <- summary(lm(scale(phenotype[i]) ~ scale(phenotype$AgeAccelGrim) + as.factor(phenotype$sex)))
}
coefs<-lapply(grim_v_blood,function(x)coef(x)[2,c(1:4)])
coefficient_df <- as.data.frame(coefs) 
coefficient_df <- t(coefficient_df)
coefficient_df <- as.data.frame(coefficient_df)
coefficient_df <- coefficient_df[order(coefficient_df$`Pr(>|t|)`),]
coefficient_df$FDR <- p.adjust(coefficient_df$`Pr(>|t|)`, method = "BH")

grim_v_blood <- coefficient_df
names(grim_v_blood)

a = which(grim_v_blood$FDR < 0.05)
blood_plot <- grim_v_blood[a,]
blood_plot$variable <- rownames(blood_plot)
blood_plot$upper_95 <- blood_plot$Estimate + 1.96*blood_plot$`Std. Error`
blood_plot$lower_95 <- blood_plot$Estimate - 1.96*blood_plot$`Std. Error`

blood_plot$variable <- factor(blood_plot$variable, levels=blood_plot$variable[order(blood_plot$Estimate, decreasing=F)])

ggplot(blood_plot,aes(y=blood_plot$Estimate, x=blood_plot$variable)) + 
  geom_point(position=position_dodge(width=0.5), size = 3, colour = "#009cae")+
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("standardised beta")+ xlab ("")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), legend.position = "none")+
  coord_flip()

## Plot of Grim Age versus Blood Traits - adjusted age 11 IQ 


grim_blood_11 <- list()

for(i in colnames(phenotype)[c(81:103, 105:126)]){ 
  grim_blood_11[[i]] <- summary(lm(scale(phenotype[i]) ~ scale(phenotype$AgeAccelGrim) + scale(phenotype$age) + as.factor(phenotype$sex) + scale(phenotype$age11IQ)))
}
coefs<-lapply(grim_blood_11,function(x)coef(x)[2,c(1:4)])
coefficient_df <- as.data.frame(coefs) 
coefficient_df <- t(coefficient_df)
coefficient_df <- as.data.frame(coefficient_df)
coefficient_df <- coefficient_df[order(coefficient_df$`Pr(>|t|)`),]
coefficient_df$FDR <- p.adjust(coefficient_df$`Pr(>|t|)`, method = "BH")

grim_v_blood_11 <- coefficient_df
names(grim_v_blood_11)

a = which(grim_v_blood_11$FDR < 0.05)
blood_plot_11 <- grim_v_blood_11[a,]
blood_plot_11$variable <- rownames(blood_plot_11)
blood_plot_11$upper_95 <- blood_plot_11$Estimate + 1.96*blood_plot_11$`Std. Error`
blood_plot_11$lower_95 <- blood_plot_11$Estimate - 1.96*blood_plot_11$`Std. Error`

blood_plot_11$variable <- factor(blood_plot_11$variable, levels=blood_plot_11$variable[order(blood_plot_11$Estimate, decreasing=F)])

ggplot(blood_plot_11,aes(y=blood_plot_11$Estimate, x=blood_plot_11$variable)) + 
  geom_point(position=position_dodge(width=0.5), size = 3, colour = "#009cae")+
  geom_errorbar(aes(ymin = lower_95, ymax = upper_95),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("standardised beta")+ xlab ("")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), legend.position = "none")+
  coord_flip()

## Plot of Grim Age versus Physical Traits - not adjusted for age 11 IQ 

physical_plot <- physical_traits[physical_traits$`FDR-adjusted P Value` < 0.05, ]


physical_plot$Trait <- factor(physical_plot$Trait, levels=physical_plot$Trait[order(physical_plot$`Standarised Beta`, decreasing=F)])


ggplot(physical_plot,aes(y=physical_plot$`Standarised Beta`, x=physical_plot$Trait)) + 
  geom_point(position=position_dodge(width=0.5), size = 3, colour = "#009cae")+
  geom_errorbar(aes(ymin = physical_plot$`2.5% CI`, ymax = physical_plot$`97.5% CI`),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("standardised beta")+ xlab ("")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), legend.position = "none")+
  coord_flip()


## Plot of Grim Age versus Physical Traits - adjusted for age 11 IQ 

physical_plot <- physical_traits_11[physical_traits_11$`FDR-adjusted P Value` < 0.05, ]


physical_plot$Trait <- factor(physical_plot$Trait, levels=physical_plot$Trait[order(physical_plot$`Standarised Beta`, decreasing=F)])


ggplot(physical_plot,aes(y=physical_plot$`Standarised Beta`, x=physical_plot$Trait)) + 
  geom_point(position=position_dodge(width=0.5), size = 3, colour = "#009cae")+
  geom_errorbar(aes(ymin = physical_plot$`2.5% CI`, ymax = physical_plot$`97.5% CI`),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("standardised beta")+ xlab ("")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), legend.position = "none")+
  coord_flip()


## Forest Plot of Hazard Ratios - Different Predictors  

forest = read.csv("Forest_Plots.csv")
head(forest)
forest$Shape <- ifelse(forest$Shape %in% "Before", "Unadjusted", "Adjusted")
forest$Trait <- factor(forest$Trait, levels=unique(as.character(forest$Trait)[order(forest$HR)]))
fp = ggplot(data=forest, aes(x = forest$Trait, y=forest$HR, shape = forest$Shape)) +  
  geom_point(position=position_dodge(width=0.5), size = 3) +
  geom_errorbar(aes(ymin = forest$Low, ymax = forest$High),position = position_dodge(0.5), width = 0.05, colour = "black")  +
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Trait") + ylab("Hazard Ratio (95% CI)") + ## Rename axes
  labs(shape="Age 11 IQ") + theme(legend.title.align=0.5) + ## Title and centre legend
  scale_fill_discrete(guide = guide_legend(reverse=TRUE)) ## Flip order of fill legend 

print(fp)



## Plot of Effect Size Comparison Between AgeAccelGrim and DNAmPACKYEARS 

effect <- read.csv("Effect Sizes Shared.csv", check.names = F)
effect <- effect[,c(1:6)]
cor(effect$`Standarised Beta`, effect$`Standardised Beta Pack`)
cor(effect$`Standarised Beta Adjusted`, effect$`Standarised Beta Adjusted Pack Years`, use = "complete.obs")

## Before Adjustment for Age 11 IQ  

beta = ggplot(effect, aes(effect$`Standarised Beta`, effect$`Standardised Beta Pack`))
beta1 = beta+geom_point()
beta2 = beta1 + xlab("Betas for DNAm GrimAge") + 
  ylab("Betas for DNAm PackYears") 
beta2 + geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_hline(yintercept = 0, color = "dimgrey") + 
  geom_vline(xintercept = 0, color = "dimgrey") + 
  ggtitle("Correlation Before Adjustment for Age 11 IQ") + theme(plot.title = element_text(hjust = 0.5))

## After Adjustment for Age 11 IQ  

beta = ggplot(effect, aes(effect$`Standarised Beta Adjusted`, effect$`Standarised Beta Adjusted Pack Years`))
beta1 = beta+geom_point()
beta2 = beta1 + xlab("Betas for DNAm GrimAge") + 
  ylab("Betas for DNAm PackYears") 
beta2 + geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_hline(yintercept = 0, color = "dimgrey") + 
  geom_vline(xintercept = 0, color = "dimgrey") + 
  ggtitle("Correlation After Adjustment for Age 11 IQ") + theme(plot.title = element_text(hjust = 0.5))



## Correlation Plot between Predictors

cor<-phenotype[,c(18, 10:17)]
cor <- cor(as.matrix(cor)) 
round(cor,2)
corrplot(cor, type = "upper", method = "circle", tl.col = "black", tl.cex = 0.8, addCoef.col = "black")



## Automation of Descriptive Statistics 

mean <- list()
 for(i in colnames(phenotype)[42:73]){ 
   mean[[i]] <- mean(phenotype[,i], na.rm = T)
   }
mean1 = as.data.frame(t(as.data.frame(mean)))
names(mean1)[1] <- "mean"

sd <- list()
for(i in colnames(phenotype)[42:73]){ 
  sd[[i]] <- sd(phenotype[,i], na.rm = T)
}
sd1 = as.data.frame(t(as.data.frame(sd)))
names(sd1)[1] <- "sd"

n <- list()
for(i in colnames(phenotype)[42:73]){ 
  n[[i]] <- length(which(!is.na(phenotype[,i])))
}
n1 = as.data.frame(t(as.data.frame(n)))
names(n1)[1] <- "n"

mean1$Trait <- rownames(mean1)
sd1$Trait <- rownames(sd1)
n1$Trait <- rownames(n1)
cog_sum <- merge(mean1, sd1, by = "Trait")
cog_sum <- merge(cog_sum, n1, by = "Trait")
View(cog_sum)
cog_sum$Trait <- gsub("_w2", "", cog_sum$Trait)
write.csv(cog_sum, "Cognitive_Descriptive.csv")
