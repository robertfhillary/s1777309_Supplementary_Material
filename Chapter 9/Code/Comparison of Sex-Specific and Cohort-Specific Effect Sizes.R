setwd("./Epigenetic Clocks WHO")

## Installing Requisite Packages 

if(!require(ggplot2)){
  install.packages("ggplot2")
}

library(ggplot2)

## Comparison of Effect Sizes 

## Discovery Cohort - Phenotypes 

males_basic <- read.csv("Effect_Sizes/Discovery_Basic_Phenotypes_Males.csv")
females_basic <- read.csv("Effect_Sizes/Discovery_Basic_Phenotypes_Females.csv")
which(males_basic$Clock %in% "AgeAccelGrim" & males_basic$Variable %in% "Pack Years")
males_basic = males_basic[-54,]
which(females_basic$Clock %in% "AgeAccelGrim" & females_basic$Variable %in% "Pack Years")
females_basic = females_basic[-54,]
cor(males_basic$Beta, females_basic$Beta)


names(males_basic)[names(males_basic) == "Beta"] <- "Betas - Males"
names(females_basic)[names(females_basic) == "Beta"] <- "Betas - Females"
females_basic = females_basic[,c(1,4)]
cor = cbind(males_basic, females_basic)
colnames(cor) <- make.unique(names(cor))


disc_phen = ggplot(cor, aes(cor$`Betas - Males`, cor$`Betas - Females`))
disc_phen = disc_phen+geom_point(aes(colour = factor(cor$Variable), 
                                     shape = factor(cor$Clock)))
disc_phen = disc_phen + xlab("Betas for Males") + 
  ylab("Betas for Females") +
  labs(colour="Trait", shape = "Clock") +
 theme(legend.title.align=0.5) +
 annotate("text", x = 0.49, y = 0.58, label = "", 
           size =4) + 
  annotate("text", x = 0.375, y = 0.3, label = expression(paste(italic("r"), "= 0.92")), 
           size =4.5) + ggtitle("Discovery Cohort - Basic Model - Continuous Variables") + 
  theme(plot.title = element_text(hjust =0.5)) 
disc_phen + geom_abline(intercept = 0, slope = 1, linetype = "dashed") + geom_hline(yintercept = 0, color = "dimgrey") + geom_vline(xintercept = 0, color = "dimgrey")



## Discovery Cohort - Diseases 

males_basic <- read.csv("Effect_Sizes/Discovery_Basic_Diseases_Males.csv")
females_basic <- read.csv("Effect_Sizes/Discovery_Basic_Diseases_Females.csv")
cor(males_basic$OR, females_basic$OR)

names(males_basic)[names(males_basic) == "OR"] <- "Odds Ratio - Males"
names(females_basic)[names(females_basic) == "OR"] <- "Odds Ratio - Females"
females_basic = females_basic[,c(1,5)]
cor = cbind(males_basic, females_basic)
colnames(cor) <- make.unique(names(cor))

cor = cor[-which(cor$Variable %in% c("Lung Cancer", "Bowel Cancer", "Breast Cancer")),] 
cor(cor$`Odds Ratio - Males`, cor$`Odds Ratio - Females`)

disc_phen = ggplot(cor, aes(cor$`Odds Ratio - Males`, cor$`Odds Ratio - Females`))
disc_phen = disc_phen+geom_point(aes(colour = factor(cor$Variable), 
                                     shape = factor(cor$Clock)))
disc_phen = disc_phen + xlab("Log Odds for Males") + 
  ylab("Log Odds for Females") +
  labs(colour="Trait", shape = "Clock") +
  theme(legend.title.align=0.5) +  
  annotate("text", x = 0.6, y = 0.3, label = expression(paste(italic("r"), "= 0.81")), 
           size =4.5) + ggtitle("Discovery Cohort - Basic Model - Disease States") + 
  theme(plot.title = element_text(hjust =0.5)) + 
  guides(shape = guide_legend(order = 1), colour = guide_legend(order = 2))
disc_phen + geom_abline(intercept = 0, slope = 1, linetype = "dashed") + geom_hline(yintercept = 0, color = "dimgrey") + geom_vline(xintercept = 0, color = "dimgrey")



 ## W/ Lung Cancer: 0.41
 ## W/o Lung Cancer: 0.73



## Replication Cohort - Phenotypes 

males_basic <- read.csv("Effect_Sizes/Replication_Basic_Phenotypes_Males.csv") 
females_basic <- read.csv("Effect_Sizes/Replication_Basic_Phenotypes_Females.csv")
which(males_basic$Clock %in% "AgeAccelGrim" & males_basic$Variable %in% "Pack Years")
males_basic = males_basic[-51,]
which(females_basic$Clock %in% "AgeAccelGrim" & females_basic$Variable %in% "Pack Years")
females_basic = females_basic[-51,]
cor(males_basic$Beta, females_basic$Beta)

  ## r = 0.88

names(males_basic)[names(males_basic) == "Beta"] <- "Betas - Males"
names(females_basic)[names(females_basic) == "Beta"] <- "Betas - Females"
females_basic = females_basic[,c(1,4)]
cor = cbind(males_basic, females_basic)
colnames(cor) <- make.unique(names(cor))


disc_phen = ggplot(cor, aes(cor$`Betas - Males`, cor$`Betas - Females`))
disc_phen = disc_phen+geom_point(aes(colour = factor(cor$Variable), 
                                     shape = factor(cor$Clock)))
disc_phen = disc_phen + xlab("Betas for Males") + 
  ylab("Betas for Females") +
  labs(colour="Trait", shape = "Clock") +
  theme(legend.title.align=0.5) +
  annotate("text", x = 0.45, y = 0.58, label = "", 
           size =4) + 
  annotate("text", x = 0.375, y = 0.3, label = expression(paste(italic("r"), "= 0.84")), 
           size =4.5) + ggtitle("Replication Cohort - Basic Model - Continuous Variables") + 
  theme(plot.title = element_text(hjust =0.5)) 
disc_phen + geom_abline(intercept = 0, slope = 1, linetype = "dashed") + geom_hline(yintercept = 0, color = "dimgrey") + geom_vline(xintercept = 0, color = "dimgrey")


## Replication Cohort - Diseases 

males_basic <- read.csv("Effect_Sizes/Replication_Basic_Diseases_Males.csv")
females_basic <- read.csv("Effect_Sizes/Replication_Basic_Diseases_Females.csv")
cor(males_basic$OR, females_basic$OR)

names(males_basic)[names(males_basic) == "OR"] <- "Odds Ratio - Males"
names(females_basic)[names(females_basic) == "OR"] <- "Odds Ratio - Females"
females_basic = females_basic[,c(1,5)]
cor = cbind(males_basic, females_basic)
colnames(cor) <- make.unique(names(cor))
cor = cor[-which(cor$Variable %in% c("Lung Cancer", "Bowel Cancer", "Breast Cancer")),] 
cor(cor$`Odds Ratio - Males`, cor$`Odds Ratio - Females`)

disc_phen = ggplot(cor, aes(cor$`Odds Ratio - Males`, cor$`Odds Ratio - Females`))
disc_phen = disc_phen+geom_point(aes(colour = factor(cor$Variable), 
                                     shape = factor(cor$Clock)))
disc_phen = disc_phen + xlab("Log Odds for Males") + 
  ylab("Log Odds for Females") +
  labs(colour="Trait", shape = "Clock") +
  theme(legend.title.align=0.5) +  
  annotate("text", x = 0.7, y = 0.3, label = expression(paste(italic("r"), "= 0.65")), 
           size =4.5) + ggtitle("Replication Cohort - Basic Model - Disease States") + 
  theme(plot.title = element_text(hjust =0.5)) +
  guides(shape = guide_legend(order = 1), colour = guide_legend(order = 2))
disc_phen + geom_abline(intercept = 0, slope = 1, linetype = "dashed") + geom_hline(yintercept = 0, color = "dimgrey") + geom_vline(xintercept = 0, color = "dimgrey")

cor(cor$`Odds Ratio - Males`, cor$`Odds Ratio - Females`)

## W/ Lung Cancer: 0.41
## W/o Lung Cancer: 0.57



## Discovery vs Replication - Phenotypes 


basic_males <- read.csv("Effect_Sizes/Discovery_Phenotypes.csv")
basic_females <- read.csv("Effect_Sizes/Replication_Phenotypes.csv")

 ## r = 0.96 

names(basic_males)[names(basic_males) == "Beta"] <- "Betas - Discovery"
names(basic_females)[names(basic_females) == "Beta"] <- "Betas - Replication"
basic_females = basic_females[,c(1,4)]
cor = cbind(basic_males, basic_females)
colnames(cor) <- make.unique(names(cor))
cor(cor$`Betas - Discovery`, cor$`Betas - Replication`)

disc_phen = ggplot(cor, aes(cor$`Betas - Discovery`, cor$`Betas - Replication`))
disc_phen = disc_phen+geom_point(aes(colour = factor(cor$Variable), 
                                     shape = factor(cor$Clock)))
disc_phen = disc_phen + xlab("Betas for Discovery Cohort") + 
  ylab("Betas for Replication Cohort") +
  labs(colour="Trait", shape = "Clock") +
  theme(legend.title.align=0.5) +  
  annotate("text", x = 0.38, y = 0.3, label = expression(paste(italic("r"), "= 0.95")), 
           size =4.5) +
  annotate("text", x = 0.4, y = 0.55, label = "", size = 4) + 
  ggtitle("Discovery vs Replication - Basic Model - Continuous Variables") + 
  theme(plot.title = element_text(hjust =0.5)) 
disc_phen + geom_abline(intercept = 0, slope = 1, linetype = "dashed") + geom_hline(yintercept = 0, color = "dimgrey") + geom_vline(xintercept = 0, color = "dimgrey")




## Discovery vs Replication - Diseases 


basic_males <- read.csv("Effect_Sizes/Discovery_Diseases.csv")
basic_females <- read.csv("Effect_Sizes/Replication_Diseases.csv")

   # w/lung cancer: 0.72 
   # w/o lung cancer: 0.81 

names(basic_males)[names(basic_males) == "OR"] <- "Odds Ratio - Discovery"
names(basic_females)[names(basic_females) == "OR"] <- "Odds Ratio - Replication"
basic_females = basic_females[,c(1,5)]
cor = cbind(basic_males, basic_females)
colnames(cor) <- make.unique(names(cor))
cor(cor$`Odds Ratio - Discovery`, cor$`Odds Ratio - Replication`, use = "complete.obs") 
cor = cor[-which(cor$Variable %in% c("Lung Cancer", "Bowel Cancer")),] 
cor$`Odds Ratio - Discovery` <- log(cor$`Odds Ratio - Discovery`)

disc_phen = ggplot(cor, aes(cor$`Odds Ratio - Discovery`, cor$`Odds Ratio - Replication`))
disc_phen = disc_phen+geom_point(aes(shape = factor(cor$Clock), 
                                     colour = factor(cor$Variable)))
disc_phen = disc_phen + xlab("Log Odds for Discovery Cohort") + 
  ylab("Log Odds for Replication Cohort") +
  labs(shape = "Clock", colour = "Trait") +
  theme(legend.title.align=0.5) +  
  annotate("text", x = 0.7, y = 0.43, label = expression(paste(italic("r"), "= 0.83")), 
           size =4.5) +
  annotate("text", x = 0.73, y = 1.17, label = "DNAm GrimAge - COPD", size = 4) + 
  ggtitle("Discovery vs Replication - Basic Model - Disease States") + 
  theme(plot.title = element_text(hjust =0.5)) + 
  guides(shape = guide_legend(order = 1), colour = guide_legend(order = 2))
disc_phen + geom_abline(intercept = 0, slope = 1, linetype = "dashed") + geom_hline(yintercept = 0, color = "dimgrey") + geom_vline(xintercept = 0, color = "dimgrey")

