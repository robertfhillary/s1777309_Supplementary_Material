setwd("/Neurology pQTL/") 

## Installing Requisite Packages 


if(!require(ggplot2)){
  install.packages("ggplot2")
}

if(!require(coloc)){
  install.packages("coloc")
}

if(!require(devtools)){
  install.packages("devtools")
}


if(!require(remotes)){
  install.packages("remotes")
}

remotes::install_github("MRCIEU/TwoSampleMR")

library(ggplot2)
library(coloc)

library(devtools)
install_github("phenoscanner/phenoscanner")
library(phenoscanner)

library(TwoSampleMR)


# Beta Coefficient vs Minor Allele Frequency 
x1 <- read.csv("pqtls.csv")
x1$Type[x1$Type == "CIS"] <- "Cis"
x1$Type[x1$Type == "TRANS"] <- "Trans"
beta = ggplot(x1, aes(x1$MAF, abs(x1$beta)))
beta1 = beta + geom_point(aes(colour = factor(x1$Type))) + xlab("Minor Allele Frequency") + ylab("Effect Size") +
  scale_y_continuous(breaks = c(0.5, 1.0, 1.5, 2.0, 2.5)) + theme(legend.title = element_blank()) + theme(legend.text = element_text(face = "italic")) + theme(legend.position=c(0.88, 0.84)) + theme(legend.key.size = unit(0.6, "cm")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                                                                                      panel.background = element_blank(), axis.line = element_line(colour = "black"))
beta1 + theme(axis.title.x = element_text(size=13.5), axis.title.y = element_text(size=13.5), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))



## Distance of Sentinel Variant from TSS vs Significance of Association 

x2 <- x1[x1$Sentinel == 1,]
x2$diff <- x2$diff/1e6
q <- ggplot(data = x2, aes(x2$diff, -log10(x2$p))) 
q1 = q + geom_point() + xlab("Distance of Sentinel Variant to TSS (Mb)") + labs(y = expression(-log[10](P))) + ylim(8.55,100) + xlim(-3, 3) +
  theme(legend.title = element_blank()) + theme(legend.text = element_text(face = "italic")) + theme(legend.position=c(0.88, 0.84)) + theme(legend.key.size = unit(0.6, "cm")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))                      

q1 + theme(axis.title.x = element_text(size=13.5), axis.title.y = element_text(size=13.5), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))



## Number of Variants per Protein 
y <- matrix(ncol = 2, nrow = 5)
colnames(y) <- c("Group", "Freq")
y <- data.frame(y)   
y$Group <- c("1", "2","3","4",">4")
y$Freq <- c(25, 7, 2, 0, 3)

cbPalette <- c("#9e66ab", "#f1595f", "#79c36a", "#599ad3", "#f9a65a")

bar = ggplot(data = y, aes(y$Group, y$Freq)) + geom_bar(stat="identity", fill = cbPalette) + xlab("Number of Loci") + ylab("Number of Proteins") 
bar + ylim(0,25) + theme(axis.title.x = element_text(size=13.5), axis.title.y = element_text(size=13.5), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))



## Genetic Atlas Figure 

x1 <- read.csv("pqtls.csv")

p = ggplot(x1, aes(x1$relative_chr_position_gwas_hit, x1$relative_chr_position_gene))
q = p + geom_jitter(aes(colour = as.factor(x1$Type))) + xlab("pQTL Position") + ylab("Protein Position")
q = q + scale_x_continuous(limits = c(1,23), breaks = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)) + scale_y_continuous(limits = c(1,23), breaks = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)) + theme_bw()
q + theme(legend.title = element_blank()) + theme(legend.text = element_text(face = "italic")) + theme(panel.background = element_rect(colour = "black")) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 0.1, color = "darkslategray4")



## Epigenetic Atlas Figure 

ewas3 <- read.csv("cpgs.csv")

p = ggplot(ewas3, aes(ewas3$Relative_Hit_Length, ewas3$Relative_Gene_Length))
q = p + geom_point(aes(colour = factor(ewas3$Type))) + xlab("CpG Position") + ylab("Protein Position")
q = q + scale_x_continuous(limits = c(1,23), breaks = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)) + scale_y_continuous(limits = c(1,23), breaks = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)) + theme_bw()
q = q + theme(legend.title = element_blank()) + theme(legend.text = element_text(face = "italic")) + theme(panel.background = element_rect(colour = "black")) + geom_vline(xintercept = 0)
q + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", size = 0.1)



## Comparison between Effect Sizes in Literature and Hillary et al 
betas <- read.csv("betas.csv")

beta = ggplot(betas, aes(betas$Literature, betas$Hillary))
beta1 = beta + geom_point(aes(colour = factor(betas$Protein), 
                              shape = factor(betas$Study)),size=2.35) + 
  scale_shape_manual(values=c(15,17,18))
beta2 = beta1 + 
  geom_smooth(method = "lm", se = FALSE) + 
  xlab("Betas from Literature") + 
  ylab(expression(paste("Betas from Hillary ", italic("et al.")))) +
  ggtitle("Correlation of Betas for pQTLs") +
  theme(plot.title= element_text(hjust =0.5)) + 
  labs(colour="Protein", shape="Study")

beta2 + geom_abline(intercept = 0, slope = 1) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)



## Bayesian Colocalisation Analysis Example - colocalisation of expression QTL and protein QTL data

  ## Example Protein: SIGLEC9

  ## Take 400kb region around sentinel cis variant 

  ## Extract known eQTLs in this region from phenoscanner 

res_eqtl <- phenoscanner(regionquery="chr19:51422766-51827766", proxies = "EUR", catalogue = "eQTL") 

  ## Fix Phenoscanner Results so no missing MAF, BETA present and eqtl results are in correct format

snps_from_phenoscanner <- res_eqtl$regions
results_from_phenoscanner <- res_eqtl$results

## Missing MAF 
if(length(which(results_from_phenoscanner$eur %in% "-")) >= 1) { 
  results_from_phenoscanner <- results_from_phenoscanner[-which(results_from_phenoscanner$eur %in% "-"),]
} 

## Missing BETA 
if(length(which(results_from_phenoscanner$beta %in% "NA")) >= 1) { 
  results_from_phenoscanner <- results_from_phenoscanner[-which(results_from_phenoscanner$beta %in% "NA"),]
} 

## Subset to Gene of Interest 
results_from_phenoscanner1 <- results_from_phenoscanner[which(results_from_phenoscanner$exp_gene %in% "SIGLEC9"),]
eqtls1 <- results_from_phenoscanner1

eqtls1$eur <- as.numeric(eqtls1$eur)
eqtls1$MAF <- ifelse(eqtls1$eur > 0.5, 1-eqtls1$eur, eqtls1$eur)
eqtls1$beta <- as.numeric(eqtls1$beta)

  ## Import GWAS Summary Statistics subset to 400 kb (+/- 200 kb) Region of Interest 

gwas_hits1 <- read.csv("siglec9_gwas.csv")

  
  ## Coloc analysis - calculation of abf (approximate Bayes Factor) and posterior probabilities 

dataset1 = list(snp = as.character(eqtls1$rsid), N = 31684, pvalues = as.numeric(eqtls1$p), MAF = as.numeric(eqtls1$MAF), type = "quant")
dataset2 = list(snp = as.character(gwas_hits1$marker), N = 750, pvalues = gwas_hits1$p, MAF = gwas_hits1$freq, type = "quant")
coloc = coloc.abf(dataset2, dataset1)


## Mendelian Randomisation analysis 
   
   ## swap exposure and outcome for bidirectional analysis 


## Open Exposure - replace with relevant datasets and assign column names appropriately 

exposure_dat <- read_exposure_data(
  filename = '/Proteomics/MR/ad_snps.txt',
  sep = ',',
  snp_col = 'SNP',
  beta_col = 'BETA',
  se_col = 'SE',
  effect_allele_col = 'A1',
  phenotype_col = '',
  units_col = '',
  other_allele_col = 'A2',
  eaf_col = '',
  samplesize_col = '',
  ncase_col = '',
  ncontrol_col = '',
  gene_col = '',
  pval_col = 'P'
)

## Clump Exposures

exposure_dat <- clump_data(exposure_dat)

## Open Outcome 

outcome_dat <- read_outcome_data(
  filename = '/Proteomics/MR/pvr1.txt',
  sep = ',',
  snp_col = 'marker',
  beta_col = 'beta',
  se_col = 'stderr',
  effect_allele_col = 'effect_allele',
  phenotype_col = '',
  units_col = '',
  other_allele_col = 'other_allele',
  eaf_col = 'freq',
  samplesize_col = 'n',
  ncase_col = '',
  ncontrol_col = '',
  gene_col = '',
  pval_col = 'p'
)

## Harmonisation of Alleles

dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)

## Run MR

mr_results <- mr(dat)

