## Preparation of BayesR+ Summary Statistic Files 

setwd("/Beta/")
cpgs = readRDS("../CpG_Names_BayesR.rds")
loop = list.files(".", ".csv")
for(i in loop){ 
tmp = fread(i)
tmp = as.data.frame(tmp)
print(i)
mean = apply(tmp, 2, mean)
mean_LCI = apply(tmp, 2, function(x) quantile(x, probs =  0.025))
mean_HCI = apply(tmp, 2, function(x) quantile(x, probs =  0.975))
names(mean) <- cpgs
names(mean_LCI) <- cpgs
names(mean_HCI) <- cpgs
mean = as.matrix(mean)
mean_LCI = as.matrix(mean_LCI)
mean_HCI = as.matrix(mean_HCI)
betas <- cbind(mean,cbind(mean_LCI,mean_HCI))
betas <- as.data.frame(betas)
betas$CpG <- row.names(betas)
betas <- betas[,c(4,1,2,3)]
names(betas)[2:4] <- c("Beta", "Beta_LCI", "Beta_HCI")
comp = fread(paste0("../Comp/", i))
comp = as.data.frame(comp)
pip <- lapply(comp,function(x){length(x[which(x>0)])/length(x)})
pip <- as.data.frame(reshape2::melt(unlist(pip)))
pip <- setDT(pip, keep.rownames = TRUE) 
names(pip) <- c("Marker", "PIP") 
pip$Marker <- cpgs
betas1 <- cbind(betas, pip[,2])
A <- gsub(".csv.*", "", i)
A <- gsub("\\.", "_", A)
print(A)
write.table(betas1, paste0("/BayesR_GWAS/", A, ".txt"), row.names = F, sep = ",", quote = F)
}
