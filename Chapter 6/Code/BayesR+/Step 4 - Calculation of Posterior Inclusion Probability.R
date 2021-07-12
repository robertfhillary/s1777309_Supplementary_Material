
## Calculate Posterior Inclusion Probability of CpGs - over 95% inclusion means significant i.e.  PIP > 0.95 

setwd("Some_path/to/PIP") 
loop = list.files("../Comp/", pattern = ".csv") 
install.packages("data.table")
library(data.table) 
names = readRDS("../CpG_Names.rds") 
for(i in loop){ 
  comp <- fread(paste("../Comp/", i, sep=""))  
  comp<-as.data.frame(comp) 

  pip <- lapply(comp,function(x){length(x[which(x>0)])/length(x)})
  pip <- as.data.frame(reshape2::melt(unlist(pip)))
  pip <- setDT(pip, keep.rownames = TRUE) 
  names(pip) <- c("Marker", "PIP") 
  pip$Marker <- names
  A <-gsub(".csv*","",i)
  pip$Biomarker <- A 
  
write.csv(pip, file = paste(A, "_pip.csv", sep = ""), row.names = F) 
} 


## Find All CpGs that had PIP > 0.95

L <- list.files(".", ".csv")

pip_files <- lapply(L, read.csv, header = T) 
names <- as.character(L) 
names <- gsub("_.*", "", names) 
names(pip_files) <- names 
pip_files = Map(cbind, pip_files, "Biomarker"=names) 

pip_top <- do.call(rbind, lapply(pip_files, function(x)x[x$PIP > 0.95,]))
