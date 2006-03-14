source("../perfplot.R")
library("multcomp")

rdafiles <- list.files(pattern = "error")

JournalsRMSE <- c()
JournalsNPAR <- c()
model <- c()

for (f in rdafiles) {
    load(f)    
    JournalsRMSE <- cbind(JournalsRMSE, sqrt(error))
    JournalsNPAR <- cbind(JournalsNPAR, np)
    model <- c(model, strsplit(f, "_")[[1]][2])
}

if(identical(model, c("M5P", "ctree", "guide", "lm", "mob", "rpart"))) {
  model <- c("M5'", "CTree", "GUIDE", "lm", "MOB", "RPart")
  o <- c(5, 3, 1, 2, 6)
} else {
  o <- 1:length(model)
}

colnames(JournalsRMSE) <- model
colnames(JournalsNPAR) <- model
JournalsRMSE <- as.data.frame(JournalsRMSE[,o])
JournalsNPAR <- as.data.frame(JournalsNPAR[,o])

## median performance
sink(file = "Journals_results.txt")
print(round(cbind(MSE = sapply(JournalsRMSE, function(x) median(x^2)),
  Complexity = sapply(JournalsNPAR, median)), digits = 3))
sink()

### graphical comparison
perfplot(JournalsRMSE, file = "Journals_RMSE.pdf", lab = "RMSE")
perfplot(JournalsNPAR, file = "Journals_NPAR.pdf", lab = "Complexity")
