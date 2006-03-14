source("../perfplot.R")
library("multcomp")

rdafiles <- list.files(pattern = "error")

BostonHousingRMSE <- c()
BostonHousingNPAR <- c()
model <- c()

for (f in rdafiles) {
    load(f)    
    BostonHousingRMSE <- cbind(BostonHousingRMSE, sqrt(error))
    BostonHousingNPAR <- cbind(BostonHousingNPAR, np)
    model <- c(model, strsplit(f, "_")[[1]][2])
}

if(identical(model, c("M5P", "ctree", "guide", "lm", "mob", "rpart"))) {
  model <- c("M5'", "CTree", "GUIDE", "lm", "MOB", "RPart")
  o <- c(5, 3, 1, 2, 6)
} else {
  o <- 1:length(model)
}

colnames(BostonHousingRMSE) <- model
colnames(BostonHousingNPAR) <- model
BostonHousingRMSE <- as.data.frame(BostonHousingRMSE[,o])
BostonHousingNPAR <- as.data.frame(BostonHousingNPAR[,o])

## median performance
sink(file = "BostonHousing_results.txt")
print(round(cbind(MSE = sapply(BostonHousingRMSE, function(x) median(x^2)),
  Complexity = sapply(BostonHousingNPAR, median)), digits = 3))
sink()

### graphical comparison
perfplot(BostonHousingRMSE, file = "BostonHousing_RMSE.pdf", lab = "RMSE")
perfplot(BostonHousingNPAR, file = "BostonHousing_NPAR.pdf", lab = "Complexity")
