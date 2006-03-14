source("../perfplot.R")
library("multcomp")

rdafiles <- list.files(pattern = "error")

PID_MC <- c()
PID_NPAR <- c()
model <- c()

for (f in rdafiles) {
    load(f)    
    PID_MC <- cbind(PID_MC, error)
    PID_NPAR <- cbind(PID_NPAR, np)
    model <- c(model, strsplit(f, "_")[[1]][2])
}

if(identical(model, c("J48", "LMT", "ctree", "glm", "mob", "quest", "rpart"))) {
  model <- c("J4.8", "LMT", "CTree", "GLM", "MOB", "QUEST", "RPart")
  o <- c(5, 2, 3, 6, 1, 7)
} else {
  o <- 1:length(model)
}

colnames(PID_MC) <- model
colnames(PID_NPAR) <- model
PID_MC <- as.data.frame(PID_MC[,o])
PID_NPAR <- as.data.frame(PID_NPAR[,o])
PID_MC <- na.omit(PID_MC)
PID_NPAR <- na.omit(PID_NPAR)

## median performance
sink(file = "PID_results.txt")
print(round(cbind(Misclassification = sapply(PID_MC, median),
  Complexity = sapply(PID_NPAR, median)), digits = 3))
sink()

### graphical comparison
perfplot(PID_MC, file = "PID_MC.pdf", lab = "Misclassification")
perfplot(PID_NPAR, file = "PID_NPAR.pdf", lab = "Complexity")
