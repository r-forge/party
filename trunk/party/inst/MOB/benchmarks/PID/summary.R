
source("../perfplot.R")
library("multcomp")

rdafiles <- list.files(pattern = "error")

PIDmc <- c()
model <- c()

for (f in rdafiles) {

    load(f)
    PIDmc <- cbind(PIDmc, error)
    model <- c(model, strsplit(f, "_")[[1]][2])
}

PIDmc <- as.data.frame(PIDmc)
colnames(PIDmc) <- model

summary(PIDmc)

PIDmc <- PIDmc[complete.cases(PIDmc),]

perfplot(PIDmc, file = "PID_MC.pdf", boxplot = TRUE, 
         ylab = "Misclass")

### alignment
tmp <- PIDmc - apply(PIDmc, 1, mean)

tmp <- data.frame(error = unlist(tmp),
                  model = factor(rep(colnames(tmp),
                                 rep(nrow(tmp),    
                                     ncol(tmp))))) 

si <- simint(error ~ model, data = tmp, type = "Dunnett",
             base = which(levels(tmp$model) == "mob"))   

pdf("PIDmc_CI.pdf")
plot(si)
dev.off()

