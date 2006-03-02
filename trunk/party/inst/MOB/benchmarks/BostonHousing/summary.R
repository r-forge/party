
source("../perfplot.R")
library("multcomp")

rdafiles <- list.files(pattern = "error")

BostonHousingMSE <- c()
model <- c()

for (f in rdafiles) {

    load(f)
    print(f)
    print(summary(error))
    BostonHousingMSE <- cbind(BostonHousingMSE, sqrt(error))
    model <- c(model, strsplit(f, "_")[[1]][2])
}

BostonHousingMSE <- as.data.frame(BostonHousingMSE)
colnames(BostonHousingMSE) <- model

summary(BostonHousingMSE)

### re-order and omit lm
BostonHousingMSE <- BostonHousingMSE[,c(5, 1, 3, 2, 6)]
perfplot(BostonHousingMSE, file = "BostonHousing_MSE.pdf", boxplot = TRUE, 
         ylab = "RMSE")


### alignment
tmp <- BostonHousingMSE - apply(BostonHousingMSE, 1, mean)

tmp <- data.frame(error = unlist(tmp),
                  model = factor(rep(colnames(tmp), 
                                 rep(nrow(tmp), 
                                     ncol(tmp)))))

si <- simint(error ~ model, data = tmp, type = "Dunnett", 
             base = which(levels(tmp$model) == "mob"))
rownames(si$estimate) <- c("M5'", "GUIDE", "ctree", "rpart")

pdf("BostonHousing_CI.pdf")
plot(si, xlim = c(-0.1, 1.2), main = "", xlab = "RMSE difference")
dev.off()
