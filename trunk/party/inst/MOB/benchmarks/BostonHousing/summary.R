
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

perfplot(BostonHousingMSE, file = "BostonHousing_MSE.pdf", boxplot = TRUE, 
         ylab = "RMSE")

tmp <- data.frame(error = unlist(BostonHousingMSE),
                  model = factor(rep(colnames(BostonHousingMSE), 
                                 rep(nrow(BostonHousingMSE), 
                                     ncol(BostonHousingMSE)))))

si <- simint(error ~ model, data = tmp, type = "Dunnett", base = which(levels(tmp$model) == "mob"))

pdf("BostonHousing_CI.pdf")
plot(si, xlim = c(-0.1, 1.2))
dev.off()




