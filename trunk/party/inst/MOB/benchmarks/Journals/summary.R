
source("../perfplot.R")
library("multcomp")

rdafiles <- list.files(pattern = "error")

JournalsMSE <- c()
model <- c()

for (f in rdafiles) {

    load(f)
    print(f)
    print(summary(error))
    JournalsMSE <- cbind(JournalsMSE, sqrt(error))
    model <- c(model, strsplit(f, "_")[[1]][2])
    if (exists("np")) {
        print(summary(np))
        rm(np)
    }
}

JournalsMSE <- as.data.frame(JournalsMSE)
colnames(JournalsMSE) <- model

summary(JournalsMSE)

### re-order and omit lm
JournalsMSE <- JournalsMSE[,c(5, 1, 3, 2, 6)]
perfplot(JournalsMSE, file = "Journals_MSE.pdf", boxplot = TRUE, 
         ylab = "RMSE")


### alignment
tmp <- JournalsMSE - apply(JournalsMSE, 1, mean)

tmp <- data.frame(error = unlist(tmp),
                  model = factor(rep(colnames(tmp), 
                                 rep(nrow(tmp), 
                                     ncol(tmp)))))

si <- simint(error ~ model, data = tmp, type = "Dunnett", 
             base = which(levels(tmp$model) == "mob"))
rownames(si$estimate) <- c("M5'", "ctree", "GUIDE", "rpart")

pdf("Journals_CI.pdf")
plot(si, main = "", xlab = "RMSE difference")
dev.off()
