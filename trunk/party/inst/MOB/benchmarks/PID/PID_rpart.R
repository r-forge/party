
library("rpart")
source("../npar.R")

load("PIDBootstrap.rda")
n <- nrow(bdata$data)
B <- ncol(bdata$boot)

error <- numeric(B)
np <- numeric(B)

for (b in 1:B) {

    bs <- bdata$boot[,b]
    oob <- bs == 0
    ldata <- bdata$data[rep(1:n, bs),]
    tdata <- bdata$data[oob,]

    fm <- rpart(bdata$fm, data = ldata)
    if (nrow(fm$cptable) > 1) {
        cpopt <- fm$cptable[which.min(fm$cptable[,"xerror"]), 1]
        fm <- prune(fm, cp = cpopt)
    }

    yhat <- predict(fm, newdata = tdata, type = "class")
    error[b] <- mean(yhat != tdata[,bdata$response])
    np[b] <- npar(fm)
    cat("b: ", b, " error: ", error[b], " #par: ", np[b], "\n")
}

save(error, np, file = "PID_rpart_error.rda")

