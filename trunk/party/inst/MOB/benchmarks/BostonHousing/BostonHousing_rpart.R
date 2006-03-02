
library("rpart")

load("BostonHousingBootstrap.rda")
n <- nrow(bdata$data)
B <- ncol(bdata$boot)

error <- numeric(B)

for (b in 1:B) {

    bs <- bdata$boot[,b]
    oob <- bs == 0
    ldata <- bdata$data[rep(1:n, bs),]
    tdata <- bdata$data[oob,]

    fm <- rpart(bdata$fm, data = ldata)
    cpopt <- fm$cptable[which.min(fm$cptable[,"xerror"]), 1]
    fm <- prune(fm, cp = cpopt)

    yhat <- predict(fm, newdata = tdata)
    error[b] <- mean((yhat  - tdata[,bdata$response])^2)
    cat("b: ", b, " error: ", error[b], "\n")
}

save(error, file = "BostonHousing_rpart_error.rda")

