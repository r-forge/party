
load("BostonHousingBootstrap.rda")
source("../npar.R")
n <- nrow(bdata$data)
B <- ncol(bdata$boot)

error <- numeric(B)
np <- numeric(B)

for (b in 1:B) {

    bs <- bdata$boot[,b]
    oob <- bs == 0
    ldata <- bdata$data[rep(1:n, bs),]
    tdata <- bdata$data[oob,]

    fm <- lm(bdata$fm, data = ldata)
    yhat <- predict(fm, newdata = tdata)
    error[b] <- mean((yhat  - tdata[,bdata$response])^2)
    np[b] <- npar(fm)
    cat("b: ", b, " error: ", error[b], " #par: ", np[b], "\n")
}

save(error, np, file = "BostonHousing_lm_error.rda")

