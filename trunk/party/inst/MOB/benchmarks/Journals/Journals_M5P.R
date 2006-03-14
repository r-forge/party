
library("RWeka")
source("../npar.R")

load("JournalsBootstrap.rda")
n <- nrow(bdata$data)
B <- ncol(bdata$boot)

error <- numeric(B)
np <- numeric(B)

for (b in 1:B) {

    bs <- bdata$boot[,b]
    oob <- bs == 0
    ldata <- bdata$data[rep(1:n, bs),]
    tdata <- bdata$data[oob,]

    fm <- M5P(bdata$fm, data = ldata)
    yhat <- predict(fm, newdata = tdata)
    error[b] <- mean((yhat  - tdata[,bdata$response])^2)
    cat("b: ", b, " error: ", error[b], "\n")
    np[b] <- npar(fm)
}

save(error, np, file = "Journals_M5P_error.rda")

