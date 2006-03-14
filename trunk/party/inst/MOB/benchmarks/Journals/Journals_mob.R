
library("party")
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

    fm <- try(mob(bdata$mobfm, data = ldata, control = mob_control(minsplit = 10), 
              model = linearModel))
    if (!inherits(fm, "try-error")) {
        yhat <- predict(fm, newdata = tdata)
        error[b] <- mean((yhat  - tdata[,bdata$response])^2)
        np[b] <- npar(fm)
    } else {
        error[b] <- NA
        np[b] <- NA
    }
    cat("b: ", b, " error: ", error[b], " #par: ", np[b], "\n")
}

save(error, np, file = "Journals_mob_error.rda")
