
library("party")

load("BostonHousingBootstrap.rda")
n <- nrow(bdata$data)
B <- ncol(bdata$boot)

bdata$data$rad <- factor(bdata$data$rad, ordered = TRUE)

error <- numeric(B)

for (b in 1:B) {

    bs <- bdata$boot[,b]
    oob <- bs == 0
    ldata <- bdata$data[rep(1:n, bs),]
    tdata <- bdata$data[oob,]

    fm <- try(mob(bdata$mobfm, data = ldata, control = mob_control(minsplit = 40), 
              model = linearModel))
    if (!inherits(fm, "try-error")) {
        yhat <- predict(fm, newdata = tdata)
        error[b] <- mean((yhat  - tdata[,bdata$response])^2)
    } else {
        error[b] <- NA
    }
    cat("b: ", b, " error: ", error[b], "\n")
}

save(error, file = "BostonHousing_mob_error.rda")
