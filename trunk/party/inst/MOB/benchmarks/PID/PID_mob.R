
library("party")

load("PIDBootstrap.rda")
n <- nrow(bdata$data)
B <- ncol(bdata$boot)

error <- numeric(B)

for (b in 1:B) {

    bs <- bdata$boot[,b]
    oob <- bs == 0
    ldata <- bdata$data[rep(1:n, bs),]
    tdata <- bdata$data[oob,]

    fm <- try(mob(bdata$mobfm, data = ldata, 
                  control = mob_control(minsplit = 40),
                  family = binomial(), model = glinearModel))
    if (!inherits(fm, "try-error")) {
        yhat <- factor(predict(fm, newdata = tdata, type = "response") > 0.5,
                       labels = levels(tdata[,bdata$response]))              
        error[b] <- mean(yhat != tdata[,bdata$response])
    } else {
        error[b] <- NA
    }
    cat("b: ", b, " error: ", error[b], "\n")
}

save(error, file = "PID_mob_error.rda")

