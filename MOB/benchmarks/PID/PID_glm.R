
load("PIDBootstrap.rda")
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

    fm <- glm(bdata$fm, data = ldata, family = binomial())
    yhat <- factor(predict(fm, newdata = tdata, type = "response") > 0.5, 
                   labels = levels(tdata[,bdata$response]))
    error[b] <- mean(yhat != tdata[,bdata$response])
    np[b] <- npar(fm)
    cat("b: ", b, " error: ", error[b], " #par: ", np[b], "\n")
}

save(error, np, file = "PID_glm_error.rda")

