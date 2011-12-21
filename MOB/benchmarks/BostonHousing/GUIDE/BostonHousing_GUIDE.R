
source("R2guide.R")

load("../BostonHousingBootstrap.rda")

rval <- foo(bdata$data, bdata$response, bdata$boot)
error <- rval$error
np <- rval$npar
summary(error)
summary(np)

save(error, np, file = "../BostonHousing_guide_error.rda")
