
source("R2guide.R")

load("../JournalsBootstrap.rda")

rval <- foo(bdata$data, bdata$response, bdata$boot)
error <- rval$error
np <- rval$npar
summary(error)
summary(np)

save(error, np, file = "../Journals_guide_error.rda")
