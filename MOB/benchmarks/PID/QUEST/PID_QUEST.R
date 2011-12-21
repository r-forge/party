
source("R2quest.R")

load("../PIDBootstrap.rda")

rval <- foo(bdata$data, bdata$response, bdata$boot)
error <- rval$error
np <- rval$npar
summary(error)
summary(np)

save(error, np, file = "../PID_quest_error.rda")
     
