
source("R2quest.R")

load("../PIDBootstrap.rda")

error <- foo(bdata$data, bdata$response, bdata$boot)
print(error)
summary(error)

save(error, file = "../PID_quest_error.rda")
     