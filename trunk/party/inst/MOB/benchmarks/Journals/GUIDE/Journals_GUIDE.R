
source("R2guide.R")

load("../JournalsBootstrap.rda")

error <- foo(bdata$data, bdata$response, bdata$boot)
print(error)
summary(error)

save(error, file = "../Journals_guide_error.rda")
