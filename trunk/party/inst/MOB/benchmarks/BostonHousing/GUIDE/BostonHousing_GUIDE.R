
source("R2guide.R")

load("../BostonHousingBootstrap.rda")

error <- foo(bdata$data, bdata$response, bdata$boot)
print(error)
summary(error)

save(error, file = "../BostonHousing_guide_error.rda")
