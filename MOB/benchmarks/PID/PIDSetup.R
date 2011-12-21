
set.seed(290875)

data("PimaIndiansDiabetes2", package = "mlbench")
PimaIndiansDiabetes <- na.omit(PimaIndiansDiabetes2[,-c(4, 5)])

n <- nrow(PimaIndiansDiabetes)
B <- 250

PimaIndiansDiabetes_bs <- rmultinom(B, n, rep(1, n) / n)

bdata <- list(data = PimaIndiansDiabetes,
              mobfm = diabetes ~ glucose | pregnant + pressure + mass + pedigree + age,
              fm = diabetes ~ glucose + pregnant + pressure + mass + pedigree + age,
              boot = PimaIndiansDiabetes_bs,
              response = "diabetes")

save(bdata, file = "PIDBootstrap.rda")

