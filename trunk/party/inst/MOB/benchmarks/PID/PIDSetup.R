
set.seed(290875)

data("PimaIndiansDiabetes", package = "mlbench")

n <- nrow(PimaIndiansDiabetes)
B <- 250

PimaIndiansDiabetes_bs <- rmultinom(B, n, rep(1, n) / n)

bdata <- list(data = PimaIndiansDiabetes,
              mobfm = diabetes ~ glucose | pregnant + pressure + triceps + 
                      insulin + mass + pedigree + age,
              fm = diabetes ~ glucose + pregnant + pressure + triceps +
                   insulin + mass + pedigree + age,
              boot = PimaIndiansDiabetes_bs,
              response = "diabetes")

save(bdata, file = "PIDBootstrap.rda")

