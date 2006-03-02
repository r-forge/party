
set.seed(290875)

data("BostonHousing", package = "mlbench")

BostonHousing$lstat <- log(BostonHousing$lstat)
BostonHousing$rm <- BostonHousing$rm^2

## as well as partitioning variables (for fluctuation testing)
BostonHousing$chas <- factor(BostonHousing$chas, levels = 0:1, labels = c("no", "yes"))

n <- nrow(BostonHousing)
B <- 100

BostonHousing_bs <- rmultinom(B, n, rep(1, n) / n)

bdata <- list(data = BostonHousing,
              mobfm = medv ~ lstat + rm | zn + indus + chas + nox + age +
                      dis + rad + tax + crim + b + ptratio,
              fm = medv ~ lstat + rm + zn + indus + chas + nox + age + dis +
                   rad + tax + crim + b + ptratio,
              boot = BostonHousing_bs,
              response = "medv")

save(bdata, file = "BostonHousingBootstrap.rda")
