
set.seed(290875)

load("../../journals.rda")

n <- nrow(journals)
B <- 100

journals_bs <- rmultinom(B, n, rep(1, n) / n)

bdata <- list(data = journals,
              mobfm = subs ~ citeprice | society + citations + age + chars + price,
              fm = subs ~ citeprice + society + citations + age + chars + price,
              boot = journals_bs,
              response = "subs")

save(bdata, file = "JournalsBootstrap.rda")
