
set.seed(290875)
gctorture(on = FALSE)
library(party)
if (!require(maxstat))
    stop("cannot load package maxstat")
gctorture(on = GCtorture)

### load additional R code which is only partially in arty/R'
source(file.path(.find.package("party"), "Rcode", "TestCode.R"))

### 
###
###    Regression tests for cutpoint search
###    
###    functions defined in file `./src/Splits.c'    

### tests for function C_Split
x <- rnorm(6000)
y <- rnorm(6000)
weights <- rep(1, length(x))
splitctrl <- new("SplitControl")
split <- Split(x, y, weights, splitctrl)
mydata <- data.frame(y, x)
ms <- maxstat.test(y ~ x, data = mydata, smethod = "Data", pmethod = "none")
stopifnot(isequal(split[[1]], ms$estimate))
stopifnot(isequal(split[[2]], ms$statistic))
stopifnot(isequal(max(split[[3]]), ms$statistic))

### Hohnloser data
data(hohnloser)
ms <-  maxstat.test(Surv(month, cens) ~ EF, data = hohnloser,
                    smethod = "LogRank")
splitctrl <- new("SplitControl")
splitctrl@minprob <- 0.1
splitctrl@minsplit <- as.integer(5)

split <- Split(hohnloser$EF, cscores(Surv(hohnloser$month, hohnloser$cens)),
               rep(1, nrow(hohnloser)), splitctrl)
stopifnot(isequal(split[[1]], ms$estimate))
stopifnot(isequal(split[[2]], ms$statistic))
stopifnot(isequal(max(split[[3]]), ms$statistic))

### categorical splits
n <- 100
xf <- gl(5, 100/5)
yf <- gl(4, 100/4)[sample(1:length(xf))]
weights <- rep(1, length(xf))
splitctrl <- new("SplitControl")
splitctrl@minprob <- 0.1
splitctrl@minsplit <- as.integer(5)
split <- Split(xf, yf, weights, splitctrl)
split

### Check if the statistic used for selecting the split is
### correct: For the ranks of a continuous response the statistic
### needs to be equal to the standardized Wilcoxon statistic

y <- rnorm(100) + c(rep(0, 25), rep(1, 25), rep(0, 25), rep(1, 25))
x <- gl(4, 25)
weights <- rep(1, length(y))
split <- Split(x, rank(y), weights, splitctrl)
levelset <- levels(x)[split[[4]] == 1]
tstat <- split[[2]]
p <- wilcox.test(y ~ I(x %in% levelset),corr = FALSE,
                alternative = "less")$p.value
stopifnot(isequal(round(abs(qnorm(p)), 6), round(tstat, 6)))

y <- rnorm(100) + c(rep(0, 25), rep(1, 25), rep(0, 25), rep(1, 25))
x <- rnorm(100)
weights <- rep(1, length(y))
split <- Split(x, rank(y), weights, splitctrl)
tstat <- split[[2]]
p <- wilcox.test(y ~ I(x <= split[[1]]), corr = FALSE,
                alternative = "less")$p.value
stopifnot(isequal(round(abs(qnorm(p)), 6), round(tstat, 6)))
