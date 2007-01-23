#########################
## auxiliary functions ##
#########################

sim <- function(n = 500, ...) {
  rvalMOB <- matrix(rep(NA, n * 3), ncol = 3)
  rvalGUIDE <- matrix(rep(NA, n * 3), ncol = 3)
  colnames(rvalMOB) <- colnames(rvalGUIDE) <- c("nsplit", "firstvar", "firstpoint")
  for(i in 1:n) {
    dat <- dgp(...)
    rvalMOB[i,] <- mymob(dat)
    rvalGUIDE[i,] <- myGUIDE(dat)
  }

  rvalMOB <- as.data.frame(rvalMOB)
  rvalGUIDE <- as.data.frame(rvalGUIDE)
  rvalMOB$nsplit <- factor(rvalMOB$nsplit)
  rvalGUIDE$nsplit <- factor(rvalGUIDE$nsplit)
  rvalMOB$firstvar <- factor(rvalMOB$firstvar)
  rvalGUIDE$firstvar <- factor(rvalGUIDE$firstvar)
  
  return(list(mob = rvalMOB, GUIDE = rvalGUIDE))
}

dgp <- function(nobs = 500, point = 0, coef = c(0, 1))
{
  coef <- rep(coef, length.out = 2)
  rval <- data.frame(x = rnorm(nobs),  z1 = rnorm(nobs), z2 = rnorm(nobs),
                     z3 = rnorm(nobs), z4 = round(runif(nobs), digits = 1),
		     z5 = factor(sample(1:2, replace = TRUE)),
		     z6 = factor(sample(1:5, replace = TRUE)))
  rval$y <- if(isTRUE(all.equal(coef, c(0, 0)))) rnorm(nobs)
    else ifelse(rval$z1 > point, coef[1] + coef[2] * rval$x, 0) + rnorm(nobs)
  return(rval)
}

mymob <- function(data) {
  fm <- try(mob(y ~ x | z1 + z2 + z3 + z4 + z5 + z6, data = data, model = linearModel,
    control = mob_control(alpha = 0.05, minsplit = 50, breakties = TRUE)), silent = TRUE)
  if(inherits(fm, "try-error")) {
    rval <- c(0, NA, NA)
  } else {  
    rval <- if(firstvar(fm) > 1) c(nsplit(fm), firstvar(fm), NA)
      else c(nsplit(fm), 1, firstpoint(fm))
  }
  names(rval) <- c("nsplit", "firstvar", "firstpoint")
  return(rval)
}

myGUIDE <- function(data) {
  fm <- GUIDE(y ~ x | z1 + z2 + z3 + z4 + z5 + z6, data = data)
  if(nsplit(fm) < 1) {
    rval <- c(0, NA, NA)
  } else {
    rval <- if(firstvar(fm) > 1) c(nsplit(fm), firstvar(fm), NA)
      else c(nsplit(fm), 1, firstpoint(fm))
  }
  names(rval) <- c("nsplit", "firstvar", "firstpoint")
  return(rval)
}

nsplit <- function(obj) {
  if(class(obj) == "GUIDE") {
    obj <- obj$tree
    nl <- obj[grep("Number of terminal nodes of final tree", obj)]
    nl <- as.numeric(strsplit(nl, ": ")[[1]][2])
    nl <- nl - 1
  } else {
    nl <- NROW(coef(obj)) - 1
  }
  return(nl)
}

firstvar <- function(obj) {
  if(nsplit(obj) < 1) return(NA)
  if(class(obj) == "GUIDE") {
    obj <- obj$tree
    rval <- obj[grep("Regression tree:", obj) + 2]
    rval <- gsub("<", "", rval)
    rval <- strsplit(rval, "=")[[1]][1]
    rval <- as.numeric(strsplit(rval, "1: z")[[1]][2])
  } else {
    rval <- as.numeric(obj@tree$psplit$variableID)
  }
  return(rval)
}

firstpoint <- function(obj) {
  if(nsplit(obj) < 1) return(NA)
  if(class(obj) == "GUIDE") {
    obj <- obj$tree
    rval <- obj[grep("Regression tree:", obj) + 2]
    rval <- as.numeric(strsplit(rval, "= ")[[1]][2])
  } else {
    rval <- as.numeric(obj@tree$psplit$splitpoint)
  }
  return(rval)
}

###########################
## size/power simulation ##
###########################

library("party")
source("LohTools")
source("dgp2.R")

set.seed(20012007)

sim00 <- sim(coef = c(0, 0))
print(summary(sim00$mob))
print(summary(sim00$GUIDE))

sim01 <- sim(coef = c(0, 1))
print(summary(sim01$mob))
print(summary(sim01$GUIDE))

sim10 <- sim(coef = c(1, 0))
print(summary(sim10$mob))
print(summary(sim10$GUIDE))

sim11 <- sim(coef = c(1, 1))
print(summary(sim11$mob))
print(summary(sim11$GUIDE))

save(sim00, sim01, sim10, sim11, file = "size-power.rda")
