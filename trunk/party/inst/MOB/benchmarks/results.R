###############################
## collect bootstrap results ##
###############################
setwd("Journals/")
rdafiles <- list.files(pattern = "error")
JournalsRMSE <- c()
JournalsNPAR <- c()
model <- c()
for (f in rdafiles) {
    load(f)    
    JournalsRMSE <- cbind(JournalsRMSE, sqrt(error))
    JournalsNPAR <- cbind(JournalsNPAR, np)
    model <- c(model, strsplit(f, "_")[[1]][2])
}
if(identical(model, c("M5P", "ctree", "guide", "lm", "mob", "rpart"))) {
  model <- c("M5'", "CTree", "GUIDE", "lm", "MOB", "RPart")
  o <- c(5, 3, 1, 2, 6)
} else {
  o <- 1:length(model)
}
colnames(JournalsRMSE) <- model
colnames(JournalsNPAR) <- model
JournalsRMSE <- as.data.frame(JournalsRMSE[,o])
JournalsNPAR <- as.data.frame(JournalsNPAR[,o])

setwd("../BostonHousing/")
rdafiles <- list.files(pattern = "error")
BostonHousingRMSE <- c()
BostonHousingNPAR <- c()
model <- c()
for (f in rdafiles) {
    load(f)    
    BostonHousingRMSE <- cbind(BostonHousingRMSE, sqrt(error))
    BostonHousingNPAR <- cbind(BostonHousingNPAR, np)
    model <- c(model, strsplit(f, "_")[[1]][2])
}
if(identical(model, c("M5P", "ctree", "guide", "lm", "mob", "rpart"))) {
  model <- c("M5'", "CTree", "GUIDE", "lm", "MOB", "RPart")
  o <- c(5, 3, 1, 2, 6)
} else {
  o <- 1:length(model)
}
colnames(BostonHousingRMSE) <- model
colnames(BostonHousingNPAR) <- model
BostonHousingRMSE <- as.data.frame(BostonHousingRMSE[,o])
BostonHousingNPAR <- as.data.frame(BostonHousingNPAR[,o])

setwd("../PID/")
rdafiles <- list.files(pattern = "error")
PID_MC <- c()
PID_NPAR <- c()
model <- c()
for (f in rdafiles) {
    load(f)    
    PID_MC <- cbind(PID_MC, error)
    PID_NPAR <- cbind(PID_NPAR, np)
    model <- c(model, strsplit(f, "_")[[1]][2])
}
if(identical(model, c("J48", "LMT", "cruise", "ctree", "glm", "mob", "quest", "rpart"))) {
  model <- c("J4.8", "LMT", "CRUISE", "CTree", "GLM", "MOB", "QUEST", "RPart")
  o <- c(6, 2, 4, 7, 3, 1, 8)
} else {
  o <- 1:length(model)
}
colnames(PID_MC) <- model
colnames(PID_NPAR) <- model
PID_MC <- as.data.frame(PID_MC[,o])
PID_NPAR <- as.data.frame(PID_NPAR[,o])
PID_MC <- na.omit(PID_MC)
PID_NPAR <- na.omit(PID_NPAR)

setwd("../")

##########################
## fit empirical models ##
##########################

library("party")
library("RWeka")
library("rpart")
library("LohTools")
load("../journals.rda")
fmJ_M5P <- M5P(subs ~ citeprice + society + citations + age + chars + price, data = journals)
fmJ_CTree <- ctree(subs ~ citeprice + society + citations + age + chars + price, data = journals)
fmJ_RPart <- rpart(subs ~ citeprice + society + citations + age + chars + price, data = journals)
fmJ_RPart <- prune(fmJ_RPart, cp = fmJ_RPart$cptable[which.min(fmJ_RPart$cptable[,"xerror"]), 1])
fmJ_GUIDE <- GUIDE(subs ~ citeprice + society + citations + age + chars + price, data = journals)
fmJ <- mob(subs ~ citeprice | society + citations + age + chars + price, data = journals,
  control = mob_control(minsplit = 10), model = linearModel)

data("BostonHousing", package = "mlbench")
BostonHousing$lstat <- log(BostonHousing$lstat)
BostonHousing$rm <- BostonHousing$rm^2
BostonHousing$chas <- factor(BostonHousing$chas, levels = 0:1, labels = c("no", "yes"))
fmBH_M5P <- M5P(medv ~ lstat + rm + zn + indus + chas + nox + age + dis + rad + tax + crim + b + ptratio, data = BostonHousing)
fmBH_CTree <- ctree(medv ~ lstat + rm + zn + indus + chas + nox + age + dis + rad + tax + crim + b + ptratio, data = BostonHousing)
fmBH_RPart <- rpart(medv ~ lstat + rm + zn + indus + chas + nox + age + dis + rad + tax + crim + b + ptratio, data = BostonHousing)
fmBH_RPart <- prune(fmBH_RPart, cp = fmBH_RPart$cptable[which.min(fmBH_RPart$cptable[,"xerror"]), 1])
fmBH_GUIDE <- GUIDE(medv ~ lstat + rm + zn + indus + chas + nox + age + dis + rad + tax + crim + b + ptratio, data = BostonHousing)
BostonHousing$rad <- factor(BostonHousing$rad, ordered = TRUE)
fmBH <- mob(medv ~ lstat + rm | zn + indus + chas + nox + age + dis + rad + tax + crim + b + ptratio,
  data = BostonHousing, control = mob_control(minsplit = 40), model = linearModel)

data("PimaIndiansDiabetes2", package = "mlbench")
PimaIndiansDiabetes <- na.omit(PimaIndiansDiabetes2[,-c(4, 5)])
fmPID_J48 <- J48(diabetes ~ glucose + pregnant + pressure + mass + pedigree + age, data = PimaIndiansDiabetes)
fmPID_LMT <- LMT(diabetes ~ glucose + pregnant + pressure + mass + pedigree + age, data = PimaIndiansDiabetes)
fmPID_CTree <- ctree(diabetes ~ glucose + pregnant + pressure + mass + pedigree + age, data = PimaIndiansDiabetes)
fmPID_RPart <- rpart(diabetes ~ glucose + pregnant + pressure + mass + pedigree + age, data = PimaIndiansDiabetes)
fmPID_RPart <- prune(fmPID_RPart, cp = fmPID_RPart$cptable[which.min(fmPID_RPart$cptable[,"xerror"]), 1])
setwd("PID/QUEST/")
source("R2quest.R")
fmPID_QUEST <- foo(PimaIndiansDiabetes, "diabetes", matrix(1, ncol = 1, nrow = NROW(PimaIndiansDiabetes)), obvious = TRUE)
setwd("../../")
fmPID_CRUISE <- CRUISE(diabetes ~ glucose + pregnant + pressure + mass + pedigree + age, data = PimaIndiansDiabetes)
fmPID <- mob(diabetes ~ glucose | pregnant + pressure + mass + pedigree + age,
  data = PimaIndiansDiabetes, control = mob_control(minsplit = 40),
  family = binomial(), model = glinearModel)

data("GBSG2", package = "ipred")
nloglik <- function(x) -logLik(x)
GBSG2$time <- GBSG2$time/365
GBSG2$pTh <- (as.numeric(GBSG2$horTh)-1) * GBSG2$pnodes
fmGBSG2 <- mob(Surv(time, cens) ~ horTh + pnodes | progrec + menostat + estrec + menostat + age + tsize + tgrade,
  data = GBSG2, model = survReg, control = mob_control(objfun = nloglik, minsplit = 40))


###############################
## compute empirical results ##
###############################

source("npar.R")
rmse <- function(obj, obs) sqrt(mean((obs - predict(obj))^2))
mc <- function(obj, obs, type = "response") {
  yhat <- predict(obj, data = PimaIndiansDiabetes, type = type)
  if(is.numeric(yhat)) yhat <- yhat > 0.5
  if(!is.factor(yhat)) yhat <- factor(yhat, labels = levels(obs))
  mean(yhat != obs)
}

JournalsObvious <- sapply(list(fmJ, fmJ_GUIDE, fmJ_M5P, fmJ_CTree, fmJ_RPart),
  function(obj) c(rmse(obj, journals$subs), npar(obj)))
colnames(JournalsObvious) <- names(JournalsRMSE)

BostonHousingObvious <- sapply(list(fmBH, fmBH_GUIDE, fmBH_M5P, fmBH_CTree, fmBH_RPart),
  function(obj) c(rmse(obj, BostonHousing$medv), npar(obj)))
colnames(BostonHousingObvious) <- names(BostonHousingRMSE)

PID_Obvious <- rbind(c(mc(fmPID, PimaIndiansDiabetes$diabetes),
  mc(fmPID_LMT, PimaIndiansDiabetes$diabetes, type = "class"),
  mc(fmPID_CTree, PimaIndiansDiabetes$diabetes),
  fmPID_QUEST$error,
  mean(predict(fmPID_CRUISE) != PimaIndiansDiabetes$diabetes),
  mc(fmPID_J48, PimaIndiansDiabetes$diabetes, type = "class"),
  mc(fmPID_RPart, PimaIndiansDiabetes$diabetes, type = "class")),
  c(npar(fmPID), npar(fmPID_LMT), npar(fmPID_CTree), fmPID_QUEST$npar, npar(fmPID_CRUISE), 
    npar(fmPID_J48), npar(fmPID_RPart)))
colnames(PID_Obvious) <- names(PID_MC)


######################
## store everything ##
######################

save(JournalsRMSE, JournalsNPAR, BostonHousingRMSE, BostonHousingNPAR, PID_MC, PID_NPAR,
  JournalsObvious, BostonHousingObvious, PID_Obvious, file = "results.rda")
save(fmJ, fmBH, fmPID, fmGBSG2, file = "MOB-fit.rda")
