
library(ipred)
library(party)
library(exactRankTests)
source("./dendrogram.R")

data(BostonHousing)
data(Ozone)
data(GlaucomaM)

MFF = ModelFrame(medv ~ ., data = BostonHousing)

VarList = treedesign(MFF)
VarList@control = new("GrowControl", minsplit = 20, 
                       minprob = 0.1,
                       minstat = qnorm(1 - 0.05/ncol(MFF@input)),
                       varnull = 1e-10)

Rprof("tree")
x = stree(VarList)
Rprof(NULL)
show(x)
dx = as.dendrogram(x)
save(dx, file = "dendrogram.Rda")
plot(dx, center = TRUE)
plot(dx)

index = c(51,100,120,29,150)
pr = treepredict(x, VarList, MFF@input[index,])

a = rpart(medv ~ ., data = BostonHousing)
print(a)
predict(a, newdata = BostonHousing[index,])
pr
MFF@response[index,]

pr = treepredict(x, VarList, MFF@input)[,1]
pro = predict(a)        
re = range(c(pr, pro))
#par(mfrow = c(1,3))
plot(pr + rnorm(length(pr)), pro + rnorm(length(pro)), 
     xlim = re, ylim = re)

plot(pr + rnorm(length(pr)), MFF@response[[1]] + rnorm(length(pro)), 
     xlim = re, ylim = re)

plot(pro + rnorm(length(pro)), MFF@response[[1]] + rnorm(length(pro)), 
     xlim = re, ylim = re)


myOzone = Ozone[complete.cases(Ozone$V4),]
MFF = ModelFrame(V4 ~ ., data = myOzone, na.action = na.pass)

VarList = treedesign(MFF)
VarList@control = new("GrowControl", minsplit = 20, 
                       minprob = 0.1,
                       minstat = qnorm(1 - 0.05/ncol(MFF@input)),
                       varnull = 1e-10)

# Rprof("tree")
x = stree(VarList)
# Rprof(NULL)
show(x)
dx = as.dendrogram(x)
plot(dx, center = TRUE)
plot(dx)

index = c(51,100,120,29,150)
pr = treepredict(x, VarList, MFF@input[index,])

a = rpart(V4 ~ ., data  = myOzone)
print(a)
predict(a, newdata = myOzone[index,])
pr
MFF@response[index,]

MFF = ModelFrame(Class ~ ., data = GlaucomaM)


VarList = treedesign(MFF)
VarList@control = new("GrowControl", minsplit = 20, 
                       minprob = 0.1,
                       minstat = qnorm(1 - 0.05/ncol(MFF@input)),
                       varnull = 1e-10)

# Rprof("tree")
x = stree(VarList)
# Rprof(NULL)
dx = as.dendrogram(x)
plot(dx, center = TRUE)
plot(dx)

show(x)

index = c(51,100,120,29,150)
pr = treepredict(x, VarList, MFF@input[index,])

a = rpart(Class ~ ., data = GlaucomaM)
print(a)
predict(a, newdata = GlaucomaM[index,])
pr
MFF@response[index,]


data(GBSG2)
res = cscores(Surv(GBSG2$time, GBSG2$cens))
dat = GBSG2[,-(9:10)]
MFF = ModelFrame(res ~ ., data = dat)


VarList = treedesign(MFF)
VarList@control = new("GrowControl", minsplit = 20, 
                       minprob = 0.1,
                       minstat = qnorm(1 - 0.05/ncol(MFF@input)),
                       varnull = 1e-10)

# Rprof("tree")
x = stree(VarList)
# Rprof(NULL)

dx = as.dendrogram(x)
plot(dx, center = TRUE)
plot(dx)

show(x)

Rprof("pred")
pr = treepredict(x, VarList, MFF@input)
Rprof(NULL)

VarList@control@minstat = 2.4
show(stree(VarList))

