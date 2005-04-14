
set.seed(290875)
gctorture(on = FALSE)
library(party)
if (!require(ipred))
    stop("cannot load package ipred")
gctorture(on = GCtorture)

### load additional R code which is only partially in arty/R'
source(file.path(.find.package("party"), "Rcode", "TestCode.R"))

gtctrl <- new("GlobalTestControl")
tlev <- levels(gtctrl@testtype)

data(GlaucomaM, package = "ipred")
inp <- party:::initVariableFrame(GlaucomaM[,-63,drop = FALSE]) #, fun = rank)
resp <- party:::initVariableFrame(GlaucomaM[,"Class",drop = FALSE])
ls <- new("LearningSample", inputs = inp, responses = resp,
          weights = rep(1, inp@nobs), nobs = nrow(GlaucomaM), 
          ninputs = inp@ninputs)
tm <- TreeFitMemory(ls, TRUE)
ctrl <- new("TreeControl")
ctrl@gtctrl@testtype <- factor("Bonferroni", levels = tlev)
ctrl@varctrl@teststattype <- factor("quadform", levels = c("maxabs", "quadform"))
where <- rep(0, ls@nobs)
storage.mode(where) <- "integer"
tree <- .Call("R_TreeGrow", ls, ls@weights, tm, ctrl, where)
stopifnot(isequal(tree[[5]][[3]], 0.059))

# print(tree)

ctrl@tgctrl@stump = TRUE
stump <- .Call("R_TreeGrow", ls, ls@weights, tm, ctrl, where)
# print(stump)

data(treepipit)

tr <- ctree(counts ~ ., data = treepipit, teststattype = "quadform")
tr
plot(tr)


data(GlaucomaM, package = "ipred")

tr <- ctree(Class ~ ., data = GlaucomaM, teststattype = "quadform")
tr
plot(tr)

data(GBSG2, package = "ipred")  

GBSG2tree <- ctree(Surv(time, cens) ~ ., data = GBSG2, teststattype = "quadform")
GBSG2tree
plot(GBSG2tree)
plot(GBSG2tree, terminal = survNode(GBSG2tree))
survfit(Surv(time, cens) ~ as.factor(GBSG2tree@where), data = GBSG2)
names(GBSG2)

tr <- ctree(Surv(time, cens) ~ ., data = GBSG2, teststattype = "maxabs", testtype = "Raw")
tr
plot(tr)

data(mammoexp)
attr(mammoexp$ME, "scores") <- 1:3   
attr(mammoexp$SYMPT, "scores") <- 1:4
attr(mammoexp$DECT, "scores") <- 1:3 
names(mammoexp)[names(mammoexp) == "SYMPT"] <- "symptoms"
names(mammoexp)[names(mammoexp) == "PB"] <- "benefit"

names(mammoexp)
tr <- ctree(ME ~ ., data = mammoexp, teststattype = "quadform")
tr
plot(tr)

