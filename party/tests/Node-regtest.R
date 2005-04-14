
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
node <- .Call("R_Node", ls, ls@weights, tm, ctrl)
stopifnot(isequal(node[[5]][[3]], 0.059))

### and now with ranked inputs -> Wilcoxon-Mann-Whitney tests
inp <- party:::initVariableFrame(GlaucomaM[,-63,drop = FALSE], fun = rank)
resp <- party:::initVariableFrame(GlaucomaM[,"Class",drop = FALSE])
ls <- new("LearningSample", inputs = inp, responses = resp,
          weights = rep(1, inp@nobs), nobs = nrow(GlaucomaM), 
          ninputs = inp@ninputs)
tm <- TreeFitMemory(ls, TRUE)
ctrl <- new("TreeControl")
ctrl@gtctrl@testtype <- factor("Bonferroni", levels = tlev)
ctrl@varctrl@teststattype <- factor("quadform", levels = c("maxabs", "quadform"))
node <- .Call("R_Node", ls, ls@weights, tm, ctrl)
stopifnot(isequal(node[[5]][[3]], 0.059))
