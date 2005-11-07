
library("Rweka")
data(iris)

j48tree <- j48(Species ~ ., data = iris)
a <- as.BinaryTree(j48tree)
a@data <- ModelEnvFormula(Species ~ ., data = iris)


### `new' input variables
inp <- party:::initVariableFrame(j48tree@data, trafo = NULL)

### predictions
.Call("R_predict", a@tree, inp, 0.0, PACKAGE = "party")
