
setClass("Variable", representation = representation(
    name = "character",
    columns = "integer",
    whichNA = "INTEGERorNULL"))

setClass("OrderedVariable", representation = representation(
    order = "integer"), contains = "Variable")

setClass("CategoricalVariable", representation = representation(
    levels = "character"), contains = "Variable")

setClass("Split", representation = representation(
    variable = "integer", ### should be: a function
    totheleft = "logical",
    criterium = "numeric"))

setClass("OrderedSplit", representation = representation(
    cutpoint = "numeric"), contains = "Split")

setClass("CategoricalSplit", representation = representation(
    levelset = "numeric"), contains = "Split")

setClass("TerminalNode", representation = representation(
    number = "numeric",
    weights = "vector"))

setClass("Node", representation = representation(
    primarysplit = "Split",
    criterion = "vector",  
    surrogatesplits = "list"), contains = "TerminalNode")
  

setClass("GrowControl", representation = representation(
    minsplit = "numeric",
    minprob = "numeric",
    minstat = "numeric",
    varnull = "numeric"))
    
setClass("TreeGrow", representation = representation(
    inputs = "list",
    response = "Variable",
    p = "numeric",
    Scores = "matrix",
    Weights = "matrix",
    control = "GrowControl"))

setClass("BinaryTree", representation = representation(
    nodes = "list",
    nodeindex = "matrix",
    treegrow = "TreeGrow"))



