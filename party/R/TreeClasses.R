
setClass("Variable", representation = representation(
    name = "character",
    columns = "integer",
    whichNA = "INTEGERorNULL"))

setClass("OrderedVariable", representation = representation(
    order = "integer"), contains = "Variable")

setClass("CategoricalVariable", representation = representation(
    levels = "character"), contains = "Variable")


setClass("Split", representation = representation(
    variable = "integer",
    cutpoint = "numeric",
    totheleft = "logical"))

setClass("TerminalNode", representation = representation(
    number = "numeric",
    weights = "vector"))

setClass("Node", representation = representation(
    primarysplit = "Split",
    criterion = "vector",  
    surrogatesplits = "list"), contains = "TerminalNode")
  

setClass("GrowControl", representation = representation(
    minsplit = "numeric",
    minstat = "numeric"))

    
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



