

###### Classes for dealing with variables ######


### A general class for univariate variables
setClass(Class = "Variable", 
    representation = representation(
        name = "character",                # variable name 
        values = "matrix",                 # numeric representation:
                                           # either original values
                                           # or dummy encoded
        whichNA = "INTEGERorNULL"          # position of NA elements
    )
)

### A continuous variable
setClass(Class = "ContinuousVariable", 
    representation = representation(
        order = "integer"                  # the ordering: once and for all
    ), 
    contains = "Variable"
)

### A categorical variable: unordered (corresponds to `factor' objects)
setClass(Class = "CategoricalVariable", 
    representation = representation(
        levels = "character",              # the set of levels and
        coding = "integer"                 # the levels of the observations
    ), 
    contains = "Variable"
)

### An ordered categorical varibale (corresponds to `ordered' objects)
setClass(Class = "OrderedCategoricalVariable", 
    representation = representation(
        scores = "numeric",                # numeric scores giving the 
                                           # "distances"
        order = "integer"                  # the ordering
    ), 
    contains = "CategoricalVariable"
)


###### Classes for dealing with splits ######


### A general class for splits
setClass(Class = "Split", 
    representation = representation(
        variable = "integer",              # number of a certain variable
                                           # should be: a function of all
                                           # variables returning a scalar
                                           # value
        totheleft = "logical",             # condition sends observations 
                                           # to the left or right
        criterium = "numeric"              # a measure for goodness of
                                           # separation
    )
)

### Split in a continuous variable
setClass(Class = "ContinuousSplit", 
    representation = representation(
        cutpoint = "numeric"               # condition: x <= cutpoint 
    ), 
    contains = "Split"
)

### Split in a ordered categorical variable
setClass(Class = "OrderedSplit", 
    representation = representation(
        levelset = "numeric"
    ), 
    contains = "Split"
)

### Split in a categorical variable
setClass(Class = "CategoricalSplit", 
    representation = representation(
        levelset = "numeric"               # condition: x %in% levelset
    ), 
    contains = "Split"
)


###### Classes for dealing with binary trees ######


### A class for terminal nodes
setClass(Class = "TerminalNode", 
    representation = representation(
        number = "numeric",                # a unique number identifying
                                           # this node
        weights = "vector"                 # case weights in this node
    )
)

### A class for non-terminal nodes
setClass(Class = "Node", 
    representation = representation(
        primarysplit = "Split",            # the condition for splitting
        criterion = "vector",              # a measure for goodness of 
                                           # selection
        surrogatesplits = "list"           # a list of `Split' objects
                                           # considered when NA's occur
    ), 
    contains = "TerminalNode"
)
  
### A class for some parameters controlling the tree growing
setClass(Class = "PartyControl", 
    representation = representation(
        minsplit = "numeric",              # try to split when > minsplit
                                           # obs. are available only
        minprob = "numeric",               # at least minprob*100% 
                                           # obs. in each branch of a node
        minstat = "numeric",               # stop splitting if the best 
                                           # test statistic is < minstat
        varnull = "numeric"                # variances are interpreted as
                                           # zero is < varnull
    ),
    prototype = list(minsplit = 20, 
                     minprob = 0.1, 
                     minstat = 1.96, 
                     varnull = 1e-10)
)

### A class for handling all ingredients needed to start the party    
setClass(Class = "PartyStarters", 
    representation = representation(
        inputs = "list",                   # a list of input `Variables' 
        p = "numeric",                     # length(inputs)
        response = "Variable",             # a single response `Variable'
        workingresponse = "Variable",      # a single (transformed) response
                                           # `Variable'
        nobs = "numeric",                  # number of observations
        control = "PartyControl"           # some parameters
    ), 
)

### A class for binary trees
setClass(Class = "BinaryTree", 
    representation = representation(
        nodes = "list",                    # a list of objects inheriting
                                           # from class `TerminalNode'
        nodeindex = "matrix",              # a matrix of "pointers"
        treegrow = "PartyStarters"         # everything used to grow the
                                           # tree
    )
)
