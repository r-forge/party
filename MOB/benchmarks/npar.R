npar <- function(obj) {  

  switch(class(obj)[1],

  "lm" = {
    ## number of leaves
    nl <- 1
    ## number of coefficients
    nc <- length(coef(obj))
  },
  
  "glm" = {
    ## number of leaves
    nl <- 1
    ## number of coefficients
    nc <- length(coef(obj))
  },
  
  "mob" = {
    ## number of leaves
    nl <- party:::nterminal(obj@tree)
    ## number of coefficients
    nc <- prod(dim(coef(obj)))
  },
  
  "BinaryTree" = {
    ## number of leaves
    nl <- party:::nterminal(obj@tree)
    ## number of coefficients
    nc <- nl
  },

  "rpart" = {
    nl <- length(unique(obj$where))
    nc <- nl
  },
  
  "J48" = {
    ## number of leaves
    nl <- strsplit(.jcall(obj$classifier, "S", "toString"), "\n")[[1]]
    nl <- nl[grep("Number of Leaves", nl)]
    nl <- as.numeric(strsplit(nl, "\t")[[1]][2]) 
    ## number of coefficients
    nc <- nl
  },
  
  "LMT" = {
    ## number of leaves
    nl <- strsplit(.jcall(obj$classifier, "S", "toString"), "\n")[[1]]
    nl <- nl[grep("Number of Leaves", nl)]
    nl <- as.numeric(strsplit(nl, "\t")[[1]][2]) 
    ## number of coefficients
    nc <- strsplit(.jcall(obj$classifier, "S", "toString"), "Number of Leaves")[[1]][2]
    nc <- strsplit(nc, "LM_", extended = FALSE, fixed = TRUE)[[1]][-1]
    nc <- strsplit(nc, "\n")
    if(length(nc) != nl) warning("something went wrong in computing the number of leaves")
    nc <- sapply(nc, function(x) length(grep("[", x, fixed = TRUE, extended = FALSE))/2 + 1)
    nc <- sum(nc)    
  },
  
  "M5P" = {
    ## number of leaves
    nl <- strsplit(.jcall(obj$classifier, "S", "toString"), "\n")[[1]]
    nl <- nl[grep("Number of Rules", nl)]
    nl <- as.numeric(strsplit(nl, "Number of Rules : ")[[1]][2])
    ## number of coefficients
    nc <- .jcall(obj$classifier, "S", "toString")
    nc <- strsplit(nc, "LM num:", extended = FALSE, fixed = TRUE)[[1]][-1]
    nc <- strsplit(nc, "\n")
    if(length(nc) != nl) warning("something went wrong in computing the number of leaves")
    nc <- sapply(nc, function(x) length(grep("\t", x, fixed = TRUE, extended = FALSE)))
    nc <- sum(nc)    
  },

  "GUIDE" = {
    ## read output
    nl <- nc <- obj$tree
    ## number of leaves
    nl <- nl[grep("Number of terminal nodes of final tree", nl)]
    nl <- as.numeric(strsplit(nl, ": ")[[1]][2]) 
    ## number of parameters
    if(nl > 1) {
      nc <- diff(tail(which(nc == " ----------------------------"), nl+1)) - 5
    } else {
      nc <- which(nc == " ----------------------------") - which(nc == " Node 1: Terminal node") - 4
    }
    if(length(nc) != nl) warning("something went wrong in computing the number of leaves")
    nc <- sum(nc)
  },
  
  "CRUISE" = {
    nl <- obj$tree[grep("Number of terminal nodes in final tree", obj$tree)]
    nl <- as.numeric(strsplit(nl, " = ")[[1]][2])
    nc <- nl  
  },
  
  "QUEST" = {
    nl <- obj$tree
    nl <- nl[grep("Number of terminal nodes of final tree", nl)]
    nl <- as.numeric(strsplit(nl, " = ")[[1]][2]) 
    nc <- nl
  })
  
  ## number of coefficients + number of splits
  rval <- nc + (nl - 1)
  return(rval)
}
