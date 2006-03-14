npar <- function(obj) {  

  switch(class(obj)[1],

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
  })

  ## number of coefficients + number of splits
  rval <- nc + (nl - 1)
  return(rval)
}
