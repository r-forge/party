
R2guide <- function(x, response, name) {

    datname <- paste(name, ".dat", sep = "")
    write.table(x, sep = " ", file = datname, row.names = FALSE, col.names = FALSE)

    dscname <- paste(name, ".dsc", sep = "")

    cat(datname, "\n", file = dscname, append = FALSE)
    cat("NA\n", file = dscname, append = TRUE)
    cat("column, varname, vartype\n", file = dscname, append = TRUE)

    for (n in 1:ncol(x)) {

        cat(n, names(x)[n], " ", file = dscname, append = TRUE)
  
        if (names(x)[n] == response) {
            cat("d\n", file = dscname, append = TRUE)
        } else {
            if (names(x)[n] == "weight") {
                cat("w\n", file = dscname, append = TRUE)
            } else {
                cat(ifelse(is.factor(x[[n]]), "c", "n"), "\n",
                    file = dscname, append = TRUE)
            }
        }

    }
}

complexity_guide <- function(file = "results.out") {
  ## read output
  nl <- nc <- readLines(file)
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
  rval <- c(nl = nl, nc = nc)
  return(rval)  
}

foo <- function(x, response, bs, obvious = FALSE) {

    err <- rep(0, ncol(bs))
    npar <- rep(0, ncol(bs))
    for (i in 1:ncol(bs)) {

        learn <- x
        learn$weight <- bs[,i]

        R2guide(learn, response, "learn")

        a <- system("./guide_linux < reg.run", intern = TRUE)
        pred <- read.table("predict.txt", header = TRUE)$fitted

        diff <- learn[[response]] - pred
        err[i] <- if(obvious) (mean(diff[bs[,i] > 0]^2)) else (mean(diff[bs[,i] == 0]^2))
        npar[i] <- sum(complexity_guide()) - 1
        cat(i, " ", err[i], " ", npar[i], "\n")

        system("rm results.out")
        system("rm predict.txt")

    }
    rval <- list(error = err, npar = npar)
    return(rval)
}

