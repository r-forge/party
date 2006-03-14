
R2guide <- function(x, response, name) {

    dscname <- paste(name, ".dsc", sep = "")
    datname <- paste(name, ".dat", sep = "")

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

    write.table(x, sep = " ", file = datname, row.names = FALSE, col.names = FALSE)
}

foo <- function(x, response, bs) {

    err <- rep(0, ncol(bs))
    for (i in 1:ncol(bs)) {

        learn <- x
        learn$weight <- bs[,i]

        R2guide(learn, response, "learn")

        a <- system("./guide_linux < reg.run", intern = TRUE)
        pred <- read.table("predict.txt", header = TRUE)$fitted

        diff <- learn[[response]] - pred
        err[i] <- (mean(diff[bs[,i] == 0]^2))
        cat(i, " ", err[i], "\n")

        system("rm results.out")
        system("rm predict.txt")

    }
    err
}

