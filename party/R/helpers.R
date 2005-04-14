
# $Id$

dostep <- function(x, y) {

    ### create a step function based on x, y coordinates
    ### modified from `survival:print.survfit'
    if (is.na(x[1] + y[1])) {
        x <- x[-1]
        y <- y[-1]
    }
    n <- length(x)
    if (n > 2) {  
        # replace verbose horizonal sequences like
        # (1, .2), (1.4, .2), (1.8, .2), (2.3, .2), (2.9, .2), (3, .1)
        # with (1, .2), (3, .1).  They are slow, and can smear the looks
        # of the line type.
        dupy <- c(TRUE, diff(y[-n]) !=0, TRUE)
        n2 <- sum(dupy)

        #create a step function
        xrep <- rep(x[dupy], c(1, rep(2, n2-1)))
        yrep <- rep(y[dupy], c(rep(2, n2-1), 1))
        RET <- list(x = xrep, y = yrep)
    } else {
        if (n == 1) {
            RET <- list(x = x, y = y)
        } else {
            RET <- list(x = x[c(1,2,2)], y = y[c(1,1,2)])
        }
    }
    return(RET)
}

### Normalized mutual information, Stehl & Gosh (2002) JMLR
nmi <- function(x, y)
{
    x <- table(x, y)
    x <- x / sum(x)                     # ???

    m_x <- rowSums(x)
    m_y <- colSums(x)
    y <- outer(m_x, m_y)

    i <- which((x > 0) & (y > 0))
    out <- sum(x[i] * log(x[i] / y[i]))
    e_x <- sum(m_x * log(ifelse(m_x > 0, m_x, 1)))
    e_y <- sum(m_y * log(ifelse(m_y > 0, m_y, 1)))

    out / sqrt(e_x * e_y)
}

### Median Survival Time, see survival:::print.survfit

mst <- function(x) {
    minmin <- function(y, xx) {

        ### don't complain about Inf 
        ww <- getOption("warn")
        on.exit(options(warn = ww))
        options(warn = -1)

        if (any(!is.na(y) & y==.5)) {
            if (any(!is.na(y) & y <.5))
                .5*(min(xx[!is.na(y) & y==.5]) + min(xx[!is.na(y) & y<.5]))
            else
                .5*(min(xx[!is.na(y) & y==.5]) + max(xx[!is.na(y) & y==.5]))
        } else   min(xx[!is.na(y) & y<=.5])
    }
    med <- minmin(x$surv, x$time)
    return(med)
}

### Logrank scores, see exactRankTests:::cscores.Surv

logrank_scores <- function(x) {

  if (class(x) != "Surv")
      stop(sQuote("x"), " is not of class ", sQuote("Surv"))

  time <- x[,1]
  event <- x[,2]
  N <- length(time)
  ot <- order(time)
  rt <- rank(time, ties.method = "max")
  fact <- event/(N - rt + 1)
  RET <- event - cumsum(fact[ot])[rt]
  attr(RET, "scores") <- "LogRank"
  return(RET)
}
