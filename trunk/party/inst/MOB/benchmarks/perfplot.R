pcbplot <- function(x, file, boxplot = TRUE, alpha = NULL, border = 1,
  height = 5, width = 7, ylab = "Performance") {
  x <- as.matrix(x)
  nc <- NCOL(x)
  nr <- NROW(x)
  nam <- colnames(x)
  nam[nam == "party"] <- "ctree"
  if(is.null(alpha)) alpha <- 0.17 - 0.0002 * nr

  pdf <- !missing(file)

  if(pdf) pdf(file = file, height = height, width = width, version = "1.4")
  plot(rep(1:nc, 2), rep(range(x), nc), type = "n", axes = FALSE,
       ylab = ylab, xlab = "", xlim = c(0.6, nc + 0.4))
  axis(1, at = 1:nc, labels = nam)
  axis(2)
  box()

  if(!boxplot) abline(v = 1:4, col = grey(0.7))

  matlines(t(x), col = rgb(0,0,0, alpha), lty = 1)

  if(boxplot) boxplot(as.vector(x) ~ factor(rep(nam, rep(nr, nc)), levels = nam),
          add = TRUE, xlab = "", ylab = "", axes = FALSE, border = border)

  if(pdf) dev.off()
}

ciplot <- function(x, file, height = 5, width = 7, xlab = "Performance difference") {
  tmp <- x - apply(x, 1, mean)
  tmp <- data.frame(error = unlist(tmp), model = factor(rep(colnames(tmp),
    rep(nrow(tmp), ncol(tmp))), levels = names(x)))
  si <- simint(error ~ model, data = tmp, type = "Dunnett")
  rownames(si$estimate) <- names(x)[-1]

  pdf <- !missing(file)
  if(pdf) pdf(file = file, height = height, width = width, version = "1.4")
  plot(si, main = "", xlab = xlab)
  if(pdf) dev.off()
  invisible(si)
}

perfplot <- function(x, file, height = 5, width = 11, lab = "Performance", ...) {
  pdf <- !missing(file)
  if(pdf) pdf(file = file, height = height, width = width, version = "1.4")
  par(mfrow = c(1, 2))
  pcbplot(x, ylab = lab, ...)
  ciplot(x, xlab = paste(lab, "difference"))
  if(pdf) dev.off()
}
