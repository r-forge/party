perfplot <- function(x, file, boxplot = FALSE, alpha = NULL, border = 1,
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
