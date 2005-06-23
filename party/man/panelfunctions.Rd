\name{Panel Functions}
\alias{panel_inner_node}
\alias{panel_terminal_node}
\alias{panel_simple_edge}
\alias{panel_surv}
\alias{panel_bar}
\alias{panel_box}
\alias{panel_hist}
\alias{panel_dens}
\title{ Panel Functions for the Visualization of Binary Regression Trees }
\description{
    The plot method for \code{BinaryTree} objects is rather flexible and can
    be extended by panel functions. Some pre-defined panel functions for the
    most important cases are documented here.
}
\usage{
panel_inner_node(ctreeobj, digits = 3, abbreviate = FALSE, 
                 fill = "white", id = TRUE)
panel_terminal_node(ctreeobj, digits = 3, abbreviate = FALSE, 
         fill = c("lightgray", "white"), id = TRUE)
panel_simple_edge(ctreeobj, digits = 3, abbreviate = FALSE)
panel_surv(ctreeobj, ylines = 2, id = TRUE, \dots)
panel_bar(ctreeobj, col = "black", fill = "lightgray", ymax = NULL,
        ylines = 3, widths = 1, gap = NULL, id = TRUE)
panel_box(ctreeobj, col = "black", fill = "lightgray",  width = 0.5,
        yscale = NULL, ylines = 3, id = TRUE)
panel_hist(ctreeobj, col = "black", fill = "lightgray", freq = FALSE,
         horizontal = TRUE, xscale = NULL, ymax = NULL, ylines = 3,
         id = TRUE, \dots)
panel_dens(ctreeobj, col = "black", rug = TRUE, horizontal = TRUE,
         xscale = NULL, yscale = NULL, ylines = 3, id = TRUE)
}
\arguments{
  \item{ctreeobj}{ an object of class \code{BinaryTree}.}
  \item{digits}{ integer, used for formating numbers. }
  \item{abbreviate}{ logical indicating whether strings should be 
                     abbreviated. }
  \item{col}{ a color for points and lines. }
  \item{fill}{ a color to filling rectangles. }
  \item{id}{ a logical indicating whether the node ID should be plotted.}
  \item{ylines}{ number of lines in y-direction. }
  \item{widths}{ widths in barplots. }
  \item{width}{ width in boxplots. }
  \item{gap}{ gap between bars in a barplot (\code{panel_bar}). }
  \item{yscale}{ limits in y-direction}
  \item{xscale}{ limits in x-direction}
  \item{ymax}{ upper limit in y-direction}
  \item{horizontal}{ logical indicating if the plots should be horizontal. }
  \item{freq}{logical; if \code{TRUE}, the histogram graphic is a representation
              of frequencies. If \code{FALSE}, probabilities are plotted.}
  \item{rug}{logical indicating if a rug representation should be added. }
  \item{\dots}{ additional arguments passed to callies.}
}
\details{

  This \code{plot} method for \code{BinaryTree} objects provides an
  extensible framework for the visualization of binary regression trees. The
  user is allowed to specify panel functions for plotting terminal and inner
  nodes as well as the corresponding edges. For convenience, the panel
  generating functions \code{panel_inner_node} and \code{panel_simple_edge} 
  return panel functions to draw inner nodes and left and right edges. 
  For drawing terminal nodes, the functions returned by the other panel 
  functions can be used. The panel function \code{panel_terminal_node} is
  a terse text-based representation of terminal nodes.

  Graphical representations of terminal nodes are available and depend on
  the measurement scale of the response variable the tree was fitted to.
  \code{panel_surv} returns a functions that plots Kaplan-Meier curves in each 
  terminal node. \code{panel_bar}, \code{panel_box}, \code{panel_hist} and
  \code{panel_dens} can be used to plot barplots, boxplots, histograms and 
  estimated densities into the terminal nodes.
}
\examples{

  data(airquality)
  airq <- subset(airquality, !is.na(Ozone))
  airct <- ctree(Ozone ~ ., data = airq)

  ### very simple; the mean is given in each terminal node
  plot(airct, type = "simple")

  ### density estimates
  plot(airct, terminal_panel = panel_dens)
    
  ### histograms 
  plot(airct, terminal_panel = panel_hist)

  ### boxplots
  plot(airct, terminal_panel = panel_box, drop_terminal = TRUE)

}
\keyword{hplot}
