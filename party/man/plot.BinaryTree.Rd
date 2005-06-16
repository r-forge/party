\name{plot.BinaryTree}
\alias{plot.BinaryTree}
\alias{innerNode}
\alias{survNode}
\alias{barNode}
\alias{boxNode}
\alias{histNode}
\alias{densNode}
\alias{termNode}
\alias{simpleEdge}
\title{ Visualization of Binary Regression Trees }
\description{
  This method plots objects of class \code{BinaryTree}.
}
\usage{
\method{plot}{BinaryTree}(x, main = NULL,
     type = c("extended", "simple"), terminal.panel = NULL,
     inner.panel = innerNode, edge.panel = simpleEdge,
     drop.terminal = (type[1] == "extended"), 
     tnex = (type[1] == "extended") + 1, new = TRUE,  
     pop = FALSE, ...)
innerNode(ctreeobj, digits = 3, abbreviate = FALSE, 
          fill = "white", id = TRUE)
survNode(ctreeobj, ylines = 2, id = TRUE, ...)
barNode(ctreeobj, col = "black", fill = "lightgray", ymax = NULL,
        ylines = 3, widths = 1, gap = NULL, id = TRUE)
boxNode(ctreeobj, col = "black", fill = "lightgray",  width = 0.5,
        yscale = NULL, ylines = 3, id = TRUE)
histNode(ctreeobj, col = "black", fill = "lightgray", freq = FALSE,
         horizontal = TRUE, xscale = NULL, ymax = NULL, ylines = 3,
         id = TRUE)
densNode(ctreeobj, col = "black", rug = TRUE, horizontal = TRUE,
         xscale = NULL, yscale = NULL, ylines = 3, id = TRUE)
termNode(ctreeobj, digits = 3, abbreviate = FALSE, 
         fill = c("lightgray", "white"), id = TRUE)
simpleEdge(ctreeobj, digits = 3, abbreviate = FALSE)

}
\arguments{
  \item{x}{ an object of class \code{BinaryTree}.}
  \item{main}{ an optional title for the plot.}
  \item{type}{ a character specifying the complexity of the plot:
               \code{extended} tries to visualize the distribution of the
               response variable in each terminal node whereas \code{simple} 
               only gives some summary information.}
  \item{terminal.panel}{ an optional panel function of the form 
                         \code{function(node)} plotting the terminal nodes.
                         Alternatively, a panel generating function of the form
                         \code{function(ctreeobj, ...)} can be supplied.}
  \item{inner.panel}{ an optional panel function of the form 
                         \code{function(node)} plotting the inner nodes.}
  \item{edge.panel}{ an optional panel function of the form 
                    \code{function(split, ordered = FALSE, left = TRUE)}
                    plotting the edges.}
  \item{drop.terminal}{ a logical indicating whether all terminal nodes
                        should be plotted at the bottom.}
  \item{tnex}{ a numeric value giving the terminal node extention. }
  \item{new}{ a logical indicating whether a new plot should be created. }
  \item{pop}{ a logical whether the viewport tree should be popped before
              return. }
  \item{...}{ additional arguments passed to callies.}
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
  \item{gap}{ gap between bars in a barplot (\code{barNode}). }
  \item{yscale}{ limits in y-direction}
  \item{xscale}{ limits in x-direction}
  \item{ymax}{ upper limit in y-direction}
  \item{horizontal}{ logical indicating if the plots should be horizontal. }
  \item{freq}{logical; if \code{TRUE}, the histogram graphic is a representation
              of frequencies. If \code{FALSE}, probabilities are plotted.}
  \item{rug}{logical indicating if a rug representation should be added. }
}
\details{

  This \code{plot} method for \code{BinaryTree} objects provides an
  extensible framework for the visualization of binary regression trees. The
  user is allowed to specify panel functions for plotting terminal and inner
  nodes as well as the corresponding edges. For convenience, the panel
  generating functions \code{innerNode} and \code{simpleEdge} return 
  panel functions to draw inner nodes 
  and left and right edges. For drawing terminal nodes, the functions
  returned by the other panel functions can be used. 

  \code{survNode} returns a functions that plots Kaplan-Meier curves in each 
  terminal node. \code{barNode, boxNode, histNode} and \code{densNode} can
  be used to plot barplots, boxplots, histograms and density estimators into 
  the terminal nodes.
}
\examples{

  data(airquality)
  airq <- subset(airquality, !is.na(Ozone))
  airct <- ctree(Ozone ~ ., data = airq)

  ### very simple; the mean is given in each terminal node
  plot(airct, type = "simple")

  ### density estimates
  plot(airct, terminal.panel = densNode)
    
  ### histograms 
  plot(airct, terminal.panel = histNode)

  ### histograms
  plot(airct, terminal.panel = boxNode, drop.terminal = TRUE)

}
\keyword{hplot}
