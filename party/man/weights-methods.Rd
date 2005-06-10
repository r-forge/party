\name{weights-methods}
\docType{methods}
\alias{weights}
\alias{weights-methods}
\alias{weights,BinaryTree-method}
\alias{where}
\alias{where-methods}
\alias{where,BinaryTree-method}
\alias{nodes}
\alias{nodes-methods}
\alias{nodes,BinaryTree,integer-method}
\alias{nodes,BinaryTree,numeric-method}
\title{ Nodes of BinaryTree Objects }
\description{

    Extract node numbers, the nodes itself or weights associated with
    nodes of BinaryTree objects, for the learning sample or new test
    samples.

}
\section{Methods}{
\describe{

\item{object = "BinaryTree"}{ the fitted tree }

\item{object = "BinaryTree", where = "integer"}{ a vector of node numbers. }

\item{object = "BinaryTree", where = "numeric"}{ a vector of node numbers. }


}}
\section{Details}{

    The \code{where} method determines the number of the terminal node for
    each observation in the learning sample or, if \code{newdata} is
    supplied, for each observation in a test sample.

    The \code{weights} method returns the vector of case weights for each
    observation in the learning sample or, if \code{newdata} is
    supplied, for each observation in a test sample.
 
    The \code{nodes} method can be used to extract nodes (terminal or inner
    nodes) from a tree \code{object}.

} 
\examples{

  data(airquality)
  airq <- subset(airquality, !is.na(Ozone))
  airct <- ctree(Ozone ~ ., data = airq, 
                 controls = ctree.control(maxsurrogate = 3))

  ### distribution of responses in the terminal nodes
  plot(airq$Ozone ~ as.factor(where(airct)))

  ### get all terminal nodes from the tree
  nodes(airct, unique(where(airct)))

  ### extract weights and compute predictions
  pmean <- sapply(weights(airct), function(w) weighted.mean(airq$Ozone, w))

  ### the same as 
  drop(Predict(airct))

  ### don't use the mean but the median as prediction in each terminal node
  pmedian <- sapply(weights(airct), function(w) 
                                        median(airq$Ozone[rep(1:nrow(airq), w)]))

  plot(airq$Ozone, pmean, col = "red")
  points(airq$Ozone, pmedian, col = "blue")

}
\keyword{methods}
