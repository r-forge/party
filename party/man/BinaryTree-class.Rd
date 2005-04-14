\name{BinaryTree-class}
\docType{class}
\alias{BinaryTree-class}

\title{Class "BinaryTree"}
\description{A class for representing binary trees and some information
             about the data used for fitting.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("BinaryTree", ...)}.
The most important slot is \code{tree}, a list with components
\describe{
  \item{nodeID}{ an integer giving the number of the node,}
  \item{weights}{ the case weights (of the learning sample) corresponding to
                this node,}
  \item{criterion}{ a list with test statistics and P-values for each partial
                  hypothesis,}
  \item{terminal}{ a logical specifying if this is a terminal node,}
  \item{psplit}{ primary split: a list with elements \code{variableID} (the
               number of the input variable splitted), \code{ordered} (a
               logical whether the input variable is ordered),
               \code{splitpoint} (the cutpoint or set of levels to the left)
               and \code{splitstatistics} (some details about the selection),}
  \item{ssplits}{ surrogate splits (not implemented, yet),}
  \item{prediction}{ the prediction of this node,}
  \item{left}{ the left node and }
  \item{right}{ the right node.}
}

Please note that this data structure may be subject to change in future
releases of the package.

}
\section{Slots}{
  \describe{
    \item{\code{inputnames}:}{ a character vector with the names of
                               the input variables used for fitting the tree.}
    \item{\code{responses}:}{ an object of class \code{"VariableFrame"}
                              storing the values of the response variable(s). }
    \item{\code{levels}:}{ a list of character vectors storing the levels of
                           nominal or ordered input variables. }
    \item{\code{cond_distr_response}:}{ a function computing the conditional
                                        distribution of the response. } 
    \item{\code{predict_response}:}{ a function for computing predictions. }
    \item{\code{tree}:}{ a recursive list representing the tree. See above. }
    \item{\code{where}:}{ an integer vector of length n (number of
                          observations in the learning sample) giving the
                          number of the terminal node the corresponding
                          observations is element of. }
  }
}
\section{Extends}{
Class \code{"BinaryTreePartition"}, directly.
}
\section{Methods}{
A \code{print}, \code{plot} and \code{predict} method are available. 
}
\keyword{classes}
