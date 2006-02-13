\name{mob}
\encoding{latin1}

\alias{mob}
\alias{coef.mob}
\alias{fitted.mob}
\alias{logLik.mob}
\alias{plot.mob}
\alias{predict.mob}
\alias{print.mob}
\alias{residuals.mob}
\alias{sctest.mob}
\alias{summary.mob}

\title{Model-based Recursive Partitioning}

\description{
  MOB is an algorithm for model-based recursive partitioning yielding
  a tree with fitted models associated with each terminal node.
}

\usage{
mob(formula, weights, data = list(), model = glinearModel,
  control = mob_control(), \dots)

\method{predict}{mob}(object, newdata = NULL, type = c("response", "node"), \dots)
\method{summary}{mob}(object, node = NULL, \dots)
\method{coef}{mob}(object, node = NULL, \dots)
\method{sctest}{mob}(object, node = NULL, \dots)
}

\arguments{
  \item{formula}{A symbolic description of the model to be fit. This
    should be of type \code{y ~ x1 + \dots + xk | z1 + \dots + zl} where
    the variables before the \code{|} are passed to the \code{model} and
    the variables after the \code{|} are used for partitioning.}
  \item{weights}{An optional vector of weights to be used in the fitting
    process. Only non-negative integer valued weights are allowed (default = 1).}
  \item{data}{A data frame containing the variables in the model.}
  \item{model}{A model of class \code{\link{StatModel-class}}. See details
    for requirements.}
  \item{control}{A list with control parameters as returned by
    \code{\link{mob_control}}.}
  \item{\dots}{Additional arguments passed to the \code{\link{fit}} call for
    the \code{model}.}
  \item{object}{A fitted \code{mob} object.}
  \item{newdata}{A data frame with new inputs, by default the learning data
    is used.}
  \item{type}{A character string specifying whether the response should be
    predicted (inherited from the \code{predict} method for the \code{model})
    or the ID of the associated terminal node.}
  \item{node}{A vector of node IDs for which the corresponding method should
    be applied.}
}

\details{
Model-based partitioning fits a model tree using the following algorithm:
  \enumerate{
    \item \code{fit} a \code{model} (default: a generalized linear model
          \code{\link{glinearModel}} with formula \code{y ~ x1 + \dots + xk}
	  for the observations in the current node.
    \item Assess the stability of the model parameters with respect to each
          of the partitioning variables \code{z1}, \dots, \code{zl}. If
	  there is some overall instability, choose the variable \code{z}
	  associated with the smallest \eqn{p} value for partitioning, otherwise
	  stop. For performing the parameter instability fluctuation test,
	  a \code{\link{estfun}} method and a \code{\link{weights}} method is
	  needed.
    \item Search for the locally optimal split in \code{z} by minimizing the
          objective function of the \code{model}. Typically, this will be
	  something like \code{\link{deviance}} or the negative \code{\link{logLik}}
	  and can be specified in \code{\link{mob_control}}.
    \item Re-fit the \code{model} in both children, using \code{\link{reweight}}
          and repeat from step 2.
  }
  
More details on the conceptual design of the algorithm can be found in 
Zeileis, Hothorn, Hornik (2005) and some illustrations are provided in
\code{vignette("MOB")}.  
  
For the fitted MOB tree, several standard methods are inherited if they are
available for fitted \code{model}s, such as \code{print}, \code{predict},
\code{residuals}, \code{logLik}, \code{coef} and \code{summary}. By default, the latter
two return the coefficients and summary for all terminal nodes,
but take a \code{node} argument that can be set to any node ID. The \code{sctest}
method extracts the results of the parameter stability tests (aka structural 
change tests) for any given node, by default for all nodes. Some examples are
given below.
}

\value{
  An object of class \code{mob} inheriting from \code{\link{BinaryTree-class}}.
}

\references{ 

   Achim Zeileis, Torsten Hothorn and Kurt Hornik (2005). Model-based
   Recursive Partitioning. \emph{Report 19}, Department of Statistics and Mathematics, 
   \enc{Wirtschaftsuniversität}{Wirtschaftsuniversitaet} Wien, Research Report Series.
   \url{http://epub.wu-wien.ac.at/dyn/openURL?id=oai:epub-wu-01_86e}

}

\seealso{\code{\link{mob_control}}}

\examples{

if(require("mlbench")) {

## recursive partitioning of a linear regression model
## load data
data("BostonHousing", package = "mlbench")
## and transform variables appropriately (for a linear regression)
BostonHousing$lstat <- log(BostonHousing$lstat)
BostonHousing$rm <- BostonHousing$rm^2
## as well as partitioning variables (for fluctuation testing)
BostonHousing$chas <- factor(BostonHousing$chas, levels = 0:1, labels = c("no", "yes"))
BostonHousing$rad <- factor(BostonHousing$rad, ordered = TRUE)

## partition the linear regression model medv ~ lstat + rm
## with respect to all remaining variables:
fm <- mob(medv ~ lstat + rm | zn + indus + chas + nox + age + dis + rad + tax + crim + b + ptratio,
  control = mob_control(minsplit = 40), data = BostonHousing, model = linearModel)

## print the resulting tree
fm
## or better visualize it
plot(fm, tnex = 4)

## extract coefficients in all terminal nodes
coef(fm)
## look at full summary, e.g., for node 7
summary(fm, node = 7)
## results of parameter stability tests for that node
sctest(fm, node = 7)
## -> no further significant instabilities (at 5\% level)

## compute mean squared error (on training data)
mean((BostonHousing$medv - fitted(fm))^2)
mean(residuals(fm)^2)

## evaluate AIC
AIC(fm)
## (Note that this penalizes estimation of error variances, which
## were treated as nuisance parameters in the fitting process.)

}
}
\keyword{tree}
