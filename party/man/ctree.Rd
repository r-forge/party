\name{ctree}
\alias{ctree}
\alias{conditionalTree}
\title{ Conditional Trees }
\description{
  Recursive partitioning for continuous, censored, ordered, nominal and
  multivariate response variables in a conditional inference framework. 
}
\usage{
ctree(formula, data, subset = NULL, weights = NULL, 
      teststattype = c("quadform", "maxabs"), 
      testtype = c("Bonferroni", "MonteCarlo", "Raw"), 
      mincriterion = 0.95, minsplit = 20, stump = FALSE, 
      nresample = 9999, maxsurrogate = 0, ifun = NULL, rfun = NULL)
}
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. }
  \item{data}{ an data frame containing the variables in the model. }
  \item{subset}{ an optional vector specifying a subset of observations to be
                 used in the fitting process.}
  \item{weights}{ an optional vector of weights to be used in the fitting
                  process. Only non-negative integer valued weights are
                  allowed.}
  \item{teststattype}{ a character specifying the type of the test statistic 
                       to be applied. }
  \item{testtype}{ a character specifying how to compute the distribution of
                   the test statistic. }
  \item{mincriterion}{ the value of the test statistic or 1 - P-value that
                       must be exceeded in order to implement a split. }
  \item{minsplit}{ the minimum sum of the weights in a node. }
  \item{stump}{ a logical determining whether a stump (a tree with three
                nodes only) is to be computed. }
  \item{nresample}{ number of Monte-Carlo replications to use when the 
                    distribution of the test statistic is simulated.}
  \item{maxsurrogate}{ number of surrogate splits to evaluate. Note the
                       currently only surrogate splits in ordered
                       covariables are implemented. }
  \item{ifun}{ a function to be applied to all input variables. }
  \item{rfun}{ a function to be applied to all response variables. }
}
\details{

  Conditional trees estimate a regression relationship by binary recursive
  partitioning in a conditional inference framework. Roughly, the algorithm
  works as follows: 1) Test the global null hypothesis of independence between
  any of the input variables and the response (which may be multivariate as well). 
  Stop if this hypothesis cannot be rejected. Otherwise select the input
  variable with strongest association to the resonse. This
  association is measured by a P-value corresponding to a test for the
  partial null hypothesis of a single input variable and the response.
  2) Implement a binary split in the selected input variable. 3) Recursively repeate
  steps 1) and 2). 

  The implementation utilizes a unified framework for conditional inference,
  or permutation tests, developed by Strasser and Weber (1999). The stop
  criterion in step 1) is either based on a P-value (\code{teststattype} is
  \code{Bonferroni} or \code{MonteCarlo}) or on the raw (standardized) test
  statistic (\code{teststattype} is \code{Raw}). In both cases, the
  criterion is maximized, i.e., 1 - P-value is used. A split is implemented 
  when the criterion exceeds the value given by \code{mincriterion}. For
  example, when \code{mincriterion = 0.95}, the P-value must be smaller than
  $0.05$ in order to split this node. This statistical approach ensures that
  the right sized tree is grown and no form of pruning or cross-validation
  or whatsoever is needed. The selection of the input variable to split in
  is based on the univariate P-values (when \code{teststattype} is
  \code{Bonferroni} or \code{MonteCarlo}) avoiding a variable selection bias
  towards input variables with many possible cutpoints.

  For a general description of the methodology see Hothorn, Hornik and
  Zeileis (2004).

  The current handling of missing values is very rudimentary. An error is
  thrown when the response has at least one missing value. Observations
  with missing values in an input variable defining a split go with the
  majority.

}
\value{
  An object of class \code{BinaryTree}.
}
\references{ 

   Torsten Hothorn, Kurt Hornik and Achim Zeileis (2004). Unbiased Recursive
   Partitioning: A Conditional Inference Framework. Technical Report Nr. 8, 
   Research Report Series / Department of Statistics and Mathematics, WU Wien.
   \url{http://epub.wu-wien.ac.at/dyn/virlib/wp/showentry?ID=epub-wu-01_756&back=/}

   Helmut Strasser and Christian Weber (1999). On the asymptotic theory of permutation
   statistics. \emph{Mathematical Methods of Statistics}, \bold{8}, 220--250.

}
\examples{

### regression
data(airquality)
airq <- subset(airquality, !is.na(Ozone))
airct <- ctree(Ozone ~ ., data = airq, maxsurrogate = 3)
airct
plot(airct)
mean((airq$Ozone - predict(airct))^2)

### classification
irisct <- ctree(Species ~ .,data = iris)
irisct
plot(irisct)
table(predict(irisct), iris$Species)

### ordinal regression
data(mammoexp)
### attach scores
attr(mammoexp$ME, "scores") <- 1:3
attr(mammoexp$SYMPT, "scores") <- 1:4
attr(mammoexp$DECT, "scores") <- 1:3 
mammoct <- ctree(ME ~ ., data = mammoexp)
plot(mammoct)

### survival analysis
if (require(ipred)) {
    data(GBSG2, package = "ipred")
    GBSG2ct <- ctree(Surv(time, cens) ~ .,data = GBSG2)
    plot(GBSG2ct)
}
}
\keyword{tree}
