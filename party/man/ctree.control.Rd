\name{ctree.control}
\alias{ctree.control}
\title{ Control for Conditional Tree Models }
\description{

  Various parameters that control aspects of the `ctree' fit.

}
\usage{
ctree.control(teststattype = c("quadform", "maxabs"), 
              testtype = c("Bonferroni", "MonteCarlo", "Raw"), 
              mincriterion = 0.95, minsplit = 20, stump = FALSE, 
              nresample = 9999, maxsurrogate = 0, mtry = 0, 
              savesplitstats = TRUE)
}
\arguments{
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
  \item{mtry}{ number of input variables randomly sampled as candidates 
               at each node for random forest like algorithms. The default
               \code{mtry = 0} means that no random selection takes place.}
  \item{savesplitstats}{ a logical determining if the process of standardized
                         two-sample statistics for split point estimate
                         is saved for each primary split.}
}
\details{

The stop
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


}
\value{
  An object of class \code{\link{TreeControl}}.
}
\keyword{misc}
