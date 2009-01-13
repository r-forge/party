\name{readingSkills}
\alias{readingSkills}
\docType{data}
\title{ Reading Skills }
\description{
  A toy data set to illustrate the spurious correlation
  between reading skills and shoe size
}
\usage{data("readingSkills")}
\format{
  A data frame with 200 observations on the following 4 variables.
  \describe{
    \item{\code{nativeSpeaker}}{a factor with levels \code{no} and \code{yes},
                                where \code{yes} indicates that the child
                                is a native speaker of the reading test.}
    \item{\code{age}}{age of the child in years.}
    \item{\code{shoeSize}}{shoe size of the child in cm.}
    \item{\code{score}}{raw score on the reading test.}
  }
}
\examples{

   set.seed(290875)
   readingSkills.cf <- cforest(score ~ ., data = readingSkills,
       control = cforest_unbiased(mtry = 2, ntree = 50))

   varimp(readingSkills.cf, conditional = TRUE)


}
\keyword{datasets}
