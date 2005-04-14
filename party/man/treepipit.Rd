\name{treepipit}
\alias{treepipit}
\docType{data}
\title{ Tree Pipit (Anthus trivialis) Forest Data }
\description{
  Data on the population density of tree pipits in Franconian oak forests
  including variables describing the forest ecosystem.
}
\usage{data(treepipit)}
\format{
  A data frame with 86 observations describing 86 stands with 
  the following 10 variables:
  \describe{
    \item{counts}{the number of tree pipits observed.}
    \item{age}{The ages of the overstorey oaks taken from
               forest data.}
    \item{coverstorey}{Percentage of cover of canopy overstorey.
                       The crown cover is described relative to a fully
                       stocked stand. Very dense overstorey with multiple 
                       crown cover could reach values greater than 100
                       percent.}
    \item{coverregen}{Percentage of cover of regeneration and shrubs.}
    \item{meanregen}{Mean height of regeneration and shrubs.}
    \item{coniferous}{Percentage of coniferous trees.}
    \item{deadtree}{Number of dead trees per grid.}
    \item{cbpiles}{Number of crowns and branch piles.
                   All laying crowns and branch piles per hectare were 
                   counted. These were induced by logging and the creation 
                   of wind breaks.}
    \item{ivytree}{Number of ivied trees. All ivied trees were counted per
                   hectare.}
    \item{fdist}{Distance to the forest edge.
                 The closest distance to the forest edge was measured 
                 from the centre of each grid.}
  }
}
\details{
  This study is based on fieldwork conducted in three lowland oak forests in
  the Franconian region of northern Bavaria close to Uffenheim, Germany.
  Diurnal breeding birds were samples from March to 
  June 2002 five times, using a
  quantitative grid mapping. Each grid was a 1-ha square. In total, 86 sample
  sites weres established. All individuals were counted in time intervals of
  7 min/grid during slow walks along the middle of the grid with a stop in the 
  centre. Environmental factors were measured for each grid.
}
\references{

  Joerg Mueller and Torsten Hothorn (2005). Maximally Selected Two-Sample 
  Statistics as a new Tool for the Identification and Assessment of Habitat 
  Factors with an Application to Breeding Bird Communities in Oak Forests. 
  \emph{Journal of Forest Research} (accepted).

}
\examples{
data(treepipit)
ctree(counts ~ ., data = treepipit)
}
\keyword{datasets}
