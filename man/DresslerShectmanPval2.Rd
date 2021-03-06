% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FeaturesEstimationFunctions.R
\name{DresslerShectmanPval2}
\alias{DresslerShectmanPval2}
\title{DresslerShectmanPval2}
\usage{
DresslerShectmanPval2(cat)
}
\arguments{
\item{cat}{data frame with a the angular coordinates (ra and dec in radians) and the redshift. The columns must be named 'ra', 'dec' and 'z'.}
}
\value{
The p-value of the Dressler-Shectman test corresponding to the galaxy cluster.
}
\description{
This function returns the p-value of the Dressler-Shectman statistic for a galaxy cluster with ngal = floor(sqrt(ngal)).
}
\examples{
DresslerShectmanPval2(cat)
}
