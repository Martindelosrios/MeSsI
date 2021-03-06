% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FeaturesEstimationFunctions.R
\name{WeightedMclust}
\alias{WeightedMclust}
\title{WeightedMclust}
\usage{
WeightedMclust(group)
}
\arguments{
\item{group}{Data frame with the galaxies of the cluster.}
}
\value{
data frame with the same info that the input data frame plus the substructure id.
}
\description{
This function performs a wegihted mixture of gaussians.
}
\examples{
WeightedMclust(cat)
}
