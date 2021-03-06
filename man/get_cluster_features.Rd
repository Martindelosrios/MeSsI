% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FeaturesEstimationFunctions.R
\name{get_cluster_features}
\alias{get_cluster_features}
\title{get_cluster_features}
\usage{
get_cluster_features(
  dat,
  ntotal = 0,
  name.groups = "clustersOutput.dat",
  ngal.lim = 30
)
}
\arguments{
\item{dat}{Data frame with a the angular coordinates (ra and dec in radians), the redshift, the color (g-r), the apparent magnitude (r) and the id of the group for each galaxy of each galaxy cluster. The columns must be named 'ra', 'dec' and 'z', 'color', 'mag' and 'id'.}
}
\value{
Data frame with the features of all the galaxy clusters.
}
\description{
This function estimates the features for a catalog of galaxy clusters.
}
\examples{
get_cluster_features(dat)
}
