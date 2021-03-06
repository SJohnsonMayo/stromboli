% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/beta_diversity.R
\name{perform_betadisper_test}
\alias{perform_betadisper_test}
\title{Perform Betadispersion test}
\usage{
perform_betadisper_test(
  data.obj,
  dist.obj,
  dist.names = c("UniFrac", "GUniFrac", "WUniFrac", "BC"),
  grp.name
)
}
\arguments{
\item{data.obj}{Data object created by load_data()}

\item{dist.obj}{Distance object created by construct_distance()}

\item{dist.names}{Distance metrics to include in PERMANOVA test. Default includes "UniFrac", "GUniFrac", "WUniFrac", and "BC"}

\item{grp.name}{Variable of interest}
}
\value{
A table containing the results of the BETADISPER test
}
\description{
Perform Betadispersion test
}
\examples{
data("Constipation")
betadisper <- perform_betadisper_test(data.obj.rff, dist.obj.rff, grp.name="Visit")
}
