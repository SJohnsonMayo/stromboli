% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataobj_functions.R
\name{construct_distance}
\alias{construct_distance}
\title{Title}
\usage{
construct_distance(
  data.obj,
  unifrac.file = NULL,
  Phylum = "All",
  dist.RData = NULL,
  save.RData = NULL,
  filter.no = 0,
  rff = FALSE,
  dep = NULL,
  seed = 1234
)
}
\arguments{
\item{data.obj}{A data object created by load_data()}

\item{unifrac.file}{A file containing unifrac distances (Optional, default: NULL)}

\item{Phylum}{Phyla to include in the distance calculations (default: All)}

\item{dist.RData}{Previously generated distance object which can be loaded into the environment (Optional, default: NULL)}

\item{save.RData}{Filename to save distance object to (Optional, default: NULL)}

\item{filter.no}{Filter out OTUs with reads below this number (Default: 0)}

\item{rff}{Perform rarefaction, true or false(Default: FALSE)}

\item{dep}{Rarefaction depth, if performing rarefaction (Default: NULL)}

\item{seed}{Seed}
}
\value{
Distance object containing UniFrac, WUniFrac, GUniFrac, and Bray-Curtis distances
}
\description{
Title
}
