% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alpha_diversity.R
\name{generate_alpha_diversity}
\alias{generate_alpha_diversity}
\title{Generate alpha diversity object}
\usage{
generate_alpha_diversity(
  data.obj,
  rarefy = TRUE,
  depth = NULL,
  iter.no = 5,
  measures = c("Observed", "Chao1", "Shannon", "InvSimpson"),
  seed = 123
)
}
\arguments{
\item{data.obj}{A data object created by load_data()}

\item{rarefy}{Perform rarefaction, TRUE or FALSE (Default: TRUE)}

\item{depth}{Rarefaction depth}

\item{iter.no}{Number of iterations (Default: 5)}

\item{measures}{List of alpha diversity measures to include (Default includes "Observed", "Chao1", "Shannon", and "InvSimpson")}

\item{seed}{Random seed to use}
}
\value{
Data frame containing the alpha diversities for all samples in data.obj
}
\description{
Generate alpha diversity object
}
\examples{
data("Constipation")
alpha.obj <- generate_alpha_diversity(data.obj, measures = c("Observed", "Chao1", "Shannon", "InvSimpson"), depth=10000, iter.no=5)
}
