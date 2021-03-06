% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/beta_diversity.R
\name{generate_ordination2}
\alias{generate_ordination2}
\title{Generate ordination plots for beta diversity analysis}
\usage{
generate_ordination2(
  data.obj,
  dist.obj,
  dist.names = c("UniFrac", "GUniFrac", "WUniFrac", "BC"),
  grp.name,
  adj.name = NULL,
  emp.lev = NULL,
  strata = NULL,
  pca.method = "cmd"
)
}
\arguments{
\item{data.obj}{Data object created by load_data()}

\item{dist.obj}{Distance object created by construct_distance()}

\item{dist.names}{List of distance metrics to include in plot. Default includes "UniFrac", "GUniFrac", "WUniFrac", and "BC"}

\item{grp.name}{Variable of interest}

\item{adj.name}{List of variables to adjust for / covariates}

\item{emp.lev}{Levels of grp.name to exclude (Optional)}

\item{strata}{Variable indicating strata}

\item{pca.method}{PCA method, either "cmd" or "nmds" (Default: "cmd")}
}
\value{
A ggplot object containing ordination plots for all distance measures indicated by \code{dist.names}
}
\description{
Generate ordination plots for beta diversity analysis
}
\examples{
data("Constipation")
generate_ordination2(data.obj, dist.obj, grp.name="Visit")
}
