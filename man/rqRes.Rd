% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rqRes.R
\name{rqRes}
\alias{rqRes}
\title{Randomized Quantile Residuals}
\usage{
rqRes(y, pFun, ..., show.warnings = FALSE)
}
\arguments{
\item{y}{a matrix of observed counts, with rows = sites and columns = replicates}

\item{pFun}{CDF function of the count data. For this package it's just always Poisson}

\item{...}{additional parameters fed to pFun}
}
\value{
a matrix of residuals with the same dimensions as the observations, y.
}
\description{
Generic function for Marginal, Site-sum, or Observation Randomized-quantile
residuals, as defined in Knape et al. 2018 and available for unmarked models
in the nmixgof package. This is a generic function used in the background of
other functions.
}
