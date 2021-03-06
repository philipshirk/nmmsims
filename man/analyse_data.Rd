% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyse_data.R
\name{analyse_data}
\alias{analyse_data}
\title{Conglomerate function to analyze simulated point count and/or distance sampling data}
\usage{
analyse_data(
  simulated_data,
  reps_to_analyze = c(3, 6),
  sampling_method = c("distance", "pointcount"),
  analysis_method = c("optim", "unmarked"),
  simulate_gof_pvals = FALSE,
  simulate_gof_sims = 5,
  simulate_gof_parallel = FALSE,
  W = 20,
  return = "results"
)
}
\arguments{
\item{simulated_data}{the output from one of the sim_data functions in this package}

\item{reps_to_analyze}{the number of samples/replicates to analyze. If NA, it will analyse all replicates in the data.}

\item{sampling_method}{Analyze as distance sampling or point count. Options = 'distance' and/or 'pointcount'}

\item{analysis_method}{Analyze using optim or package unmarked. Options = 'optim' and/or 'unmarked'}

\item{simulate_gof_pvals}{Logical. Whether or not to simulate goodness-of-fit p-values, which is a VERY time-consuming process. Defaults to FALSE.}

\item{simulate_gof_sims}{Number of simulations to use for simulating goodness-of-fit p-values. Only used if simulate_gof_pvals = TRUE}

\item{simulate_gof_parallel}{Logical. Whether or not to use parallel processing to simulate goodness-of-fit p-values. Only used if simulate_gof_pvals = TRUE. Defaults to FALSE}

\item{W}{transect half-width (in meters)}

\item{return}{What to return from the function call. Currently the only option is 'results'. May change this to only analyze simulated goodness-of-fit metrics.}
}
\value{
what will be returned
}
\description{
This takes a single dataset of point count and/or distance sampling data and
analyses it in 1 to 4 ways.
}
\examples{
sd <- sim_data()
analyze_data(sd)

}
