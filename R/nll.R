#' Negative Log Likelihood Function for Estimating Sigma and Lambda
#'
#' function description
#'
#' @param par parameters to be guessed (in this case c(log(sigma), log(lambda))). The parameters are logged so that they're always positive.
#' @param K The maximum possible true abundance, N. I've been using round(lambda * 10) as a default.
#' @param n_dat matrix of observations at all sites (rows) and samples/replicates (columns)
#' @param W Transect half-width (meters). If not set, it uses the maximum observed distance.
#' @param y_dat vector of distances (meters) for all observations
#' @param distance_sampling logical. Whether or not to incorporate distance data into likelihood.
#'
#' @return total negative log likelihood of the data given the parameters.
#'
#' @importFrom matrixStats rowProds
#'
#' @export # this adds the following function to the NAMESPACE file
nll <- function(par,
                K,
                n_dat,
                W = 20,
                y_dat,
                distance_sampling = TRUE) {

   # if W is not set, use the maximum observed distance
   if(is.na(W)) W <- max(y_dat)

   # work with the non-logged parameters
   lambda <- exp(par[2])

   if(distance_sampling){
      # if distance_sampling, sigma = sd of normal (use exp & norm. dist. to convert to [0,1])
      sigma <- exp(par[1])
      # calculate detection probability from sigma
      Pa  = (pnorm(q = W,mean = 0,sd = sigma) - 0.5) / (dnorm(x = 0,mean = 0,sd = sigma) * W)
   } else Pa = exp(par[1]) / (1 + exp(par[1]))
   # else sigma = detection probability (use logit to convert to [0,1])

   # number of sites
   # just using M b/c that's the name they use in unmarked
   M <- nrow(n_dat)
   # number of pointcount per site
   # just using J b/c that's the name they use in unmarked
   J <- ncol(n_dat)
   # possible values of N (K = arbitrarily large cutoff value)
   k <- 0:K
   # every possible k value for every site
   k.km  <- rep(k, M)
   # every k value for every site and replicate
   k.kmj <- rep(k.km, J)
   # observations for every possible k value
   n.kmj <- rep(n_dat, each = (K+1))
   # detection probability for every site and every k value
   p.km  <- rep(rep(Pa, M), each = (K + 1))
   # detection probability for every site, k value, and replicate
   p.kmj <- rep(p.km, J)
   # lambda for every site and every k value
   l.km  <- rep(rep(lambda, M), each = (K + 1))
   # lambda for every site, k value, and replicate
   l.kmj <- rep(l.km, J)

   # binomial likelihood
   bin.kmj <- dbinom(x = n.kmj,    # observations, repeated over every possible k, site, and replicate
                     size = k.kmj, # every possible k value, repeated over every site and replicate
                     prob = p.kmj) # every p-value, repeated for every possible k, site, and replicate
   # just in case there were NA's in the data, but I'm not simulating it that way now
   # bin.kmj[which(is.na(bin.kmj))] <- 1

   # put the data into a matrix where each column is a replicate and each row is the same site (repeated for k rows)
   bmat <- matrix(data = bin.kmj, nrow = M * (K+1), ncol = J)

   # multiply likelihoods for each replicate (end up with 1 likelihood per K per site)
   g.km <- matrixStats::rowProds(bmat)

   # likelihood of the k values given the lambda value
   f.km <- dpois(x = k.km, lambda = l.km)

   # matrix of likelihood of (K and obs)
   dens.m.mat <- matrix(data = f.km * g.km,
                        nrow = M, # each row is a unique site
                        ncol = K + 1, # each column is a possible k-value
                        byrow = TRUE)

   # site_to_plot <- 1
   # plot(x = k, y = dens.m.mat[site_to_plot,])
   # plot(x = k, y = f.km[((site_to_plot-1)*(K+1)+1):((site_to_plot)*(K+1))])
   # plot(x = k, y = g.km[((site_to_plot-1)*(K+1)+1):((site_to_plot)*(K+1))])

   # get the likelihood of lambda at each site by summing over all k-values
   dens.m <- rowSums(dens.m.mat)

   # total negative log likelihood = sum of likelihood of obs and distances
   nll <- -sum(log(dens.m))

   if(distance_sampling){
      dll <- dnorm(x = y_dat,
                   mean = 0,
                   sd = sigma,
                   log = TRUE)
      nll <- nll - sum(dll)
   }

   return(nll)
}
