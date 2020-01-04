#' Simulate biased (via double-counting) point count and distance sampling data (Scenario 1)
#'
#' Generate point counts and/or distance sampling data that is biased by double-counting.
#' Distance sampling data is created with each dataset, but it can just be ignored if you
#' only want the point count data. This is a rudementary function. Currently it
#' assumes constant density, detection probability, transects length, etc.
#'
#' @param n_sites number of sites (transects)
#' @param n_samps number of samples (replicates) per site
#' @param lambda  mean number of organisms per site
#' @param alpha   probability of NOT being observed (i.e. alpha = 1 - detection probability; and detection probability = gamma + beta)
#' @param gamma   probability of an organism being observed exactly TWICE
#' @param sigma   detection parameter (meters). Just leave this blank and it will be calculated from det_prob
#' @param W transect half-width (meters)
#'
#' @return Returns a list of 4 items: 1. "true_N" is a matrix of true abundance values for each site (row) and sample (column). 2. "n_obs" is a matrix of the number of observed organisms at each site (row) and sample (column). 3. "y_list" is a list of vectors. Each vector holds the distance data for a single survey. 4. "inputs" is a list of input values.
#'
#' @examples
#' sim_data_nobs()
#' 
#' @importFrom truncnorm rtruncnorm
#'
#' @export
sim_data_nobs <- function(n_sites = 50, # number of sites
                          n_samps = 6,  # number of samples per site
                          lambda  = 10, # mean abundance at every site 
                          alpha   = 0.58, # probability of NOT being observed (i.e. alpha = 1 - detection probability)
                          gamma   = 0.02, # probably of being observed TWICE
                          sigma   = NA, # detection parameter, calculated from det_prob
                          W = 20){ # transect half-width
   
   # probability of being observed ONCE (beta = 1 - alpha - gamma)
   beta  <- (1 - alpha - gamma)
   
   # probability of being observed
   det_prob <- beta + gamma
   
   inputs <- list(n_sites = n_sites,
                  n_samps = n_samps,
                  lambda = lambda,
                  alpha = alpha,
                  gamma = gamma,
                  beta = beta,
                  det_prob = det_prob,
                  sigma = sigma,
                  W = W)
   
   # calculate sigma from the detection probability and W
   # but only if it's not already defined (to speed up simulations where it stays constant)
   if(is.na(sigma)){sigma <- find_sigma(W = W, p_true = det_prob)}
   
   # simulate the true number of individuals at each site
   N_sim0 <- rpois(n = n_sites, 
                   lambda = lambda)
   
   # repeat the true number for each replicate survey
   N_sim <- rep(x = N_sim0, 
                each = n_samps)
   
   # convert into a matrix just to help me keep the n_sites & n_samps straight
   Nmat <- matrix(data = N_sim, 
                  nrow = n_sites, 
                  ncol = n_samps, 
                  byrow = TRUE)
   
   # simulate the observation history with multinomial instead of binomial dist
   multi_sim <- sapply(X = N_sim, 
                       FUN = function(n) rmultinom(n = 1, 
                                                   size = n, 
                                                   prob = c(alpha, beta, gamma)))
   
   # add together number observed: 1* number obs 1x + 2*number observed twice to get total observations
   n_sim <- apply(X = multi_sim, 
                  MARGIN = 2, 
                  FUN = function(x){sum(x[2], 2*x[3])})
   
   # convert into a matrix just to help me keep the n_sites & n_samps straight
   nmat <- matrix(data = n_sim, 
                  nrow = n_sites, 
                  ncol = n_samps, 
                  byrow = TRUE)
   
   # generate distance data from a truncated normal distribution
   # but in a structure where I can easily subset it later on
   y_list <- list()
   
   for(m in 1:n_sites){
      n <- nmat[m,]
      y_list[[m]] <- lapply(X = n, 
                            FUN = function(x) truncnorm::rtruncnorm(n = x, 
                                                                    a = 0, 
                                                                    b = W, 
                                                                    mean = 0, 
                                                                    sd = sigma))
   }
   
   # y_list is now a list with n_sites primary elements and n_samps sub-lists for each primary element
   
   # naive number expected
   # sum_n_exp  <- sum(N_sim) * det_prob
   # biased number expected
   # sum_n_exp2 <- beta * sum(N_sim) + gamma * sum(N_sim) * 2
   # sum_n_obs <- sum(n_sim)
   
   return(list('true_N' = Nmat,
               'n_obs'  = nmat,
               'y_list' = y_list, 
               inputs = inputs))
}