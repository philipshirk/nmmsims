#' Simulate unbiased point count and distance sampling data
#'
#' Generate unbiased point counts and/or distance sampling data. Distance 
#' sampling data is created with each dataset, but it can just be ignored if you
#' only want the point count data. This is a rudementary function. Currently it
#' assumes constant density, detection probability, transects length, etc.
#'
#' @param n_sites number of sites (transects)
#' @param n_samps number of samples (replicates) per site
#' @param lambda  mean number of organisms per site
#' @param det_prob mean detection probability (must be between 0 and 1)
#' @param sigma detection parameter (meters). Just leave this blank and it will be calculated from det_prob
#' @param W transect half-width (meters)
#'
#' @return Returns a list of 4 items: 1. "true_N" is a matrix of true abundance values for each site (row) and sample (column). 2. "n_obs" is a matrix of the number of observed organisms at each site (row) and sample (column). 3. "y_list" is a list of vectors. Each vector holds the distance data for a single survey. 4. "inputs" is a list of input values.
#'
#' @examples
#' sim_data()
#'
#' @importFrom truncnorm rtruncnorm
#' 
#' @export
sim_data <- function(n_sites = 50, # number of sites
                     n_samps = 6,  # number of samples per site
                     lambda  = 10,
                     det_prob = 0.42, # mean detection probability
                     sigma = NA, # detection parameter, calculated from det_prob
                     W = 20){ # transect half-width
   
   # create a list of input values to include in the output
   inputs <- list(n_sites = n_sites,
                  n_samps = n_samps,
                  lambda = lambda,
                  det_prob = det_prob,
                  sigma = sigma,
                  W = W)
   
   # calculate sigma from the detection probability and W
   # but only if it's not already defined (to speed up simulations where it stays constant)
   if(is.na(sigma)){sigma <- find_sigma(W = W, p_true = det_prob)}
   
   # total number of individuals (NOT over-dispersed)
   N_sim <- rep(rpois(n = n_sites, 
                      lambda = lambda), 
                each = n_samps)
   # convert into a matrix just to help me keep the n_sites & n_samps straight
   Nmat <- matrix(data = N_sim, 
                  nrow = n_sites, 
                  ncol = n_samps, 
                  byrow = TRUE)
   
   # simulate the observation history with binomial dist
   n_sim <- rbinom(n = (n_sites * n_samps),
                   size = N_sim,
                   prob = det_prob)
   
   # convert into a matrix just to help me keep the n_sites & n_samps straight
   nmat <- matrix(data = n_sim, 
                  nrow = n_sites, 
                  ncol = n_samps, 
                  byrow = TRUE)
   
   # generate distance data from a truncated normal distribution
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
   
   # number expected   
   # sum_n_exp <- sum(N_sim) * det_prob
   # sum_n_obs <- sum(n_sim)
   
   return(list('true_N' = Nmat,
               'n_obs'  = nmat,
               'y_list' = y_list, 
               inputs = inputs))
}