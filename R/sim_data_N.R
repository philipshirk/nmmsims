#' Simulate biased (via over-dispersion in N) point count and distance sampling data
#'
#' Generate point counts and/or distance sampling data that is biased by over-dispersion
#' in the true abundance, N. Distance sampling data is created with each dataset,
#' but it can just be ignored if you only want the point count data. This is a 
#' rudementary function. Currently it assumes constant density, detection probability,
#' transects length, etc.
#'
#' @param n_sites number of sites (transects)
#' @param n_samps number of samples (replicates) per site
#' @param lambda1  partial mean abundance at every site (single draw per site that stays constant across samples)
#' @param lambda2  partial mean abundance at every replicate (new draw for every sample at every site)
#' @param det_prob mean detection probability
#' @param sigma detection parameter (meters). Just leave this blank and it will be calculated from det_prob
#' @param W transect half-width (meters)
#'
#' @return Returns a list of 4 items: 1. "true_N" is a matrix of true abundance values for each site (row) and sample (column). 2. "n_obs" is a matrix of the number of observed organisms at each site (row) and sample (column). 3. "y_list" is a list of vectors. Each vector holds the distance data for a single survey. 4. "inputs" is a list of input values.
#'
#' @examples
#' sim_data_N()
#' 
#' @importFrom truncnorm rtruncnorm
#'
#' @export
sim_data_N <- function(n_sites = 50, # number of sites
                       n_samps = 6,  # number of samples per site
                       lambda1 = 5,  # mean abundance at every site (constant across pointcount)
                       lambda2 = 5,  # mean abundance at every replicate
                       det_prob = 0.42, # mean detection probability
                       sigma = NA, # detection parameter, calculated from det_prob
                       W = 20){ # transect half-width
   
   inputs <- list(n_sites = n_sites,
                  n_samps = n_samps,
                  lambda1 = lambda1,
                  lambda2 = lambda2,
                  det_prob = det_prob,
                  sigma = sigma,
                  W = W)
   
   # calculate sigma from the detection probability and W
   # but only if it's not already defined (to speed up simulations where it stays constant)
   if(is.na(sigma)){sigma <- find_sigma(W = W, p_true = det_prob)}
   
   # simulate the constant N across pointcount in each site
   v_sim0 <- rpois(n = n_sites, lambda = lambda1)
   v_sim  <- rep(v_sim0, each = n_samps)
   
   # simulate the variable N across pointcount in each site
   A_sim <- rpois(n = (n_sites * n_samps), lambda = lambda2)
   
   # total number of individuals (over-dispersed) at all sites for all pointcount
   N_sim <- v_sim + A_sim
   
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
   
   # number expected   
   # sum_n_exp <- sum(N_sim) * det_prob
   # sum_n_obs <- sum(n_sim)
   
   return(list('true_N' = Nmat,
               'n_obs'  = nmat,
               'y_list' = y_list,
               inputs = inputs))
}