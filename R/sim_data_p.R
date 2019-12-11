#' Simulate biased (via over-dispersion in p) point count and distance sampling data
#'
#' Generate point counts and/or distance sampling data that is biased by over-dispersion
#' in the detection probability, p. Distance sampling data is created with each dataset,
#' but it can just be ignored if you only want the point count data. This is a 
#' rudementary function. Currently it assumes constant density, detection probability,
#' transects length, etc.
#'
#' @param n_sites number of sites (transects)
#' @param n_samps number of samples (replicates) per site
#' @param lambda mean abundance at every site
#' @param mean_det_prob mean detection probability across all sites and replicates (different value drawn for each site and replicate)
#' @param sigma_beta_dist_p SD of a beta distribution from which the realized detection probability is drawn (keep this below ~ 0.16)
#' @param alpha parameter of the beta distribution. This is deterministically calculated from mean_det_prob and sigma_beta_dist_p if they are provided.
#' @param beta parameter of the beta distribution. This is deterministically calculated from mean_det_prob and sigma_beta_dist_p if they are provided.
#' @param W transect half-width (meters)
#'
#' @return Returns a list of 4 items: 1. "true_N" is a matrix of true abundance values for each site (row) and sample (column). 2. "n_obs" is a matrix of the number of observed organisms at each site (row) and sample (column). 3. "y_list" is a list of vectors. Each vector holds the distance data for a single survey. 4. "inputs" is a list of input values.
#'
#' @examples
#' sim_data_p()
#' 
#' @importFrom truncnorm rtruncnorm
#'
#' @export
sim_data_p <- function(n_sites = 50, # number of sites
                       n_samps = 6,  # number of samples per site
                       lambda  = 10, # mean abundance at every site 
                       mean_det_prob = 0.42, # mean detection probability
                       sigma_beta_dist_p = 0.01, # sigma on a beta distribution for the realized detection parameter
                       alpha = NA,
                       beta = NA,
                       W = 20){ # transect half-width
   
   inputs <- list(n_sites = n_sites,
                  n_samps = n_samps,
                  lambda = lambda,
                  mean_det_prob = mean_det_prob,
                  sigma_beta_dist_p = sigma_beta_dist_p,
                  alpha = alpha,
                  beta = beta,
                  W = W)
   
   # mean of beta = a/(a+b)
   # var of beta = ab/((a+b)^2 (a+b+1)) = (pi * (1-pi)) / (a + b + 1)
   sigma_b = sigma_beta_dist_p
   
   # alpha & beta will probably be set for a given # of simulations, 
   # but if they're not passed into the function, calculate them
   if(is.na(alpha)) alpha = (mean_det_prob * (1-mean_det_prob)) / 
      (sigma_b^2 * (1 + (1 - mean_det_prob)/mean_det_prob)) - 
      1 / (1 + (1 - mean_det_prob)/mean_det_prob)
   if(is.na(beta)) beta = alpha / mean_det_prob - alpha
   
   # simulate p-values (1 per site-replicate combo)
   p_sim <- rbeta(n = (n_sites * n_samps),
                  shape1 = alpha, 
                  shape2 = beta)
   
   # calculate detection sigmas for each site-replicate
   sigma_sim <- sapply(X = p_sim, FUN = function(p) find_sigma(start = 0.1, W = W, p_true = p))
   
   # convert into a matrix just to help me keep the n_sites & n_samps straight
   sigmat <- matrix(data = sigma_sim, 
                    nrow = n_sites, 
                    ncol = n_samps, 
                    byrow = TRUE)
   
   # visualize beta distribution
   {
      # mean of beta dist
      # (pi_b <- alpha / (alpha + beta))
      # SD of beta dist
      # (sig_b <- sqrt(alpha * beta / ((alpha + beta)^2 * (alpha + beta + 1))))
      
      # visualize beta distribution(s)
      # x_plot <- seq(from = 0.001, to = .999, by = 0.001)
      # y_plot <- sapply(X = 1:6, FUN = function(x) dbeta(x = x_plot, shape1 = alpha[x], shape2 = beta[x]))
      # view <- setNames(object = as.data.frame(cbind(y_plot, x_plot)), nm = c('0.01', '0.02', '0.03', '0.04', '0.08', '0.16', 'x')) %>% pivot_longer(data = ., cols = c('0.01', '0.02', '0.03', '0.04', '0.08', '0.16'), names_to = 'sigma', values_to = 'y')
      
      # ggplot(view, aes(x = x, y = y, color=sigma))+ theme_bw() + geom_line() + xlab(label = 'Detection probability')+ ylab(label = 'Relative frequency')+ scale_x_continuous(labels = scales::percent), width = 4, height = 4, units = 'in')
   }
   
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
   
   # simulate the observation history with binomial dist
   n_sim <- rbinom(n = (n_sites * n_samps),
                   size = N_sim,
                   prob = p_sim)
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
      sigma = sigmat[m,]
      
      y_list[[m]] <- mapply(FUN = function(n, sigma) 
         truncnorm::rtruncnorm(n = n, 
                               a = 0, 
                               b = W, 
                               mean = 0, 
                               sd = sigma), 
         n = n, sigma = sigma)
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