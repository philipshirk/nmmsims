#' Bootstrap Goodness-of-fit Statistics
#'
#' This function bootstraps p-values for the test statistics calculated with fitstats_optim().
#' It uses the sim_data() function in this package to simulate unbiased datasets.
#'
#' @param n_sites Number of sites to simulate
#' @param n_samps Number of samples (replicates) per site to simulate
#' @param lambda mean number of organisms at each site (1 draw per site)
#' @param sigma SD of half-normal detection function
#' @param det_prob mean detection probability (across distances from 0 to W)
#' @param sampling_method either "pointcount' or 'distance'
#' @param W Transect half-width
#' @param nsim number of simulations
#' # progress_bar = FALSE,
#' @param parallel Logical. Whether or not to use parallel processing for the simulations. Defaults to FALSE.
#'
#' @return data.frame with 3 columns (1 for each fit statistic) and nsim rows.
#'
#' @import parallel
#' 
#' @export
boot_fitstats <- function(n_sites, 
                          n_samps, 
                          lambda, 
                          sigma, 
                          det_prob, 
                          analysis_method = 'optim',
                          sampling_method = 'distance',
                          W = 20, 
                          nsim = 10,
                          # progress_bar = FALSE, 
                          parallel = FALSE){
   # note that this function assumes that there are NO NA's in the original dataset!!!
   #     - n_sites
   #     - n_samps
   #     - lambda (estimated)
   #     - det_sigma (estimated)
   
   if(parallel) require(parallel)
   
   if(!parallel){
      # data.frame for the fit.stats
      fs <- data.frame(SSE   = rep(NA, nsim),
                       Chisq = rep(NA, nsim),
                       freemanTukey = rep(NA, nsim))
      
      # create progress bar
      # if(progress_bar) pb <- txtProgressBar(min = 0, max = nsim, style = 3)
      
      # loop over every simulation
      for(i in 1:nsim) {
         # simulate a new dataset
         nd <- sim_data(n_sites = n_sites, 
                        n_samps = n_samps, 
                        lambda = lambda, 
                        det_prob = det_prob, 
                        sigma = sigma, 
                        W = W)
         
         # if there were NA's in the original dataset, make sure they're NA's in the simulated dataset, too
         # is.na(newd) <- is.na(obs)
         
         # re-fit the model with the simulated dataset
         ad <- analyse_data(simulated_data = nd, 
                            reps_to_analyze = n_samps, 
                            sampling_method = sampling_method, 
                            analysis_method = analysis_method, 
                            W = W, 
                            simulate_gof_pvals = FALSE,
                            return = 'gof')
         
         
         # update....
         fs[i,] <- ad[[analysis_method]][[sampling_method]]
         
         # if(progress_bar){
         # update progress bar
         # setTxtProgressBar(pb, i)
         # }
      }
   }
   
   if(parallel){
      coresToUse <- detectCores() - 1
      cl <- makeCluster(coresToUse)
      on.exit(stopCluster(cl))
      
      # variables to give to the cluster so that all cores can access them
      varList <- c('n_sites', 'n_samps', 'lambda', 'sigma', 'det_prob', 
                   'sampling_method', 'W', 'sim_data', 'analyse_data', 'nll',
                   'fitstats_optim')
      clusterExport(cl, varList, envir = environment())
      clusterEvalQ(cl, library(tidyverse))
      clusterEvalQ(cl, library(unmarked))
      # clusterEvalQ(cl, list2env(dots))
      
      fs_parallel <- parLapply(cl = cl, 
                               X = 1:nsim, 
                               fun = function(i) {
                                  nd <- sim_data(n_sites = n_sites, 
                                                 n_samps = n_samps, 
                                                 lambda = lambda, 
                                                 det_prob = det_prob, 
                                                 sigma = sigma, 
                                                 W = W)
                                  
                                  # if there were NA's in the original dataset, make sure they're NA's in the simulated dataset, too
                                  # is.na(newd) <- is.na(obs)
                                  
                                  # re-fit the model with the simulated dataset
                                  ad <- analyse_data(simulated_data = nd, 
                                                     reps_to_analyze = n_samps, 
                                                     sampling_method = sampling_method, 
                                                     analysis_method = 'optim', 
                                                     W = W, 
                                                     simulate_gof_pvals = FALSE,
                                                     return = 'gof')
                                  
                                  fs <- ad$optim[[sampling_method]]
                               })
      
      
      fs <- matrix(unlist(fs_parallel), 
                   nrow = length(fs_parallel), byrow = TRUE)
      colnames(fs) <- names(fs_parallel[[1]])
      fs <- as.data.frame(fs)
   }
   
   # close the progressbar
   # if(progress_bar) close(pb)
   
   return(fs)
}