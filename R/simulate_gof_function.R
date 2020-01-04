#' Primary function for simulating & analysing datasets that are biased by over-dispersion in the number of observations (via double-counting) (Scenario 1)
#'
#' This function is one of the master functions in this package. It simulates and
#' analyzes a single dataset that is over-dispersed in observations, n_obs. 
#' To do this many times, use an apply function. It's currently simplified and 
#' assumes constant abundance, detection, and transect length. It does not (yet) 
#' simulate goodness-of-fit metrics. It simulates both point count and distance 
#' data and analyses both datasets using both unmarked and optim.
#'
#' @param filename full path and name of data file saved by the "simulate_function_nobs" function.
#' @param simulate_gof_sims The number of times to simulate goodness of fit values
#' @param simulate_gof_parallel Logical. Whether or not to do the goodness of fit simulations in parallel. The default is FALSE, since I typically run the entire function in parallel, not just parts of the internal calculations.
#' @param savefilename The simulated datasets and results ARE saved to file (currently not optional). This provides the path and filename for saving the intermediate steps in the analysis. 
#' 
#' @return if everything works well, it returns a list with the results of taking the data.frame output from simulation_function_nobs, simulating lots of datasets with set each parameter estimates, calculates a GOF metric for each one, summarizes them, and returns p-values for each one. It puts the results right back into the data.frame and returns and whole list with the simulated dataset, and original analysis. The only difference is the p-values added into the results data.frame.
#'
#' @import tidyverse
#' 
#' @export

# NOT UPDATED. I SHOULD CHANGE THIS TO READ IN SIMULATED DATA FILES AND CALCULATE THINGS BASED ON THOSE
# INSTEAD OF SIMULATING NEW DATA SETS
simulate_gof_function = function(filename = NA,
                                      simulate_gof_sims = 5,
                                      simulate_gof_parallel = F,
                                      savefilename = 'set 3/datasets/data'){
  out <- tryCatch(
    {
      # read in the data
      simdatlist <- readRDS(file = filename)
      
      # shorten some names
      # simulated data
      simdat <- simdatlist[['simulated_data']]
      
      # data.frame with simulation results (minus gof p-values)
      res.df <- simdatlist[['analyzed_data']][['res']]
      
      # which methods were used to analyze data
      analysis_methods <- as.character(res.df$LL_method)
      sampling_methods <- as.character(res.df$sampling_method)
      
      reps_to_analyze <- simdatlist$analyzed_data$inputs$reps_to_analyze
      Nmat <- simdat[['true_N']][,1:reps_to_analyze]
      nmat <- simdat[['n_obs']][,1:reps_to_analyze]
      
      # number of sites and samples per site
      n_sites <- nrow(nmat)
      n_samps <- ncol(nmat)
      
      # lambda estimates
      lambda_est <- simdatlist$analyzed_data$df$lambda_est
      # sigma estimates
      sigma_est  <- simdatlist$analyzed_data$df$sigma_est
      # estimated detection probability
      det_prob_est <- simdatlist$analyzed_data$df$det_p_est
      
      W <- simdat$inputs$W
      
      for (i in 1:nrow(res.df)) {
        gof_pvals <- switch(analysis_methods[i],
                            'optim' = {
                              fit_sims <- boot_fitstats(n_sites = n_sites, 
                                                        n_samps = n_samps, 
                                                        lambda = lambda_est[i], 
                                                        sigma = sigma_est[i], 
                                                        det_prob = det_prob_est[i], 
                                                        sampling_method = sampling_methods[i], 
                                                        analysis_method = analysis_methods[i],
                                                        W = W, 
                                                        nsim = simulate_gof_sims, 
                                                        parallel = simulate_gof_parallel)
                              
                              opt_summary(fit_data = gof_obs, fit_sims = fit_sims)
                            },
                            'unmarked' = {
                              unmarked_fit <- switch(sampling_methods[i],
                                                     'distance' = simdatlist$analyzed_data$unmarked_dist,
                                                     'pointcount' = simdatlist$analyzed_data$unmarked_pc)
                              
                              fit_sims <- parboot(unmarked_fit,
                                                  fitstats, 
                                                  nsim = simulate_gof_sims, 
                                                  parallel = simulate_gof_parallel)
                              
                              um_summary(fit_sims)
                            }
        )
       res.df[i,c('SSE_pval', 'Chisq_pval', 'Tukey_pval')] <- gof_pvals$`Pr(fit_sim > fit_obs)`
      }
      
      simdatlist$analyzed_data$res <- res.df
      simdatlist
    },
    error=function(cond) {
      folders <- c(getwd(),
                   'results',
                   'Scenario 1',
                   dirname(dirname(savefilename)),
                   'errors',
                   'gof')
      
      foldername <- do.call('file.path', as.list(folders))
      dir.create(foldername, recursive = T, showWarnings = FALSE)
      
      saveRDS(object =  list(message = cond, 
                             simdat = simdat,
                             results = res.df,
                             n_obs = simdat$n_obs,
                             lambda_est_op = lam_op), 
              file = file.path(foldername,
                               paste0(gsub(pattern = ' ', replacement = '_', 
                                           x = gsub(pattern = ':', 
                                                    replacement = '',
                                                    x = Sys.time())),
                                      paste(Sys.info()[['nodename']], Sys.getpid(), sep='_'), 
                                      '.RDS')))
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      folders <- c(getwd(),
                   'results',
                   'Scenario 1',
                   dirname(dirname(savefilename)),
                   'warnings',
                   'gof')
      
      foldername <- do.call('file.path', as.list(folders))
      dir.create(foldername, recursive = T, showWarnings = FALSE)
      
      saveRDS(object =  list(message = cond, 
                             simdat = simdat,
                             results = res.df,
                             n_obs = simdat$n_obs,
                             lambda_est_op = lam_op), 
              file = file.path(foldername,
                               paste0(gsub(pattern = ' ', replacement = '_', 
                                           x = gsub(pattern = ':', 
                                                    replacement = '',
                                                    x = Sys.time())), 
                                      '.RDS')))
      # Choose a return value in case of error
      return(res.df)
    }
  )
  return(out)
}