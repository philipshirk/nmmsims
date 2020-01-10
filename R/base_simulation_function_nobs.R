#' The guts for simulating & analysing datasets that are biased by over-dispersion in the number of observations (via double-counting) (Scenario 1)
#'
#' It's recommended to use "simulation_function_nobs" instead of this one as 
#' "simulation_function_nobs" has better error handling (hopefully).
#' This function is one of the master functions in this package. It simulates and
#' analyzes a single dataset that is over-dispersed in observations, n_obs. 
#' To do this many times, use an apply function. It's currently simplified and 
#' assumes constant abundance, detection, and transect length. It does not (yet) 
#' simulate goodness-of-fit metrics. It simulates both point count and distance 
#' data and analyses both datasets using both unmarked and optim.
#'
#' @param n_sites number of sites (transects)
#' @param n_samps number of samples (replicates) per site
#' @param lambda  mean abundance at every site (single draw per site that stays constant across samples)
#' @param alpha   probability of NOT being observed (i.e. alpha = 1 - detection probability)
#' @param gamma   probability of being observed exactly twice (gamma = 1 - alpha - beta)
#' @param beta    (optional) probability of being observed exactly once. If not specified, it's calculated as: beta = 1 - alpha - gamma
#' @param det_prob (optional) probability of being observed. If not specified, it's calculated as: det_prob = beta + gamma
#' @param sigma (optional) detection parameter (meters). Best to leave this blank and it will be calculated from det_prob.
#' @param W transect half-width (meters)
#' @param reps_to_analyze the number of samples/replicates to analyze. If NA, it will analyse all replicates in the data.
#' @param return What to return from the function call. Currently the only option is 'results'. May change this to only analyze simulated goodness-of-fit metrics.
#' @param savefilename The simulated datasets and results ARE saved to file (currently not optional). This provides the path and filename for saving the intermediate steps in the analysis. 
#'
#' @return if everything works well, it returns a data.frame with the results of simulating a single dataset, analyzing it in 4 ways, and calculating randomized quantile residuals a la Knape et al. 2018. It also saves a list with the simulated dataset, dataframe of results (minus rqr residual info), and the actual rqr residuals to a savefilename inside the folder path 'working directory'/results/Scenario 1/savefilename. If there is an error, the function returns NA and also saves a file to 'working directory'/Scenario 1/set x/errors with the simulated dataset and the results data.frame but no rq-residuals (it's basically assumed that the rq-residuals were the source of the error.) Similarly, with a warning the function returns the results data.frame and also saves a file to 'working directory'/Scenario 1/set x/warnings with the simulated dataset and the results data.frame but no rq-residuals (it's basically assumed that the rq-residuals were the source of the error.) The user will have to go back and try to calculate rq-residuals from the output later. 
#'
#' @examples
#' base_simulation_function_nobs()
#'
#' @import tidyverse
#' @import nortest
#' 
#' @export
base_simulation_function_nobs = function(n_sites = 50, # number of sites
                                    n_samps = 6,  # number of samples per site
                                    lambda  = 10, # mean abundance at every site 
                                    alpha   = (1 - 0.42), # probability of NOT being observed
                                    gamma   = 0.02,  # probably of being observed TWICE
                                    beta    = NA,    # probability of being observed ONCE (beta = 1 - alpha - gamma)
                                    det_prob = NA,   # mean detection probability
                                    sigma   = NA,    # detection parameter, calculated from det_prob
                                    W = 20,
                                    reps_to_analyze = 3, 
                                    # sampling_method = c('distance'), # , 'pointcount'
                                    # analysis_method = c('optim'), # , 'unmarked'
                                    # simulate_gof_pvals = T, 
                                    # simulate_gof_sims = 5,
                                    # simulate_gof_parallel = T, 
                                    return = 'results',
                                    savefilename = file.path('set 1', 'datasets', 'data')) {
  
  out <- tryCatch(
    {
      if(is.na(beta)) beta <- 1 - alpha - gamma
      if(is.na(det_prob)) det_prob = gamma + beta
      if( round(alpha + beta + gamma, 3) != 1.000) stop('The parameters alpha, beta, and gamma must add to 1.')
      
      simdat <- sim_data_nobs(n_sites = n_sites,
                              n_samps = n_samps,
                              lambda  = lambda,
                              alpha   = alpha,
                              gamma   = gamma,
                              sigma   = sigma,
                              W = W)
      
      resL <- analyse_data(simulated_data = simdat, 
                          reps_to_analyze = reps_to_analyze, 
                          sampling_method = c('distance', 'pointcount'),
                          analysis_method = c('optim', 'unmarked'),
                          simulate_gof_pvals = FALSE,
                          return = return,
                          W = W)
      
      # get the results data.frame
      res <- resL[['df']]
      
      # calculate several types of residuals & save everything
      {
        # Calculate randomized-quantile residuals
        # first, rename some parameters for convenience
        {
          K <- round(max(simdat$true_N)) * 10
          lam_od <- res$lambda_est[which(res$LL_method == 'optim' & 
                                           res$sampling_method == 'Distance')]
          lam_op <- res$lambda_est[which(res$LL_method == 'optim' & 
                                           res$sampling_method == 'pointcount')]
          lam_ud <- res$lambda_est[which(res$LL_method == 'unmarked' & 
                                           res$sampling_method == 'Distance')]
          lam_up <- res$lambda_est[which(res$LL_method == 'unmarked' & 
                                           res$sampling_method == 'pointcount')]
          
          sig_od <- res$sigma_est[which(res$LL_method == 'optim' & 
                                          res$sampling_method == 'Distance')]
          sig_op <- res$sigma_est[which(res$LL_method == 'optim' & 
                                          res$sampling_method == 'pointcount')]
          sig_ud <- res$sigma_est[which(res$LL_method == 'unmarked' & 
                                          res$sampling_method == 'Distance')]
          sig_up <- res$sigma_est[which(res$LL_method == 'unmarked' & 
                                          res$sampling_method == 'pointcount')]
          
          detp_od <- res$det_p_est[which(res$LL_method == 'optim' & 
                                           res$sampling_method == 'Distance')]
          detp_op <- res$det_p_est[which(res$LL_method == 'optim' & 
                                           res$sampling_method == 'pointcount')]
          detp_ud <- res$det_p_est[which(res$LL_method == 'unmarked' & 
                                           res$sampling_method == 'Distance')]
          detp_up <- res$det_p_est[which(res$LL_method == 'unmarked' & 
                                           res$sampling_method == 'pointcount')]
        }
        
        # residuals for optim & distance sampling
        {
          rqrM_od <- rqResMar_optim(n_obs = simdat$n_obs, 
                                    lambda_est = lam_od)
          
          rqrS_od <- rqResSum_optim(n_obs = simdat$n_obs, 
                                    K = K,
                                    p_est = detp_od,
                                    lambda_est = lam_od)
          
          # observation residuals are throwing some errors occassionally. I think only when sample sizes at a site are 0 and hence inadequate for estimating N distributions.
          # But anyway, I'm going to add in a trycatch to hopefully get rid of the error. 
          rqrO_od <- tryCatch(expr = {
            rqResObs_optim(n_obs = simdat$n_obs, 
                           p_est = detp_od,
                           K = K, 
                           lambda_est = lam_od,
                           sigma_est = sig_od, 
                           W = W,
                           y_list = simdat[['y_list']], 
                           sampling_method = 'distance')
          }, 
          error = function(cond) {
            print(cond)
            return(NA)
          })
        }
        
        # residuals for optim & point count sampling
        {
          rqrM_op <- rqResMar_optim(n_obs = simdat$n_obs,
                                    lambda_est = lam_op)
          
          rqrS_op <- rqResSum_optim(n_obs = simdat$n_obs, 
                                    K = K,
                                    p_est = detp_op,
                                    lambda_est = lam_op)
          
          # observation residuals are throwing some errors occassionally. I think only when sample sizes at a site are 0 and hence inadequate for estimating N distributions.
          # But anyway, I'm going to add in a trycatch to hopefully get rid of the error. 
          rqrO_op <- tryCatch(expr = {
            rqResObs_optim(n_obs = simdat$n_obs, 
                           p_est = detp_op,
                           K = K, 
                           lambda_est = lam_op, 
                           sigma_est = sig_op, 
                           W = W,
                           y_list = simdat[['y_list']], 
                           sampling_method = 'point count')  
          }, 
          error = function(cond) {
            print(cond)
            return(NA)
          })
        }
        
        # residuals for unmarked & distance sampling
        {
          rqrM_ud <- rqResMar_dist(umFit = resL$unmarked_dist)
          
          rqrS_ud <- rqResSum_dist(umFit = resL$unmarked_dist)
          
          rqrO_ud <- tryCatch(expr = {
            rqResObs_dist(umFit = resL$unmarked_dist)
          }, 
          error = function(cond) {
            print(cond)
            return(NA)
          })
        }
        
        # residuals for unmarked & point count sampling
        {
          rqrM_up <- rqResMar(umFit = resL$unmarked_pc)
          
          rqrS_up <- rqResSum(umFit = resL$unmarked_pc)
          
          rqrO_up <- tryCatch(expr = {
            rqResObs(umFit = resL$unmarked_pc)
          }, 
          error = function(cond) {
            print(cond)
            return(NA)
          })
        }
        
        res2 <- cbind(res, 
                      'n_sites'      = rep(n_sites, nrow(res)),
                      'n_samps'      = rep(n_samps, nrow(res)),
                      'lambda_true'  = rep(lambda,  nrow(res)),
                      'beta'         = rep(beta, nrow(res)), 
                      'gamma'        = rep(gamma, nrow(res)), 
                      'det_p_true'   = rep(det_prob, nrow(res)),
                      'sigma_true'   = rep(find_sigma(W = W, p_true = det_prob), nrow(res)))
        
        # save the simulated data & residuals for later
        # (saving simulated data and analysed objects so that I can calculate c-hat values later)
        savd <- list(simulated_data = simdat,
                     analyzed_data = resL,
                     rqrM_optim_distance  = rqrM_od, 
                     rqrS_optim_distance  = rqrS_od, 
                     rqrO_optim_distance  = rqrO_od,
                     rqrM_optim_pointcount = rqrM_op, 
                     rqrS_optim_pointcount = rqrS_op, 
                     rqrO_optim_pointcount = rqrO_op,
                     rqrM_unmarked_distance  = rqrM_ud, 
                     rqrS_unmarked_distance  = rqrS_ud, 
                     rqrO_unmarked_distance  = rqrO_ud,
                     rqrM_unmarked_pointcount = rqrM_up, 
                     rqrS_unmarked_pointcount = rqrS_up, 
                     rqrO_unmarked_pointcount = rqrO_up)
        savd[['analyzed_data']][['res']] <- res2
        
        # make sure the directly exists
        foldername <- file.path(getwd(),
                                'results',
                                'Scenario 1',
                                dirname(savefilename))
        
        filename <- paste0(basename(savefilename), 
                           '_',
                           gsub(pattern = ' ', replacement = '_', 
                                x = gsub(pattern = ':', 
                                         replacement = '',
                                         x = Sys.time())),
                           '_',
                           paste(Sys.info()[['nodename']], Sys.getpid(), sep='_'),
                           '.RDS')
        dir.create(path = foldername, recursive = TRUE, showWarnings = F)
        
        saveRDS(object = savd, 
                file = do.call(what = file.path, 
                               args = list(foldername,
                                        filename)))
      }
      
      # calculate pesky p-values
      {
        resid_pvals <- data.frame('rqres_Mar_pval_shapiro_wilks' = rep(NA,4),
                                  'rqres_Mar_pval_anderson_darling' = rep(NA,4),
                                  'rqres_Mar_pval_kolmogorov_Smirnov' = rep(NA,4),
                                  'rqres_SiteSum_pval_shapiro_wilks' = rep(NA,4),
                                  'rqres_SiteSum_pval_anderson_darling' = rep(NA,4),
                                  'rqres_SiteSum_pval_kolmogorov_Smirnov' = rep(NA,4),
                                  'rqres_Obs_pval_shapiro_wilks' = rep(NA,4),
                                  'rqres_Obs_pval_anderson_darling' = rep(NA,4),
                                  'rqres_Obs_pval_kolmogorov_Smirnov' = rep(NA,4))
        
        # shapiro-wilks, anderson-darling, and Lilliefors (Kolmogorov-Smirnov) tests of normality
        resid_pvals[1,] <- c(normtest(rqrM_od[,1], 's1_M_pval'), 
                             normtest(rqrS_od, 's1_S_pval'),
                             if(length(rqrO_od)==1) rep(rqrO_od, length(resid_pvals[1,])) else normtest(rqrO_od, 's1_O_pval'))
        
        # shapiro-wilks, anderson-darling, and Lilliefors (Kolmogorov-Smirnov) tests of normality
        resid_pvals[2,] <- c(normtest(rqrM_op[,1], 's2_M_pval'),
                             normtest(rqrS_op, 's2_S_pval'),
                             if(length(rqrO_op) == 1) rep(rqrO_op, length(resid_pvals[2,])) else normtest(rqrO_op, 's2_O_pval'))
        
        # shapiro-wilks, anderson-darling, and Lilliefors (Kolmogorov-Smirnov) tests of normality
        resid_pvals[3,] <- c(normtest(rqrM_ud[,1], 's3_M_pval'),
                             normtest(rqrS_ud, 's3_S_pval'),
                             if(length(rqrO_ud) == 1) rep(rqrO_ud, length(resid_pvals[3,])) else normtest(rqrO_ud, 's3_O_pval'))
        
        # shapiro-wilks, anderson-darling, and Lilliefors (Kolmogorov-Smirnov) tests of normality
        resid_pvals[4,] <- c(normtest(rqrM_up[,1], 's4_M_pval'),
                             normtest(rqrS_up, 's4_S_pval'),
                             if(length(rqrO_up) == 1) rep(rqrO_up, length(resid_pvals[4,])) else normtest(rqrO_up, 's4_O_pval'))
      }
      
      # stuff to return
      cbind(res2,
            resid_pvals)
    },
    error=function(cond) {
      folders <- c(getwd(),
                   'results',
                   'Scenario 1',
                   dirname(dirname(savefilename)),
                   'errors')
      
      foldername <- do.call('file.path', as.list(folders))
      dir.create(foldername, recursive = T, showWarnings = FALSE)
      
      saveRDS(object =  list(message = cond, 
                             simdat = if(exists(x = 'simdat')) simdat else NA,
                             results = if(exists(x = 'resL')) resL else NA,
                             n_obs = if(exists(x = 'simdat')) simdat$n_obs else NA,
                             lambda_est_op = if(exists(x = 'lam_op')) lam_op else NA), 
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
                   'warnings')
      
      foldername <- do.call('file.path', as.list(folders))
      dir.create(foldername, recursive = T, showWarnings = FALSE)
      
      saveRDS(object =  list(message = cond, 
                             simdat = if(exists(x = 'simdat')) simdat else NA,
                             results = if(exists(x = 'resL')) resL else NA,
                             n_obs = if(exists(x = 'simdat')) simdat$n_obs else NA,
                             lambda_est_op = if(exists(x = 'lam_op')) lam_op else NA), 
              file = file.path(foldername,
                               paste0(gsub(pattern = ' ', replacement = '_', 
                                           x = gsub(pattern = ':', 
                                                    replacement = '',
                                                    x = Sys.time())), 
                                      '.RDS')))
      # Choose a return value in case of error
      return(res2)
    }
  )    
  
  return(out)
}