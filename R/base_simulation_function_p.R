#' The guts for simulating & analysing datasets that are biased by over-dispersion in the detection probability (Scenario 3)
#'
#' It's recommended to use "simulation_function_p" instead of this one as 
#' "simulation_function_p" has better error handling (hopefully).
#' This function is one of the master functions in this package. It simulates and
#' analyzes a single dataset that is over-dispersed in detection probability, det_prob. 
#' To do this many times, use an apply function. It's currently simplified and 
#' assumes constant abundance, detection, and transect length. It does not (yet) 
#' simulate goodness-of-fit metrics. It simulates both point count and distance 
#' data and analyses both datasets using both unmarked and optim.
#'
#' @param n_sites number of sites (transects)
#' @param n_samps number of samples (replicates) per site
#' @param lambda  mean abundance at every site (single draw per site that stays constant across samples)
#' @param mean_det_prob average probability of being observed. Each simulation will draw a random detection probability from a beta distribution with a mean of mean_det_prob
#' @param sigma_beta_dist_p sigma on a beta distribution for the realized detection parameter. It's recommended to keep this below ~ 0.15. Above that the beta distribution becomes bimodal with modes near 0 and 1.
#' @param alpha   (optional) parameter of the beta distribution. It's recommended NOT to specify this parameter. Just specify the mean and sigma, and this value will be calculated to match the given mean and sigma. 
#' @param beta    (optional) parameter of the beta distribution. It's recommended NOT to specify this parameter. Just specify the mean and sigma, and this value will be calculated to match the given mean and sigma. 
#' @param W transect half-width (meters)
#' @param reps_to_analyze the number of samples/replicates to analyze. If NA, it will analyse all replicates in the data.
#' @param return What to return from the function call. Currently the only option is 'results'. May change this to only analyze simulated goodness-of-fit metrics.
#' @param savefilename The simulated datasets and results ARE saved to file (currently not optional). This provides the path and filename for saving the intermediate steps in the analysis. 
#'
#' @return if everything works well, it returns a data.frame with the results of simulating a single dataset, analyzing it in 4 ways, and calculating randomized quantile residuals a la Knape et al. 2018. It also saves a list with the simulated dataset, dataframe of results (minus rqr residual info), and the actual rqr residuals to a savefilename inside the folder path 'working directory'/results/Scenario 3/savefilename. If there is an error, the function returns NA and also saves a file to 'working directory'/Scenario 3/set x/errors with the simulated dataset and the results data.frame but no rq-residuals (it's basically assumed that the rq-residuals were the source of the error.) Similarly, with a warning the function returns the results data.frame and also saves a file to 'working directory'/Scenario 3/set x/warnings with the simulated dataset and the results data.frame but no rq-residuals (it's basically assumed that the rq-residuals were the source of the error.) The user will have to go back and try to calculate rq-residuals from the output later.
#'
#' @examples
#' base_simulation_function_p()
#'
#' @import tidyverse
#' @import nortest
#' 
#' @export
base_simulation_function_p = function(n_sites = 50, # number of sites
                                 n_samps = 6,  # number of samples per site
                                 lambda  = 10, # mean abundance at every site 
                                 mean_det_prob = 0.42, # mean detection probability
                                 sigma_beta_dist_p = 0.01, # sigma on a beta distribution for the realized detection parameter
                                 alpha = NA,
                                 beta = NA,
                                 W = 20,
                                 reps_to_analyze = 3, 
                                 #sampling_method = c('distance'), # , 'pointcount'
                                 #analysis_method = c('optim'), # , 'unmarked'
                                 #simulate_gof_pvals = T, 
                                 #simulate_gof_sims = 5,
                                 #simulate_gof_parallel = F, 
                                 return = 'results',
                                 savefilename = file.path('set 1', 'datasets', 'data')) {
  
  out <- tryCatch(
    {
      simdat <- sim_data_p(n_sites = n_sites,
                           n_samps = n_samps,
                           lambda  = lambda,
                           mean_det_prob = mean_det_prob, 
                           sigma_beta_dist_p = sigma_beta_dist_p, 
                           alpha = alpha, 
                           beta = beta, 
                           W = W)
      
      resL <- analyse_data(simulated_data = simdat, 
                          reps_to_analyze = reps_to_analyze, 
                          sampling_method = c('distance', 'pointcount'),
                          analysis_method = c('optim', 'unmarked'),
                          simulate_gof_pvals = FALSE,
                          return = return,
                          W = W)
      
      pco <- resL$inputs$include_pc_optim
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
        if(pco){
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
                      'sigma_beta_p' = sigma_beta_dist_p,
                      'sigma_mean'   = mean(simdat$inputs$sigmas, na.rm = T), # like sigma_true, but not quite the same since sigma varies
                      'alpha'        = simdat$inputs$alpha,
                      'beta'         = simdat$inputs$beta)
        
        # save the simulated data & residuals for later
        # (saving simulated data and analysed objects so that I can calculate c-hat values later)
        savd <- list(simulated_data = simdat,
                     analyzed_data = resL,
                     rqrM_optim_distance  = rqrM_od, 
                     rqrS_optim_distance  = rqrS_od, 
                     rqrO_optim_distance  = rqrO_od,
                     rqrM_optim_pointcount = ifelse(test = pco, yes = rqrM_op, no = NA),
                     rqrS_optim_pointcount = ifelse(test = pco, yes = rqrS_op, no = NA),
                     rqrO_optim_pointcount = ifelse(test = pco, yes = rqrO_op, no = NA),
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
                                'Scenario 3',
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
        if(pco){
          resid_pvals[2,] <- c(normtest(rqrM_op[,1], 's2_M_pval'),
                             normtest(rqrS_op, 's2_S_pval'),
                             if(length(rqrO_op) == 1) rep(rqrO_op, length(resid_pvals[2,])) else normtest(rqrO_op, 's2_O_pval'))
        }
        
        # shapiro-wilks, anderson-darling, and Lilliefors (Kolmogorov-Smirnov) tests of normality
        resid_pvals[3,] <- c(normtest(rqrM_ud[,1], 's3_M_pval'),
                             normtest(rqrS_ud, 's3_S_pval'),
                             if(length(rqrO_ud) == 1) rep(rqrO_ud, length(resid_pvals[3,])) else normtest(rqrO_ud, 's3_O_pval'))
        
        # shapiro-wilks, anderson-darling, and Lilliefors (Kolmogorov-Smirnov) tests of normality
        resid_pvals[4,] <- c(normtest(rqrM_up[,1], 's4_M_pval'),
                             normtest(rqrS_up, 's4_S_pval'),
                             if(length(rqrO_up) == 1) rep(rqrO_up, length(resid_pvals[4,])) else normtest(rqrO_up, 's4_O_pval'))
        
        # remove the point-count optim p-values if that analysis method isn't in use.
        if(! pco) resid_pvals <- resid_pvals[-2,]
      }
      
      # stuff to return
      cbind(res2,
            resid_pvals)
    },
    error=function(cond) {
      folders <- c(getwd(),
                   'results',
                   'Scenario 3',
                   dirname(dirname(savefilename)),
                   'errors')
      
      foldername <- do.call('file.path', as.list(folders))
      dir.create(foldername, recursive = T, showWarnings = FALSE)
      
      saveRDS(object =  list(message = cond, 
                             simdat = simdat,
                             results = resL,
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
                   'Scenario 3',
                   dirname(dirname(savefilename)),
                   'warnings')
      
      foldername <- do.call('file.path', as.list(folders))
      dir.create(foldername, recursive = T, showWarnings = FALSE)
      
      saveRDS(object =  list(message = cond, 
                             simdat = simdat,
                             results = resL,
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
      return(res2)
    }
  )    
  
  return(out)
}