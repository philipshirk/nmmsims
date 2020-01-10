#' Conglomerate function to analyze simulated point count and/or distance sampling data
#'
#' This takes a single dataset of point count and/or distance sampling data and
#' analyses it in 1 to 4 ways.
#'
#' @param simulated_data the output from one of the sim_data functions in this package
#' @param reps_to_analyze the number of samples/replicates to analyze. If NA, it will analyse all replicates in the data.
#' @param sampling_method Analyze as distance sampling or point count. Options = 'distance' and/or 'pointcount'
#' @param analysis_method Analyze using optim or package unmarked. Options = 'optim' and/or 'unmarked'
#' @param simulate_gof_pvals Logical. Whether or not to simulate goodness-of-fit p-values, which is a VERY time-consuming process. Defaults to FALSE.
#' @param simulate_gof_sims Number of simulations to use for simulating goodness-of-fit p-values. Only used if simulate_gof_pvals = TRUE
#' @param simulate_gof_parallel Logical. Whether or not to use parallel processing to simulate goodness-of-fit p-values. Only used if simulate_gof_pvals = TRUE. Defaults to FALSE
#' @param W transect half-width (in meters)
#' @param return What to return from the function call. Currently the only option is 'results'. May change this to only analyze simulated goodness-of-fit metrics. 
#'
#' @return what will be returned
#'
#' @examples
#' sd <- sim_data()
#' analyze_data(sd)
#'
#' @import unmarked
#' @import dplyr
#' @importFrom tidyr pivot_longer
#' 
#' @export # this adds the following function to the NAMESPACE file
analyse_data <- function(simulated_data,
                         reps_to_analyze = c(3,6), # analyze the first 3 pointcount, or all 6?
                         sampling_method = c('distance', 'pointcount'), # use distance sampling, or repeat counts?
                         analysis_method = c('optim', 'unmarked'), # analyze using optim, or unmarked?
                         include_pc_optim = FALSE,
                         simulate_gof_pvals = FALSE,
                         simulate_gof_sims = 5,
                         simulate_gof_parallel = FALSE,
                         W = 20, # transect half-width
                         return = 'results'){
   require(tidyverse)
   
   inputs = list(reps_to_analyze = reps_to_analyze,
                 sampling_method = sampling_method,
                 analysis_method =analysis_method,
                 include_pc_optim = include_pc_optim,
                 simulate_gof_pvals = simulate_gof_pvals,
                 simulate_gof_sims = simulate_gof_sims,
                 simulate_gof_parallel = simulate_gof_parallel,
                 W = W, 
                 return = return)
   
   Nmat <- simulated_data[['true_N']][,1:reps_to_analyze]
   nmat <- simulated_data[['n_obs']][,1:reps_to_analyze]
   
   # maximum number of N to "integrate" over
   # trying to parboot pcount models from unmarked seems to require me to assign this in the global environment, which is weird...
   K <<- round(max(Nmat)) * 10
   
   # transect half-width
   if(is.na(W)) W <- max(simulated_data[['y']]) # if not set, use max(y), but I'll set it so that the det. prob. is precise.
   
   # number of sites and samples per site
   n_sites <- nrow(nmat)
   n_samps <- ncol(nmat)
   
   # format y (distance) data
   {
      # y's by site
      y_list_sparse <- sapply(X = simulated_data[['y_list']], FUN = function(l) unlist(l[1:reps_to_analyze]))
      # just the y's without regard to where they came from
      y <- unlist(y_list_sparse)
      
      # convert y distance data into long format
      {
         y_wide <- data.frame(site = factor(rep(1:n_sites, 
                                                each = n_samps)),
                              repl = factor(rep(1:n_samps, 
                                                n_sites)),
                              transectName = factor(paste(as.character(formatC(x = rep(1:n_sites, 
                                                                                       each = n_samps), 
                                                                               digits = 1, 
                                                                               flag='0')),
                                                          as.character(formatC(x = rep(1:n_samps, 
                                                                                       n_sites), 
                                                                               digits = 1, 
                                                                               flag='0')), 
                                                          sep='-'))) %>% cbind(as.data.frame(matrix(data = NA, nrow = (n_sites * n_samps), ncol = max(nmat))))
         
         for(m in 1:n_sites){
            for(j in 1:n_samps){
               sub <- simulated_data[['y_list']][[m]][[j]]
               if(length(sub)>0) y_wide[which(y_wide$site==m & y_wide$repl==j),(3 + (1:nmat[m,j]))] <- sub
            }
         }
         
         # convert the (wide) matrix into a (long) data.frame
         # y_mat is by site, and I want it to be by replicate...
         y_df <- tidyr::pivot_longer(data = y_wide, 
                                     values_to = 'distance', 
                                     cols = starts_with("V")) %>%
            na.omit() %>% 
            select(transectName, distance) %>% 
            as.data.frame()
      }
      
      # define distance categories for multinomial distribution
      {
         breaks <- seq(from = 0, to = W, by = 4)
      }
      
      # multinomial distances for GOF and/or fitting model
      {
         # format data for unmarked
         (unmk.dat <- formatDistData(distData = y_df,
                                     distCol="distance",
                                     transectNameCol="transectName",
                                     dist.breaks = breaks))
         y_obs <- as.data.frame(unmk.dat)
      }
   }
   
   # object to hold results
   {
      if('results' %in% return){
         # empty data.frame for output
         out <- data.frame(LL_method = NA,
                           sampling_method = NA,
                           lambda_est = NA,
                           lambda_SE  = NA,
                           lambda_lcl = NA,
                           lambda_ucl = NA,
                           sigma_est  = NA,
                           sigma_SE   = NA,
                           sigma_lcl  = NA,
                           sigma_ucl  = NA,
                           SSE        = NA,
                           SSE_pval   = NA,
                           Chisq      = NA,
                           Chisq_pval = NA,
                           Tukey      = NA,
                           Tukey_pval = NA,
                           # c_hat      = NA,
                           # rqres_marginal = NA,
                           # rqres_sitesum  = NA,
                           # rqres_obs      = NA,
                           # c_hat_marginal = NA,
                           # c_hat_sitesum  = NA,
                           det_p_est  = NA,
                           lambda_CV  = NA,
                           sigma_CV   = NA)[0,]
         hn_null <- NA
         pc_null <- NA
      } else {
         out <- list('optim' = list('distance' = NA,
                                    'pointcount' = NA),
                     'unmarked' = list('distance' = NA,
                                       'pointcount' = NA))
      }
   }
   
   if('optim' %in% analysis_method){
      
      if(any(!is.na(match(x = tolower(sampling_method), 
                          table = c('dist', 'distance', 'distancesamp', 'distance sampling'))))){
         # might be able to speed things up slightly by 
         # splitting off a GOF version of this analysis function
         # and paring it down to optimize for GOF
         # sig_start = ifelse(test = )
         
         # if distance_sampling, sigma = sd of normal, else detection probability
         opt_out <- optim(par  = c(sigma = log(5), lambda = log(max(nmat))),
                          fn    = nll,
                          K     = K,
                          n_dat = nmat,
                          W     = W,
                          y_dat = y_df$distance, # same as y
                          distance_sampling = TRUE)
         
         # exp(opt_out$par)
         
         # Model goodness of fit
         {
            # calculate chi-squared, Sum-of-Squares, and Freeman-Tukey GOF
            gof_obs  <- fitstats_optim(optim_out = opt_out, 
                                       n_obs = nmat, 
                                       y_obs = y_obs,
                                       W = W, 
                                       sampling_method = 'distance')
         }
         
         # if only looking for gof, stop here
         if(! ('results' %in% return)){
            out$optim$distance <- gof_obs
         } else {
            # if simulating gof...
            if(simulate_gof_pvals){
               fit_sims <- boot_fitstats(n_sites = n_sites, 
                                         n_samps = n_samps, 
                                         lambda = exp(opt_out$par['lambda']), 
                                         sigma = exp(opt_out$par['sigma']), 
                                         det_prob = (pnorm(W,0,exp(opt_out$par['sigma']))-0.5)/(dnorm(0,0,exp(opt_out$par['sigma']))*W), 
                                         sampling_method = 'distance', 
                                         W = W, 
                                         nsim = simulate_gof_sims, 
                                         # progress_bar = FALSE,
                                         parallel = simulate_gof_parallel)
               
               gof_pvals_d <- opt_summary(fit_data = gof_obs, fit_sims = fit_sims)
               
               # c-hat from Kery & Royle (2016) (just ratio of chi-squared values)
               # mean chi-squared value from simulations
               sim.chisq.mean <- mean(fit_sims$Chisq, na.rm=T)
               c_hat <- t0$Chisq / sim.chisq.mean
            }
            
            # format estimates in a dataframe
            {
               out.optimd <- data.frame(LL_method = 'optim',
                                        sampling_method = 'Distance',
                                        lambda_est = exp(opt_out$par['lambda']),
                                        lambda_SE  = NA,
                                        lambda_lcl = NA,
                                        lambda_ucl = NA,
                                        sigma_est  = exp(opt_out$par['sigma']),
                                        sigma_SE   = NA,
                                        sigma_lcl  = NA,
                                        sigma_ucl  = NA,
                                        SSE        = gof_obs['SSE'],
                                        SSE_pval   = if(simulate_gof_pvals) gof_pvals_d['SSE', 'Pr(fit_sim > fit_obs)'] else NA,
                                        Chisq      = gof_obs['Chisq'],
                                        Chisq_pval  = if(simulate_gof_pvals) gof_pvals_d['Chisq', 'Pr(fit_sim > fit_obs)'] else NA,
                                        Tukey      = gof_obs['freemanTukey'],
                                        Tukey_pval = if(simulate_gof_pvals) gof_pvals_d['freemanTukey', 'Pr(fit_sim > fit_obs)'] else NA,
                                        # c_hat = if(simulate_gof_pvals) c_hat else NA,
                                        # rqres_marginal = NA,
                                        # rqres_sitesum  = NA,
                                        # rqres_obs      = NA,
                                        # c_hat_marginal = NA,
                                        # c_hat_sitesum  = NA,
                                        det_p_est = NA) %>% 
                  mutate(det_p_est  = calc_det_prob(W = W, 
                                                    sigma = exp(opt_out$par['sigma'])),
                         lambda_CV = NA,
                         sigma_CV  = NA)
               
               out <- rbind(out, out.optimd)
            }
         }
      }
      
      if(include_pc_optim & any(!is.na(match(x = tolower(sampling_method),
                          table = c('pc', 'pointcount', 'point count', 'point', 'counts'))))){

         opt_out <- optim(par  = c(sigma = log(0.5 / (1 - 0.5)), lambda = log(max(nmat))),
                          fn    = nll,
                          K     = K,
                          n_dat = nmat,
                          W     = W,
                          y_dat = y,
                          distance_sampling = FALSE)

         # exp(opt_out$par)

         # Model fit
         {
            gof_obs  <- fitstats_optim(optim_out = opt_out,
                                       n_obs = nmat,
                                       W = W,
                                       sampling_method = 'pointcount')
         }

         # if only looking for gof, stop here
         if(! ('results' %in% return)){
            out$optim$pointcount <- gof_obs
         } else {
            # if simulating gof...
            if(simulate_gof_pvals){
               fit_sims <- boot_fitstats(n_sites = n_sites,
                                         n_samps = n_samps,
                                         lambda  = exp(opt_out$par['lambda']),
                                         sigma   = 1 / (1 + exp(-opt_out$par['sigma'])),
                                         det_prob = 1 / (1 + exp(-opt_out$par['sigma'])),
                                         sampling_method = 'pointcount',
                                         W = W,
                                         nsim = simulate_gof_sims,
                                         # progress_bar = FALSE,
                                         parallel = simulate_gof_parallel)

               gof_pvals_r <- opt_summary(fit_data = gof_obs, fit_sims = fit_sims)
            }

            # format estimates in a dataframe
            {
               out.optimr <- data.frame(LL_method = 'optim',
                                        sampling_method = 'pointcount',
                                        lambda_est = exp(opt_out$par['lambda']),
                                        lambda_SE  = NA,
                                        lambda_lcl = NA,
                                        lambda_ucl = NA,
                                        sigma_est  = 1 / (1 + exp(-opt_out$par['sigma'])),
                                        sigma_SE   = NA,
                                        sigma_lcl  = NA,
                                        sigma_ucl  = NA,
                                        SSE        = gof_obs['SSE'],
                                        SSE_pval   = if(simulate_gof_pvals) gof_pvals_r['SSE', 'Pr(fit_sim > fit_obs)'] else NA,
                                        Chisq      = gof_obs['Chisq'],
                                        Chisq_pval  = if(simulate_gof_pvals) gof_pvals_r['Chisq', 'Pr(fit_sim > fit_obs)'] else NA,
                                        Tukey      = gof_obs['freemanTukey'],
                                        Tukey_pval = if(simulate_gof_pvals) gof_pvals_r['freemanTukey', 'Pr(fit_sim > fit_obs)'] else NA,
                                        # c_hat = if(simulate_gof_pvals) c_hat else NA,
                                        # rqres_marginal = NA,
                                        # rqres_sitesum  = NA,
                                        # rqres_obs      = NA,
                                        # c_hat_marginal = NA,
                                        # c_hat_sitesum  = NA,
                                        det_p_est = NA) %>%
                  mutate(det_p_est  = sigma_est,
                         lambda_CV = NA,
                         sigma_CV  = NA)

               out <- rbind(out, out.optimr)
            }
         }
      }
   }
   
   if('unmarked' %in% analysis_method & ('results' %in% return)){
      require(unmarked)
      
      if(any(!is.na(match(x = tolower(sampling_method), 
                          table = c('dist', 'distance', 'distancesamp', 'distance sampling'))))){
         # format data for distance sampling
         {
            umf <- unmarkedFrameDS(y=unmk.dat, 
                                   survey="line", 
                                   dist.breaks=breaks, 
                                   tlength=rep(100, n_sites*n_samps),
                                   unitsIn="m")
         }
         
         # Analyze data in unmarked
         {
            hn_null <- distsamp(formula = ~1~1, 
                                data = umf,
                                keyfun="halfnorm", 
                                output="abund")
         }
         
         # Model fit
         {
            gof_und  <- fitstats(fm = hn_null)
         }
         
         # if only looking for gof, stop here
         if(! ('results' %in% return)){
            out$unmarked$distance <- gof_und
         } else {
            # simulate gof
            if(simulate_gof_pvals){
               gof_pvals_ud <- um_summary(parboot(hn_null, 
                                                  fitstats, 
                                                  nsim=simulate_gof_sims, 
                                                  parallel=simulate_gof_parallel))
            }
            
            # format estimates in a data.frame
            {
               und_lambda_CI <- confint(object = hn_null, type='state') %>% as.data.frame()
               und_sigma_CI  <- confint(object = hn_null, type='det') %>% as.data.frame()
               
               out.und <- data.frame(LL_method = 'unmarked',
                                     sampling_method = 'Distance',
                                     lambda_est = exp(hn_null@estimates@estimates$state@estimates),
                                     lambda_SE  = exp(sqrt(diag((hn_null@estimates@estimates$state@covMat)))),
                                     lambda_lcl = exp(und_lambda_CI[1,'0.025']),
                                     lambda_ucl = exp(und_lambda_CI[1,'0.975']),
                                     sigma_est  = exp(hn_null@estimates@estimates$det@estimates),
                                     sigma_SE   = exp(sqrt(diag((hn_null@estimates@estimates$det@covMat)))),
                                     sigma_lcl  = exp(und_sigma_CI[1,'0.025']),
                                     sigma_ucl  = exp(und_sigma_CI[1,'0.975']),
                                     SSE        = gof_und['SSE'],
                                     SSE_pval   = if(simulate_gof_pvals) gof_pvals_ud['SSE', 'Pr(fit_sim > fit_obs)'] else NA,
                                     Chisq      = gof_und['Chisq'],
                                     Chisq_pval  = if(simulate_gof_pvals) gof_pvals_ud['Chisq', 'Pr(fit_sim > fit_obs)'] else NA,
                                     Tukey      = gof_und['freemanTukey'],
                                     Tukey_pval = if(simulate_gof_pvals) gof_pvals_ud['freemanTukey', 'Pr(fit_sim > fit_obs)'] else NA,
                                     # c_hat = if(simulate_gof_pvals) c_hat else NA,
                                     # rqres_marginal = NA,
                                     # rqres_sitesum  = NA,
                                     # rqres_obs      = NA,
                                     # c_hat_marginal = NA,
                                     # c_hat_sitesum  = NA,
                                     det_p_est = NA) %>% 
                  mutate(det_p_est = calc_det_prob(W = W, sigma = sigma_est),
                         lambda_CV = round((lambda_SE / lambda_est) * 100, 1),
                         sigma_CV  = round((sigma_SE / sigma_est) * 100, 1))
               
               out <- rbind(out, out.und)
            }
         }
      }
      
      if(any(!is.na(match(x = tolower(sampling_method), 
                          table = c('pc', 'pointcount', 'point count', 'point', 'counts'))))){
         # format data for pcount
         {
            umf <- unmarkedFramePCount(nmat, siteCovs=NULL, obsCovs=NULL)
         }
         
         # fit the model
         {
            pc_null <- pcount(~1 ~1, umf, K=K)
         }
         
         # Model fit
         {
            gof_unr  <- fitstats(fm = pc_null)
         }
         
         # if only looking for gof, stop here
         if(! ('results' %in% return)){
            out$unmarked$pointcount <- gof_unr
         } else {
            # if simulating gof...
               if(simulate_gof_pvals){
                  pbr <- parboot(pc_null,
                                 fitstats, 
                                 nsim=simulate_gof_sims, 
                                 parallel=simulate_gof_parallel)
                  
                  gof_pvals_ur <- um_summary(pbr)
               }
            
            # format estimates in a data.frame
            {
               unr_lambda_CI <- confint(object = pc_null, type='state') %>% as.data.frame()
               unr_sigma_CI  <- confint(object = pc_null, type='det') %>% as.data.frame()
               
               # copied from https://github.com/cran/unmarked/blob/master/R/utils.R
               logistic <- function(x) {
                  1/(1 + exp(-x)) # logistic = inverse of logit
               }
               logistic.grad <- function(x) {
                  exp(-x)/(exp(-x)+1)^2 # gradient = ???
               }
               
               # detection parameter estimate
               # logistic(pc_null@estimates@estimates$det@estimates)
               # detection parameter SE estimate
               # gradient^2 * coef. variance
               # logistic.grad(pc_null@estimates@estimates$det@estimates)^2 * pc_null@estimates@estimates$det@covMat
               
               out.unr <- data.frame(LL_method = 'unmarked',
                                     sampling_method = 'pointcount',
                                     lambda_est = exp(pc_null@estimates@estimates$state@estimates),
                                     lambda_SE  = exp(sqrt(diag(pc_null@estimates@estimates$state@covMat))), # sqrt(backTransform(obj = pc_null, type = 'state')@covMat)
                                     lambda_lcl = exp(unr_lambda_CI[1,'0.025']),
                                     lambda_ucl = exp(unr_lambda_CI[1,'0.975']),
                                     sigma_est  = backTransform(obj = pc_null, type = 'det')@estimate, 
                                     # = 1 / (1 + exp(-sig)) = logistic(pc_null@estimates@estimates$det@estimates)
                                     sigma_SE   = sqrt(backTransform(obj = pc_null, type = 'det')@covMat),
                                     # converting the SE back is a little tricky. See here: 
                                     # https://www.andrewheiss.com/blog/2016/04/25/convert-logistic-regression-standard-errors-to-odds-ratios-with-r/
                                     # and here: 
                                     # https://github.com/cran/unmarked/blob/master/R/utils.R
                                     # sqrt(logistic.grad(pc_null@estimates@estimates$det@estimates)^2 * pc_null@estimates@estimates$det@covMat)
                                     
                                     sigma_lcl  = logistic(unr_sigma_CI[1,'0.025']),
                                     sigma_ucl  = logistic(unr_sigma_CI[1,'0.975']),
                                     SSE        = gof_unr['SSE'],
                                     SSE_pval   = if(simulate_gof_pvals) gof_pvals_ur['SSE', 'Pr(fit_sim > fit_obs)'] else NA,
                                     Chisq      = gof_unr['Chisq'],
                                     Chisq_pval  = if(simulate_gof_pvals) gof_pvals_ur['Chisq', 'Pr(fit_sim > fit_obs)'] else NA,
                                     Tukey      = gof_unr['freemanTukey'],
                                     Tukey_pval = if(simulate_gof_pvals) gof_pvals_ur['freemanTukey', 'Pr(fit_sim > fit_obs)'] else NA,
                                     # c_hat = if(simulate_gof_pvals) c_hat else NA,
                                     # rqres_marginal = NA,
                                     # rqres_sitesum  = NA,
                                     # rqres_obs      = NA,
                                     # c_hat_marginal = NA,
                                     # c_hat_sitesum  = NA,
                                     det_p_est = NA) %>% 
                  mutate(det_p_est = sigma_est,
                         lambda_CV = round((lambda_SE / lambda_est) * 100, 1),
                         sigma_CV  = round((sigma_SE / sigma_est) * 100, 1))
               
               out <- rbind(out, out.unr)
            }
         }
      }
   }
   
   return(
      # if only looking for the gof simulations, just return that
      if(! ('results' %in% return)){
         out
      } else {
         # otherwise return a lot more info
         list(df = out, 
              unmarked_dist = hn_null,
              unmarked_pc   = pc_null,
              inputs = inputs)
      }
   ) 
}