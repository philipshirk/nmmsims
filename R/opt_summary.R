#' Summary of Goodness-of-fit Statistics
#'
#' I based this function on https://github.com/rbchan/unmarked/blob/master/R/boot.R
#' It returns a summary of the 3 fit statistics calculated by comparing fit of 
#' observed data (fitstats_optim) to fit of datasets simulated with the model 
#' estimates (boot_fitstats).
#'
#' @param fit_data output from the fitstats_optim function
#' @param fit_sims output from the boot_fitstats function
#'
#' @return data.frame with fit_obs = fit statistic of observed data; "mean(fit_obs - fit_sim)" = mean bias in the fit statistics; "StdDev(fit_obs - fit_sim)" = SD of the bias in fit statistics; "Pr(fit_sim > fit_obs)" = p-value of the fit statistics
#'
#' @export
opt_summary <- function(fit_data, fit_sims){
   # fit_data = output from fitstats_optim
   # fit_sims = output from boot_fitstats
   
   # fit stats from simulated datasets
   t.star <- fit_sims
   
   # fit stats from observed dataset
   t0     <- fit_data
   # number of simulated datasets
   nsim   <- nrow(t.star)
   
   # empty matrix for results
   biasMat <- pMat <- matrix(NA, 
                             nsim, 
                             length(t0))
   
   # cycle through each simulation 
   for(i in 1:nsim) {
      # difference between observed fit stats and simulated fit stats
      biasMat[i,] <- as.numeric((t0 - t.star[i,]))
      # are the simulated fit statistics more extreme than the observed fit statistics
      # I don't know why they subtract 1 here...
      pMat[i,] <- as.logical(abs(t.star[i,] - 1) > abs(t0 - 1))
   }
   
   # mean difference between observed fit stats and simulated fit stats
   bias <- colMeans(biasMat)
   # sd of the difference between observed fit stats and simulated fit stats
   bias.se <- apply(biasMat, 2, sd)
   # frequency with which the simulated fit stats are worse (larger) than the observed fit stats
   # I don't know why they add 1 here...
   p.val <- colSums(pMat) / (1 + nsim)
   
   stats <- data.frame("fit_obs" = t0, 
                       "mean(fit_obs - fit_sim)" = bias,
                       "StdDev(fit_obs - fit_sim)" = bias.se, 
                       "Pr(fit_sim > fit_obs)" = p.val,
                       check.names = FALSE)
   return(stats)
}