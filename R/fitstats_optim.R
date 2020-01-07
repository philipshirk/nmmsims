#' Goodness-of-fit Statistics
#'
#' I based this function on "fitstats()" from https://github.com/rbchan/unmarked/blob/master/R/boot.R
#' It's a function returning three goodness-of-fit statistics - sum-of-squares, 
#' Chi-squared, and the Freeman-Tukey test-statistic - calculated using my optim
#' estimates rather than a model fit in "unmarked".
#'
#' @param out fit from optim
#' @param n_obs vector of observations at each site & replicate
#' @param y_obs matrix of the number of observations in each distance category (from formatted unmarked data). Only used if sampling_method = "distance"
#' @param W transect half-width (in meters)
#' @param sampling_method either "pointcount" or "distance"
#'
#' @return vector of 3 fit statistics
#'
#' @import unmarked
#' 
#' @export
fitstats_optim <- function(optim_out, n_obs, y_obs = NA, W, sampling_method){
   out <- optim_out
   
   # parameter that the model estimated
   lambda    <- exp(out$par['lambda'])
   
   # calculate numbers observed and expected based on the methods
   switch(EXPR = sampling_method,
          'pointcount' = {
             # numbers observed
             obs <- n_obs
             # logistic (inverse logit) of sigma
             det_prob <- 1/(1 + exp(- out$par['sigma']))
             # numbers expected
             exp <- lambda * det_prob
          }, 
          'distance' = {
             # numbers observed
             obs <- y_obs
             det_sigma <- exp(out$par['sigma'])
             det_prob  <- (pnorm(q = W,mean = 0,sd = det_sigma) - 0.5) / (dnorm(x = 0,mean = 0,sd = det_sigma) * W)
             
             # probabilities for each of the bins
             # first get the bins from the data
             brks <- c(0, as.numeric(as.vector(sapply(X = names(y_obs), 
                                                      FUN = function(x) gsub(pattern = "[[:punct:]]([0-9]+)[[:punct:]]([0-9]+)[[:punct:]]", 
                                                                             replacement = "\\2",
                                                                             x = x)))))
             bp <- pnorm(q = brks, mean = 0, sd = det_sigma) - 0.5
             probs <- sapply(X = 2:length(brks), FUN = function(x) (bp[x] - bp[x-1]) / 0.5)
             
             # expected counts
             exp <- as.data.frame(matrix(rep( (lambda * det_prob) * probs, 
                                              nrow(n_obs) * ncol(n_obs)), 
                                         nrow = nrow(n_obs) * ncol(n_obs), 
                                         ncol = length(probs), byrow = T))
          })
   
   # residuals = observed - expected
   resids <- obs - exp
   
   # calculate fit statistics
   sse      <- sum(resids^2)
   chisq    <- sum((obs - exp)^2 / exp)
   freeTuke <- sum((sqrt(obs) - sqrt(exp))^2)
   
   output <- c(SSE=sse, 
               Chisq=chisq, 
               freemanTukey=freeTuke)
   return(output)
}