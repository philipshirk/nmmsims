#' Goodness-of-fit Statistics
#'
#' I copied this function from https://github.com/rbchan/unmarked/blob/master/R/boot.R
#' It's a function returning a summary of the 3 fit statistics calculated using
#' fitstat().
#'
#' @param object the output from the fitstats() function
#'
#' @return data.frame with fit_obs = fit statistic of observed data; "mean(fit_obs - fit_sim)" = mean bias in the fit statistics; "StdDev(fit_obs - fit_sim)" = SD of the bias in fit statistics; "Pr(fit_sim > fit_obs)" = p-value of the fit statistics
#'
#' @import unmarked
#' 
#' @export
um_summary <- function(object){
   t.star <- object@t.star
   t0 <- object@t0
   nsim <- nrow(t.star)
   biasMat <- pMat <- matrix(NA, nsim, length(t0))
   for(i in 1:nsim) {
      biasMat[i,] <- t0 - t.star[i,]
      pMat[i,] <- abs(t.star[i,] - 1) > abs(t0 - 1)
   }
   bias <- colMeans(biasMat)
   bias.se <- apply(biasMat, 2, sd)
   p.val <- colSums(pMat) / (1 + nsim)
   stats <- data.frame("fit_obs" = t0, 
                       "mean(fit_obs - fit_sim)" = bias,
                       "StdDev(fit_obs - fit_sim)" = bias.se, 
                       "Pr(fit_sim > fit_obs)" = p.val,
                       check.names = FALSE)
   return(stats)
}