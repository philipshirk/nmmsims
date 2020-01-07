#' Goodness-of-fit Statistics
#'
#' I copied this function from https://github.com/rbchan/unmarked/blob/master/R/boot.R
#' It's a function returning three goodness-of-fit statistics: sum-of-squares, Chi-squared, and the Freeman-Tukey test-statistic.
#'
#' @param fm a model fit in the unmarked package
#'
#' @return vector of 3 fit statistics
#'
#' @import unmarked
#' 
#' @export
fitstats <- function(fm) {
   observed <- getY(fm@data) # counts in each bin
   expected <- fitted(fm)    # expected in each bin
   resids <- residuals(fm)   # obs - exp
   sse <- sum(resids^2)      # 
   chisq <- sum((observed - expected)^2 / expected)
   freeTuke <- sum((sqrt(observed) - sqrt(expected))^2)
   out <- c(SSE=sse, 
            Chisq=chisq, 
            freemanTukey=freeTuke)
   return(out)
}