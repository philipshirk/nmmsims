#' Approximate Posterior Distribution of N
#'
#' Function for getting an estimate of the posterior distribution of N from my 
#' optim models. This is from the unmarked::ranef() function: https://github.com/rbchan/unmarked/blob/master/R/ranef.R
#' and is needed for calculating Observation rq residuals.
#'
#' @param lambda_est estimated lambda value
#' @param n_obs number observed
#' @param K maximum N possible (typically obtained from n_obs)
#' @param p_est detection probability. Given if sampling_method == 'pointcount'. Calculated from if_distance__sigma_est if sampling_method == 'distance'
#' @param sampling_method either 'pointcount' or 'distance'
#' @param if_distance__sigma_est Estimated detection parameter, sigma. Only needed if sampling_method == distance
#' @param if_distance__W Transect half-width. Only needed if sampling_method == distance
#' @param if_distance__y_list array of distance data with row = site, col = repl, 3rd dim = observations. Only needed if sampling_method == distance
#'
#' @return PDF of possible counts. Formatted as an array of probabilities with row = site, col = repl, 3rd dim = observations.
#'
#' @export
get_N_distribution <- function(lambda_est, 
                               n_obs, 
                               K=NA, 
                               p_est = NA, 
                               sampling_method = 'pointcount',
                               if_distance__sigma_est = NA,
                               if_distance__W = NA,
                               if_distance__y_list = NA){
   # shorten some names
   sigma_est <- if_distance__sigma_est
   W <- if_distance__W
   y_list <- if_distance__y_list
   
   # number of lambdas (sites & replicates)
   R <- nrow(n_obs)
   
   # predicted lambda from the model
   lam <- rep(lambda_est, R)
   
   # get all the detection probability estimates (sites & replicates) 
   # from the fitted model
   if(is.na(p_est) & !is.na(sigma_est)) p_est <- calc_det_prob(W = W, sigma = sigma_est)
   p <- matrix(data = p_est, 
               nrow = nrow(n_obs), 
               ncol = ncol(n_obs))
   
   # get the max K value used to fit the model
   if(is.na(K)) K <- max(n_obs)*10
   
   # all possible N values
   N <- 0:K
   
   # get the observations (all sites & replicates)
   y <- n_obs
   
   # create an array to hold results
   # rows = sites
   # columns = N = 0:K
   post <- array(NA_real_, c(R, length(N), 1))
   # name the columns
   colnames(post) <- N
   
   # get the mixture model (i.e. poisson, negative binomial, etc)
   mix <- 'P'
   
   # cycle over every location
   for(i in 1:R) {
      switch(mix,
             # for a poisson mixture, use poisson PDF to calculate the 
             # probability of observing N given the lambda value for that site
             P  = f <- dpois(N, lam[i], log=TRUE),
             
             NB = {
                alpha <- exp(coef(object, type="alpha"))
                f <- dnbinom(N, mu=lam[i], size=alpha, log=TRUE)
             },
             ZIP = {
                psi <- plogis(coef(object, type="psi"))
                f <- (1-psi)*dpois(N, lam[i])
                f[1] <- psi + (1-psi)*exp(-lam[i])
                f <- log(f)
             })
      
      # empty vector of same length of N to hold the total 
      g <- rep(0, K+1)
      
      # cycle over replicates
      for(j in 1:ncol(y)) {
         # skip replicates that have NA observations OR NA detection probability
         if(is.na(y[i,j]) | is.na(p[i,j]))
            next
         
         # probability of observing y individuals given N and detection 
         # probability p
         g <- g + dbinom(y[i,j], N, p[i,j], log=TRUE)
         
         # if distance sampling, add in horizontal distances
         if(any(!is.na(match(x = tolower(sampling_method), 
                             table = c('dist', 'distance sampling', 'd', 'distance')))) &
            !any(is.na(y_list))){
            # horizontal distances at that site & replicate
            if(is.null(y_list[[i]][[j]])) next
            g <- g + sum(dnorm(x = y_list[[i]][[j]], 
                               mean = 0, 
                               sd = sigma_est, 
                               log = TRUE))
         }
      }
      
      # total probability of observations at a site IFF true abund = N
      fudge <- exp(f+g)
      
      # fudge is a vector over all N (at 1 site)
      # divide by the total probability to get absolute probabilities instead
      # of relative probabilities (at that site)
      post[i,,1] <- if(sum(fudge)==0) fudge else (fudge / sum(fudge))
   }
   return(post)
}