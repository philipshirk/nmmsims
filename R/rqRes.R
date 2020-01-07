#' Randomized Quantile Residuals
#'
#' Generic function for Marginal, Site-sum, or Observation Randomized-quantile
#' residuals, as defined in Knape et al. 2018 and available for unmarked models
#' in the nmixgof package. This is a generic function used in the background of
#' other functions.
#'
#' @param y a matrix of observed counts, with rows = sites and columns = replicates
#' @param pFun CDF function of the count data. For this package it's just always Poisson
#' @param ... additional parameters fed to pFun
#'
#' @return a matrix of residuals with the same dimensions as the observations, y.
#'
#' @export
 rqRes = function(y, pFun, ..., show.warnings = FALSE) {
    # b is F(z)
    b = pFun(y, ...)
    # ensure that a is the same length as b
    a = b*0
    # a is F(z - 1)
    a = pFun(y-1, ...)

    # create an object of NA's the same type and dimensions as y
    if (is.null(dim(y))) {
       res = numeric(length(y)) + NA
    } else {
       res = array(data= NA, dim = dim(y))
    }

    # test which elements of y, a, and b are NOT NA
    isnna = which(!is.na(y) & !is.na(a) & !is.na(b))

    # for elements that are not NA, calculate the residuals
    res[isnna] = rqRes0(a[isnna], b[isnna])

    # throw a warning for any residuals that are infinite
    if(any(is.infinite(res)) & show.warnings){
       #browser()
       warning("Some residuals infinite.")
       saveRDS(list(n_obs = y, lambda_est = ...), file = './warnings/data_that_generated_the_most_recent_rqRes_warning.RDS')
    }

    # return the residuals
    res
 }

#' Helps calculate Randomized Quantile Residuals
#'
#' Used in the Generic function for Marginal, Site-sum, or Observation Randomized-quantile
#' residuals, as defined in Knape et al. 2018 and available for unmarked models
#' in the nmixgof package. This is a generic function used in the background of
#' other functions.
#'
#' @param a lower limit on uniform distribution = F(z-1) in Knape 2018
#' @param b upper limit on uniform distribution = F(z) in Knape 2018
#'
#' @return a single value
#'
#' @export
rqRes0 = function(a, b) {
    stopifnot(length(a) == length(b))
    qnorm(runif(length(a), a, b))
}

#' Marginal Randomized Quantile Residuals
#'
#' Calculates Marginal Randomized-quantile residuals for a model fit with optim.
#' As defined in Knape et al. 2018 and available for unmarked models in the nmixgof package.
#'
#' @param n_obs a matrix of observed counts, with rows = sites and columns = replicates
#' @param lambda_est estimated lambda (mean abundance across sites) from fitted model
#'
#' @return a matrix of residuals with the same dimensions as the observations, n_obs.
#'
#' @export
rqResMar_optim = function(n_obs, lambda_est){
   rqRes(y = n_obs, pFun = ppois, lambda = lambda_est)
}

#' Site-Sum Randomized Quantile Residuals
#'
#' Calculates Site-Sum Randomized-quantile residuals for a model fit with optim.
#' As defined in Knape et al. 2018 and available for unmarked models in the nmixgof package.
#'
#' @param n_obs a matrix of observed counts, with rows = sites and columns = replicates
#' @param K maximum possible true mean abundance, N
#' @param p_est estimated mean (across distances = 0 to W) detection probability
#' @param lambda_est estimated lambda (mean abundance across sites) from fitted model
#'
#' @return a matrix of residuals with the same dimensions as the observations, n_obs.
#'
#' @export
rqResSum_optim = function(n_obs, K, p_est, lambda_est, show.warnings = FALSE) {

   # get the predicted lambda value at each site
   lam = rep(lambda_est, nrow(n_obs)) # predict(umFit, type="state")[,1]

   # get the predicted detection probability values for every site-replicate combo
   # p = getP(umFit)
   p = matrix(data = p_est,
              nrow = nrow(n_obs),
              ncol = ncol(n_obs))

   # create an empty matrix for the counts at each site
   cumProb = matrix(0, nrow = nrow(n_obs), ncol = 2)

   # pick the PDF/density function of the appropriate distribution (used to analyze the data)
   dfun = function(N) dpois(N, lam)

   # loop over possible N values (i.e. k) to calculate the
   for (N in 0:K) {
      cumProb = cumProb +
         kronecker(dfun(N), t(c(1,1))) *
         pbinsum(n_obs,
                 rep(N, nrow(n_obs)),
                 p)[,2:3]
   }

   # calculate
   res = rqRes0(cumProb[,1], cumProb[,2])

   #throw a warning if any residuals are infinite
   if (any(is.infinite(res)) & show.warnings) {
      warning(paste(sum(is.infinite(res)), " residuals infinite."))
   }

   # return the residuals
   res
}

#' Observation Randomized Quantile Residuals
#'
#' Calculates Observation Randomized-quantile residuals for a model fit with optim.
#' As defined in Knape et al. 2018 and available for unmarked models in the nmixgof package.
#'
#' @param n_obs a matrix of observed counts, with rows = sites and columns = replicates
#' @param K maximum possible true mean abundance, N
#' @param lambda_est estimated lambda (mean abundance across sites) from fitted model
#' @param sigma_est estimated sigma from fitted model (only needed if distance sampling)
#' @param W transect half-width. (only needed if distance sampling)
#' @param p_est estimated mean (across distances = 0 to W) detection probability (only needed if point counts)
#' @param y_list list of distance values with elements = sites and replicates elements within each element
#' @param sampling_method either 'pointcount' or 'distance'
#'
#' @return a matrix of residuals with the same dimensions as the observations, n_obs.
#'
#' @export
rqResObs_optim = function(n_obs,
                          K,
                          lambda_est,
                          sigma_est = NA,
                          W = NA,
                          p_est = NA,
                          y_list = NA,
                          sampling_method,
                          show.warnings = FALSE) {

   # NA's for each site
   rN = integer(nrow(n_obs)) + NA

   # Draw/sample random N values (1 per site) from the estimated PDF of possible N values
   # sample the PDF of possible N values for each site.
   rN = apply(X = get_N_distribution(lambda_est = lambda_est,
                                     n_obs = n_obs,
                                     p_est = p_est,
                                     K = K,
                                     sampling_method = sampling_method,
                                     if_distance__W = W,
                                     if_distance__y_list = y_list,
                                     if_distance__sigma_est = sigma_est)[,,1],
              # unmarked::ranef = "Estimate posterior distributions of the random variables (latent abundance or occurrence) using empirical Bayes methods. So this returns an estimated probability that the true (latent) N == each possible k value. Basically, a PDF for N across all possible k's. They take the first slice, which corresponds to the first sampling period and just ignore all other sampling periods." These are passed to the probability weights argument in the base::sample() function.
              MARGIN = 1, # across all possible N's within a single site.
              FUN = sample,
              x = K + 1, # items from which to choose (can choose 0, so have to add 1). This is just the number or rows in X
              size = 1,        # choose only a single value
              replace = FALSE) - 1 # subtract 1 so that it includes 0 and doesn't go above K.

   # if get_N_distribution worked, continue
   if ( !all(is.na(rN)) ) {
      # getP is a function from unmarked
      # returns a matrix of estimated detection probabilities for each site (rows) and replicate (columns).
      # p = getP(umFit)
      # for my simulations, I can just use the same detection prob.
      if(is.na(p_est) & !is.na(sigma_est)) p_est <- calc_det_prob(W = W, sigma = sigma_est)
      p = matrix(data = p_est,
                 nrow = nrow(n_obs),
                 ncol = ncol(n_obs))
      
      # get the residuals
      res = rqRes(n_obs,
                  pFun = pbinom,
                  # Computes the generalised kronecker product of two arrays, X and Y
                  size = kronecker(X = rN,
                                   Y = t(rep(1, ncol(p)))),
                  prob=p)
      
      if (any(is.infinite(res)) & show.warnings) {
         #throw a warning if any residuals are infinite
         warning(paste(sum(is.infinite(res)), " residuals infinite."))
      }
   } else {
      # if get_N_distribution did NOT work, then just return NA
      res <- NA
   }

   # return the residuals
   res
}

#' Marginal Randomized Quantile Residuals
#'
#' Calculates Marginal Randomized-quantile residuals for a point count model fit with unmarked.
#' Use rqResMar_dist function for distance model.
#' As defined in Knape et al. 2018 and available for unmarked models in the nmixgof package.
#'
#' @param umFit an unmarked model fit to point count data.
#'
#' @return a matrix of residuals with the same dimensions as the observations (rows = cites, columns = replicates)
#'
#' @export
rqResMar = function(umFit) {

   # get the fitted values from the model (i.e. predicted parameter(s))
   fitval = fitted(umFit, K = umFit@K) # A bug in older versions of unmarked (fixed now) may cause incorrect fitted values for NB models unless K is supplied.

   # get the observed counts
   y = umFit@data@y

   # calculate the marginal rq residuals for each type of unmarked model
   # if the unmarked model was fit with a Poisson mixture, use ppois as the pFun here
   if (identical(umFit@mixture, "P")) {

      rqr = rqRes(y, pFun = ppois, lambda = fitval)

   } else
      # if unmarked model was fit with Negative Binomial mixture, use pnbinom as the pFun
      if (identical(umFit@mixture, "NB")) {
         if (!identical(umFit@estimates["alpha"]@invlink, "exp"))
            stop("Unknown link function.")
         # get the model-estimated size parameter
         size = exp(umFit@estimates["alpha"]@estimates)
         # calculate the rq residuals
         rqr = rqRes(y, pFun = pnbinom, mu = fitval, size = size)


      } else
         # if unmarked model was fit with Zero-Inflated Poisson mixture, use pZIP here
         if (identical(umFit@mixture, "ZIP")) {
            if (!identical(umFit@estimates["psi"]@invlink, "logistic"))
               stop("Unknown link function.")

            # define the CDF for the zero-inflated Poisson
            pZIP = function(y, lambda, psi) {
               psi*(y>=0) + (1-psi) * ppois(y, lambda/(1-psi))
            }

            # get the model-estimated ZIP parameter, psi
            psi = plogis(umFit@estimates["psi"]@estimates)

            # calculate the rq residuals

            rqr = rqRes(y, pFun = pZIP, lambda = fitval, psi = psi)
         } else {stop("Mixture not recognized.")}
   rqr
}

#' Site-Sum Randomized Quantile Residuals
#'
#' Calculates Site-Sum Randomized-quantile residuals for a point count model fit with unmarked.
#' Use rqResSum_dist function for distance model.
#' As defined in Knape et al. 2018 and available for unmarked models in the nmixgof package.
#'
#' @param umFit an unmarked model fit to point count data.
#'
#' @return a matrix of residuals with the same dimensions as the observations (rows = cites, columns = replicates)
#'
#' @importFrom unmarked getP
#'
#' @export
rqResSum = function(umFit, show.warnings = FALSE) {

   # get the predicted lambda value at each site (for the observed covariate values)
   lam = predict(umFit, type="state")[,1]

   # get the predicted detection probability values for every site-replicate combo
   p = getP(umFit)

   # get the observed counts
   if (length(umFit@sitesRemoved) > 0)
      y = umFit@data@y[-umFit@sitesRemoved,]
   else
      y = umFit@data@y

   # create an empty matrix for the counts at each site
   cumProb = matrix(0, nrow = nrow(y), ncol = 2)

   # pick the PDF/density function of the appropriate distribution (used to analyze the data)
   dfun = switch(umFit@mixture,
                 # Poisson distribution
                 P = function(N) {dpois(N, lam)},
                 # Negative Binomial
                 NB = function(N) {dnbinom(N, mu=lam, size=exp(coef(umFit, type="alpha")))},
                 # Zero-inflated Poisson
                 ZIP = function(N) {
                    psi = plogis(coef(umFit, type="psi"))
                    # adjust Poisson PDF to account for zero-inflation
                    (1-psi)*dpois(N, lam/(1-psi)) + psi*(N==0)
                 })

   # From Knape 2018:
   # "If the p_it are all the same, F_BinSum is simply the cumulative probability function of a binomial distribution with index TN but if the p_it are not all identical then F_BinSum is more complex."
   # The function Knape provides is for the more complex case:

   # loop over possible N values (i.e. k) to sum/integrate over them and calculate the F(z-1) and F(z)
   for (N in 0:umFit@K) {
      cumProb =
         cumProb + # start with previous value and add:
         # matrix with probability of that N value given each site's lambda value
         # this is the P_i(N) function in Knape
         kronecker(dfun(N), t(c(1,1))) * # multiply it by:
         # for each site, integrate the total probability of observing y_i individuals (site i) by adding across all possible y_it values (t = replicate) that could've generated that y_i value
         pbinsum(y, rep(N, nrow(y)), p)[,2:3]
      # columns 2 & 3 return (for each site)
      # 2 = total probability of observing up to (y_i - 1) individuals across T replicates
      # 3 = total probability of observing up to y_i individuals across T replicates
   }

   # calculate residuals
   res = rqRes0(cumProb[,1], cumProb[,2])


   # throw a warning if any residuals are infinite
   if (any(is.infinite(res)) & show.warnings) {
      warning(paste(sum(is.infinite(res)), " residuals infinite."))
   }

   # return the residuals
   res
}

#' Observation Randomized Quantile Residuals
#'
#' Calculates Observation Randomized-quantile residuals for a point count model fit with unmarked.
#' Use rqResSum_dist function for distance model.
#' As defined in Knape et al. 2018 and available for unmarked models in the nmixgof package.
#'
#' @param umFit an unmarked model fit to point count data.
#'
#' @return a matrix of residuals with the same dimensions as the observations (rows = cites, columns = replicates)
#'
#' @importFrom unmarked getP
#'
#' @export
rqResObs = function(umFit, show.warnings = FALSE) {

   # NA's for each site
   rN = integer(nrow(umFit@data@y)) + NA

   # Draw/sample a random N values (1 per site) from the estimated PDF of possible N values
   {
      # remove any sites that were removed from the unmarked analysis
      if (length(umFit@sitesRemoved) > 0){
         # sample the PDF of possible N values for each site.
         rN[-umFit@sitesRemoved] =
            apply(X = unmarked::ranef(umFit)@post[,,1],
                  # unmarked::ranef = "Estimate posterior distributions of the random variables (latent abundance or occurrence) using empirical Bayes methods. So this returns an estimated probability that the true (latent) N == each possible k value. Basically, a PDF for N across all possible k's. They take the first slice, which corresponds to the first sampling period and just ignore all other sampling periods." These are passed to the probability weights argument in the base::sample() function.
                  MARGIN = 1, # across all possible N's within a single site.
                  FUN = sample,
                  x = umFit@K + 1, # items from which to choose (can choose 0, so have to add 1). This is just the number or rows in X
                  size = 1,        # choose only a single value
                  replace = FALSE) - 1 # subtract 1 so that it includes 0 and doesn't go above K.
      } else {

         # draw a random N for each site (see comments above)
         rN = apply(X = unmarked::ranef(umFit)@post[,,1],
                    MARGIN = 1,
                    FUN = sample,
                    x = umFit@K + 1,
                    size = 1,
                    replace = FALSE) - 1
      }
   }

   # getP is a function from unmarked
   # returns a matrix of estimated detection probabilities for each site (rows) and replicate (columns).
   p = getP(umFit)
   # for my simulations, I can just use the same detection prob.
   # p = matrix(data = p_est, nrow = n_sites, ncol = n_samps)

   # get the residuals
   res = rqRes(umFit@data@y,
               pFun = pbinom,
               # Computes the generalised kronecker product of two arrays, X and Y
               # I think they're just using this to create a matrix out of the rN vector
               size = kronecker(rN, # vector of length n_sites
                                t(rep(1, ncol(p)))), # vector of 1's of length n_replicates
               prob=p)

   if (any(is.infinite(res)) & show.warnings) {
      #throw a warning if any residuals are infinite
      warning(paste(sum(is.infinite(res)), " residuals infinite."))
   }

   # return the residuals
   res
}

# Randomized Quantile Residuals for unmarked distance sampling
# I think I would need to do 1 of 2 things:
# 1. extend the observation rq residuals to a multinomial distribution.
# Difficulty: getting N (size of multinomial). That would require
# a crazy calculation of possibilities for each of many, many possible N's up to K.
# Benefit: could just draw a random distance from the normal distribution
# take 1 individual out of that distance category for the F(z-1) calculation.
# 2. use the generic rqres function (poisson) for each distance bin
# Difficulty: there's a F(z-1) calculation for each distance bin instead
# of 1 across all distance bins.
# Benefit: really easy

# I used option 2, but that means that the residuals for each distance
# bin are not connected, so I should probably look at each bin independently.
# function for calculating Marginal Randomized-quantile residuals

#' Marginal Randomized Quantile Residuals
#'
#' Calculates Marginal Randomized-quantile residuals for a distance model fit with unmarked.
#' Note that I currently calculate these residuals by subtracting 1 from each
#' distance bin. That's pretty hacky and likely just plain wrong. Don't use
#' residuals from more than 1 distance bin. Just use residuals from a single bin.
#' It would be better to extend the observation rq residuals to a multinomial
#' distribution, and use the pbinsum function to get the possible N's (up to K)
#' across all distance breaks. But then I could just draw a random distance from
#' the normal distribution, take 1 individual out of that distance category (for
#' the F(z-1) calculation) and calculate the residuals with those 2 values &
#' distributions.
#' Use rqResMar function for point-count model.
#' As defined in Knape et al. 2018 and available for unmarked models in the nmixgof package.
#'
#' @param umFit an unmarked model fit to point count data.
#'
#' @return a matrix of residuals with the same dimensions as the observations (rows = cites, columns = distance breaks)
#'
#' @export
rqResMar_dist = function(umFit) {
   # get the fitted values from the model (i.e. expected counts per cell)
   fitval = fitted(umFit)

   # get the observed counts
   y = umFit@data@y

   # get the model parameters
   sig <- umFit@estimates@estimates$det@estimates
   # distance breaks for the model
   db <- umFit@data@dist.breaks
   strip.widths <- diff(db)

   # calculate the marginal rq residuals for each type of unmarked model
   # distance sampling is fit with a multinomial distribution
   # calculate things entirely here instead of using the rqRes function
   {
      # y is the observed counts
      # calculate cell probabilities for the multinomial distribution
      {
         f.0 <- 2 * dnorm(0, 0, sd=sig)
         int <- 2 * (pnorm(db[-1], 0, sd=sig) -
                        pnorm(db[-(ncol(y)+1)], 0, sd=sig))
         pi  <- int / f.0 / strip.widths
      }
   }

   rqr <- y*NA
   for(d in 1:length(strip.widths)){
      rqr[,d] = rqRes(y[,d], pFun = ppois, lambda = fitval[,d])
   }
   rqr
}

#' Site-Sum Randomized Quantile Residuals
#'
#' Calculates Site-Sum Randomized-quantile residuals for a distance model fit with unmarked.
#' Note that I currently calculate these residuals by subtracting 1 from each
#' distance bin. That's pretty hacky and likely just plain wrong. Don't use
#' residuals from more than 1 distance bin. Just use residuals from a single bin.
#' It would be better to extend the observation rq residuals to a multinomial
#' distribution, and use the pbinsum function to get the possible N's (up to K)
#' across all distance breaks. But then I could just draw a random distance from
#' the normal distribution, take 1 individual out of that distance category (for
#' the F(z-1) calculation) and calculate the residuals with those 2 values &
#' distributions.
#' Use rqResSum function for point-count model.
#' As defined in Knape et al. 2018 and available for unmarked models in the nmixgof package.
#'
#' @param umFit an unmarked model fit to point count data.
#'
#' @return a matrix of residuals with the same dimensions as the observations (rows = cites, columns = distance breaks)
#'
#' @importFrom unmarked getP
#'
#' @export
rqResSum_dist = function(umFit, show.warnings = FALSE) {
   # I really don't think these are what I want, yet.

   # get the predicted lambda value at each site (for the observed covariate values)
   lam = predict(umFit, type="state")[,1]

   # get the predicted detection probability values for every site-distance-bin combo
   p = getP(umFit)

   # get the observed counts
   y = umFit@data@y

   # define K as the predicted N * 10
   K <- round(predict(umFit, type="state")[1,1]) * 10

   # create an empty matrix for the counts at each site
   cumProb = matrix(0, nrow = nrow(y), ncol = 2)

   # pick the PDF/density function of the appropriate distribution (used to analyze the data)
   dfun = switch('P', # always poisson...
                 # Poisson distribution
                 P = function(N) {dpois(N, lam)},
                 # Negative Binomial
                 NB = function(N) {dnbinom(N, mu=lam, size=exp(coef(umFit, type="alpha")))},
                 # Zero-inflated Poisson
                 ZIP = function(N) {
                    psi = plogis(coef(umFit, type="psi"))
                    # adjust Poisson PDF to account for zero-inflation
                    (1-psi)*dpois(N, lam/(1-psi)) + psi*(N==0)
                 })

   # From Knape 2018:
   # "If the p_it are all the same, F_BinSum is simply the cumulative
   # probability function of a binomial distribution with index TN but
   # if the p_it are not all identical then F_BinSum is more complex."
   # The function Knape provides is for the more complex case:

   # loop over possible N values (i.e. k) to sum/integrate over them and calculate the F(z-1) and F(z)
   for (N in 0:K) {
      cumProb =
         cumProb + # start with previous value and add:
         # matrix with probability of that N value given each site's lambda values
         # this is the P_i(N) function in Knape
         kronecker(dfun(N), t(c(1,1))) * # multiply it by:
         # for each site, integrate the total probability of observing y_i individuals (site i) by adding across all possible y_it values (t now = bin instead of replicate) that could've generated that y_i value
         pbinsum(y, rep(N, nrow(y)), p)[,2:3]
      # columns 2 & 3 return (for each site)
      # 2 = total probability of observing up to (y_i - 1) individuals across T bins
      # 3 = total probability of observing up to y_i individuals across T bins
   }

   # calculate residuals
   res = rqRes0(cumProb[,1], cumProb[,2])

   # throw a warning if any residuals are infinite
   if (any(is.infinite(res)) & show.warnings) {
      warning(paste(sum(is.infinite(res)), " residuals infinite."))
   }

   # return the residuals
   res
}

#' Observation Randomized Quantile Residuals
#'
#' Calculates Observation Randomized-quantile residuals for a distance model fit with unmarked.
#' Note that I currently calculate these residuals by subtracting 1 from each
#' distance bin. That's pretty hacky and likely just plain wrong. Don't use
#' residuals from more than 1 distance bin. Just use residuals from a single bin.
#' It would be better to extend the observation rq residuals to a multinomial
#' distribution, and use the pbinsum function to get the possible N's (up to K)
#' across all distance breaks. But then I could just draw a random distance from
#' the normal distribution, take 1 individual out of that distance category (for
#' the F(z-1) calculation) and calculate the residuals with those 2 values &
#' distributions.
#' Use rqResMar function for point-count model.
#' As defined in Knape et al. 2018 and available for unmarked models in the nmixgof package.
#'
#' @param umFit an unmarked model fit to point count data.
#'
#' @return a matrix of residuals with the same dimensions as the observations (rows = cites, columns = distance breaks)
#'
#' @export
rqResObs_dist = function(umFit, show.warnings = FALSE) {

   # NA's for each site
   rN = integer(nrow(umFit@data@y)) + NA
   # I'll set K = predicted * 10
   K = round(predict(umFit, type="state")[1,1]) * 10

   # Draw/sample a random N values (1 per site) from the estimated PDF of possible N values
   {
      rN = apply(X = unmarked::ranef(umFit, K = K)@post[,,1],
                 MARGIN = 1,
                 FUN = sample,
                 x = K + 1,
                 size = 1,
                 replace = FALSE) - 1
   }

   # getP is a function from unmarked
   # returns a matrix of estimated detection probabilities for each site (rows) and replicate (columns).
   p = getP(umFit)

   # get the residuals
   res = rqRes(umFit@data@y,
               pFun = pbinom,
               # Computes the generalised kronecker product of two arrays, X and Y
               # I think they're just using this to create a matrix out of the rN vector
               size = kronecker(rN, # vector of length n_sites
                                t(rep(1, ncol(p)))), # vector of 1's of length n_replicates
               prob=p)

   if (any(is.infinite(res)) & show.warnings) {
      #throw a warning if any residuals are infinite
      warning(paste(sum(is.infinite(res)), " residuals infinite."))
   }

   # return the residuals
   res
}
