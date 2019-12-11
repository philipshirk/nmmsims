#' Iteratively find sigma for half-normal detection function
#'
#' Because I'm too lazy to work out the math for calculating sigma (SD) precisely,
#' I created this function to iteratively find sigma (SD) for a half-normal 
#' detection function. It assumes a half-normal detection function with SD = sigma.
#'
#' @param start starting guess (default is fine)
#' @param W known transect half-width
#' @param p_true known detection probability
#'
#' @return sd of normal distribution to the nearest 0.01
#'
#' @examples
#' find_sigma()
#'
#' @export
find_sigma <- function(start = 0.1, W = 20, p_true=0.42){
   if(p_true <=0 | p_true >=1) stop("p_true must be between 0 and 1")
   if(start > .1) warning("this is a lazily-created function. If you start with a large starting guess, it may not function properly")
   
   sigma_guess <- start
   
   # calculate the average detection probability at the given sigma_guess
   pa_temp <- calc_det_prob(W = W, sigma = sigma_guess)
   
   # while pa_temp is still less than the true average detection probability, keep adding on
   while(pa_temp < (p_true)){
      sigma_guess <- sigma_guess + 0.1
      pa_temp <- calc_det_prob(W = W, sigma = sigma_guess)
   }
   # go down to just below proper sigma
   sigma_guess <- sigma_guess - 0.1
   pa_temp <- calc_det_prob(W = W, sigma = sigma_guess)
   
   # then loop through with a smaller increment
   while(pa_temp < (p_true)){
      sigma_guess <- sigma_guess + 0.01
      pa_temp <- calc_det_prob(W = W, sigma = sigma_guess)
   }
   return(sigma_guess)
}