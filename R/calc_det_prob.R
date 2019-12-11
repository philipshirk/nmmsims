#' Calculate Mean Detection Probability
#' 
#' This uses sigma (the SD of a normal distribution) and W (the transect half-width)
#' to calculate the mean detection probability. Only supports half-normal 
#' detection function at this point.
#' 
#' @param W transect half-width (in meters)
#' @param sigma Known SD of normal distribution describing the detection function (in meters)
#' 
#' @return number between 0 and 1
#' 
#' @export
calc_det_prob <- function(W, sigma){
   (pnorm(q = W, mean = 0, sd = sigma) - 0.5) / (dnorm(x = 0, mean = 0, sd = sigma)*W)
}