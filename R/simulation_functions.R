#' Primary function for simulating & analysing datasets that are biased by over-dispersion in the number of observations (via double-counting) (Scenario 1)
#'
#' This function is one of the master functions in this package. 
#' It's the same as "base_simulate_function_nobs", but wrapped in the "factory"
#' function for (hopefully) improved condition handling. It simulates and
#' analyzes a single dataset that is over-dispersed in observations, n_obs. 
#' To do this many times, use an apply function. It's currently simplified and 
#' assumes constant abundance, detection, and transect length. It does not (yet) 
#' simulate goodness-of-fit metrics. It simulates both point count and distance 
#' data and analyses both datasets using both unmarked and optim.
#'
#' @param n_sites number of sites (transects)
#' @param n_samps number of samples (replicates) per site
#' @param lambda  mean abundance at every site (single draw per site that stays constant across samples)
#' @param alpha   probability of NOT being observed (i.e. alpha = 1 - detection probability)
#' @param gamma   probability of being observed exactly twice (gamma = 1 - alpha - beta)
#' @param beta    (optional) probability of being observed exactly once. If not specified, it's calculated as: beta = 1 - alpha - gamma
#' @param det_prob (optional) probability of being observed. If not specified, it's calculated as: det_prob = beta + gamma
#' @param sigma (optional) detection parameter (meters). Best to leave this blank and it will be calculated from det_prob.
#' @param W transect half-width (meters)
#' @param reps_to_analyze the number of samples/replicates to analyze. If NA, it will analyse all replicates in the data.
#' @param return What to return from the function call. Currently the only option is 'results'. May change this to only analyze simulated goodness-of-fit metrics.
#' @param savefilename The simulated datasets and results ARE saved to file (currently not optional). This provides the path and filename for saving the intermediate steps in the analysis. 
#'
#' @return if everything works well, it returns a data.frame with the results of simulating a single dataset, analyzing it in 4 ways, and calculating randomized quantile residuals a la Knape et al. 2018. It also saves a list with the simulated dataset, dataframe of results (minus rqr residual info), and the actual rqr residuals to a savefilename inside the folder path 'working directory'/results/Scenario 1/savefilename. If there is an error, the function returns NA and also saves a file to 'working directory'/Scenario 1/set x/errors with the simulated dataset and the results data.frame but no rq-residuals (it's basically assumed that the rq-residuals were the source of the error.) Similarly, with a warning the function returns the results data.frame and also saves a file to 'working directory'/Scenario 1/set x/warnings with the simulated dataset and the results data.frame but no rq-residuals (it's basically assumed that the rq-residuals were the source of the error.) The user will have to go back and try to calculate rq-residuals from the output later. 
#'
#' @examples
#' simulation_function_nobs()
#'
#' @import tidyverse
#' 
#' @export
simulation_function_nobs <- factory(base_simulation_function_nobs)


#' Primary function for simulating & analysing datasets that are biased by over-dispersion in N (Scenario 2)
#'
#' This function is one of the master functions in this package. 
#' It's the same as "base_simulate_function_N", but wrapped in the "factory" 
#' function for (hopefully) improved condition handling. It simulates and
#' analyzes a single dataset that is over-dispersed in true abundance, N. 
#' To do this many times, use an apply function. It's currently simplified and 
#' assumes constant abundance, detection, transect length, NO simulated goodness-of-fit
#' metrics. It simulates both point count and distance data and analyses both datasets
#' using both unmarked and optim.
#'
#' @param n_sites number of sites (transects)
#' @param n_samps number of samples (replicates) per site
#' @param lambda1  partial mean abundance at every site (single draw per site that stays constant across samples)
#' @param lambda2  partial mean abundance at every replicate (new draw for every sample at every site)
#' @param det_prob mean detection probability
#' @param sigma detection parameter (meters). Just leave this blank and it will be calculated from det_prob
#' @param W transect half-width (meters)
#' @param reps_to_analyze the number of samples/replicates to analyze. If NA, it will analyse all replicates in the data.
#' @param return What to return from the function call. Currently the only option is 'results'. May change this to only analyze simulated goodness-of-fit metrics.
#' @param savefilename The simulated datasets and results ARE saved to file (currently not optional). This provides the path and filename for saving the intermediate steps in the analysis. 
#'
#' @return if everything works well, it returns a data.frame with the results of simulating a single dataset, analyzing it in 4 ways, and calculating randomized quantile residuals a la Knape et al. 2018. It also saves a list with the simulated dataset, dataframe of results (minus rqr residual info), and the actual rqr residuals to a savefilename inside the folder path 'working directory'/results/Scenario 2/savefilename. If there is an error, the function returns NA and also saves a file to 'working directory'/Scenario 2/set x/errors with the simulated dataset and the results data.frame but no rq-residuals (it's basically assumed that the rq-residuals were the source of the error.) Similarly, with a warning the function returns the results data.frame and also saves a file to 'working directory'/Scenario 2/set x/warnings with the simulated dataset and the results data.frame but no rq-residuals (it's basically assumed that the rq-residuals were the source of the error.) The user will have to go back and try to calculate rq-residuals from the output later. 
#'
#' @examples
#' simulation_function_N()
#'
#' @import tidyverse
#' 
#' @export
simulation_function_N <- factory(base_simulation_function_N)



#' Primary function for simulating & analysing datasets that are biased by over-dispersion in the detection probability (Scenario 3)
#'
#' This function is one of the master functions in this package. 
#' It's the same as "base_simulate_function_p", but wrapped in the "factory" 
#' function for (hopefully) improved condition handling. It simulates and
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
#' simulation_function_p()
#'
#' @import tidyverse
#' 
#' @export
simulation_function_p <- factory(base_simulation_function_p)