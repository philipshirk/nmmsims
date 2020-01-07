#' Test for Normality
#'
#' Run 3 different tests of normality, including shapiro.wilks (only if length(x) < 5000), 
#' AD test from package nortest, and lillie from nortest
#'
#' @param x vector of values to be tested for normality
#' @param name name for saving the results iff there's an error
#'
#' @return If everything runs smoothly, this will return a vector with 3 values for each of the 3 tests of normality.
#'
#' @examples
#' normtest(x = rnorm(100), 'name')
#'
#' @importFrom nortest ad.test lillie.test
#' 
#' @export
normtest <- function(x, name){
   out <- tryCatch(expr = {
      # if there are infinite values, just return infinity
      if(any(is.infinite(x))){
         list(p.value = x[which(is.infinite(x))][1])
      } else if(length(x) < 5000){
         # otherwise, calculate the shapiro test results
         st <- shapiro.test(x)
         ad <- nortest::ad.test(x)
         li <- nortest::lillie.test(x)
         
         c(st[['p.value']], ad[['p.value']], li[['p.value']])
      } else {
         ad <- nortest::ad.test(x)
         li <- nortest::lillie.test(x)
         
         c(NA, ad[['p.value']], li[['p.value']])
      }
   }, 
   error = function(cond){
      # if there's an error, save a file to help ID the error later on
      saveRDS(object = list(error = cond,
                            resids = x),
              file = paste0('./errors/scenario_02_set_1_Link_pars_loop_normtest_', 
                            gsub(pattern = ' ', replacement = '_', 
                                 x = gsub(pattern = '[[:punct:]]', 
                                          replacement = '',x = Sys.time())),
                            '_',
                            paste(Sys.info()[['nodename']], Sys.getpid(), sep='_'),
                            '_',
                            name,
                            '.RDS'))
      return(list(p.value = NA))
   },
   warning = function(cond){
      # if there's a warning, save a file to help ID the warning later on
      saveRDS(object = list(warning = cond,
                            resids = x),
              file = paste0('./warnings/scenario_02_set_1_Link_pars_normtest_', 
                            gsub(pattern = ' ', replacement = '_', 
                                 x = gsub(pattern = '[[:punct:]]', 
                                          replacement = '',x = Sys.time())),
                            '_',
                            paste(Sys.info()[['nodename']], Sys.getpid(), sep='_'),
                            '_',
                            name,
                            '.RDS'))
      return(list(p.value = NA))
   })
   return(out)
}