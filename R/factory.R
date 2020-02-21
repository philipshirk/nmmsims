#' Catch errors and warnings and store them for subsequent evaluation
#' from https://stackoverflow.com/a/29465795/3174566
#' 
#' Factory modified from a version written by Martin Morgan on Stack Overflow (see below).  
#' Factory generates a function which is appropriately wrapped by error handlers.  
#' If there are no errors and no warnings, the result is provided.  
#' If there are warnings but no errors, the result is provided with a warn attribute set.
#' If there are errors, the result returns a list with the elements of warn and err.
#' This is a nice way to recover from a problem that may have occurred during loop evaluation or during cluster usage.
#' Check the references for additional related functions.
#' I have not included the other factory functions included in the original Stack Overflow answer because they did not play well with the return item as an S4 object.
#' @export
#' @param fun The function to be turned into a factory
#' @return The result of the function given to turn into a "factory".  If this function was in error "An error has occurred" as a character element.  factory-error and factory-warning attributes may also be set as appropriate.
#' @references
#' \url{http://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function}
#' @author Martin Morgan; Modified by Russell S. Pierce
#' @import data.table
#' @examples 
#' f.log <- factory(log)
#' f.log("a")
#' f.as.numeric <- factory(as.numeric)
#' f.as.numeric(c("a","b",1))
factory <- function (fun) {
  errorOccurred <- FALSE
  
  function(...) {
    # set warnings and errors to NULL
    warn <- err <- NULL
    
    # with withCallingHandlers to try to run the function
    res <- withCallingHandlers(
      tryCatch(expr = fun(...), 
               error = function(e) {
                 err <<- conditionMessage(e)
                 errorOccurred <<- TRUE
                 NULL
               }),
      warning = function(w) {
        warn <<- append(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    
    # why not just put this in the error section of the trycatch above?
    if (errorOccurred) {
      res <- "An error occurred in the factory function"
    } 
    
    # set the attributes of the results
    if (is.character(warn)) {
      data.table::setattr(res,"factory-warning",warn)
    } else {
      data.table::setattr(res,"factory-warning",NULL) 
    }
    
    if (is.character(err)) {
      data.table::setattr(res,"factory-error",err)
    } else {
      data.table::setattr(res, "factory-error", NULL)
    }  
    return(res)
  }
}

#' an internal helper function for "factory()" function. It checks to see if a particular attribute is present
#' @export
#' @param x The result returned by the factory function
#' @param what the attribute to return from the result returned by the factory function
#' @return logical
.has <- function(x, what) {
  # check to make sure that attributes are present
  !is.null(attr(x,what))
}

#' Check if the "factory()" function produced a warning
#' @export
#' @param x The result returned by the factory function
#' @return logical
hasWarning <- function(x) .has(x, "factory-warning")

#' Check if the "factory()" function produced an error
#' @export
#' @param x The result returned by the factory function
#' @return logical
hasError <- function(x) .has(x, "factory-error")

#' Check if the "factory()" function produced a warning or an error
#' @export
#' @param x The result returned by the factory function
#' @return logical
isClean <- function(x) {!(hasError(x) | hasWarning(x))}