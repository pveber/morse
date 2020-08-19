#' Create a list giving data to use in Bayesian inference.
#' 
#' @param x An object of class \code{survData}
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @return A list for parameterization of priors for Bayesian inference.
#' 
#' 
modelData <- function(x, ...){
  UseMethod("modelData")
}
