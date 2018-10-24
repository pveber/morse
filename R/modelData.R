#' Create a list giving data to use in Bayesian modelling
#' 
#' @param x An object of class \code{survData}
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @return A list for parameterization of priors for Bayesian modelling
#' 
#' @importFrom dplyr group_indices_
#' 
modelData <- function(x, ...){
  UseMethod("modelData")
}
