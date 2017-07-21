#' Create a list giving Data to use in Bayesian modelling (JAGS or Stan)
#'
#' @param x An object of class 'survData'
#' @param \dots Further arguments to be passed to generic methods
#'
#' @return A list for parameterization of priors for Bayesian modelling
#'
#' @importFrom dplyr group_indices_
#'
#' @export
#' 
modelData <- function(x, ...){
  UseMethod("modelData")
}
