#' Creates a dataset for survival analysis
#'
#' Dataset for survival analysis using Bayesian inference for \code{survFitTT},
#' \code{survFitTKTD} and \code{gm_survFitTKTD} objects.
#'
#' @param x an object used to select a method
#' @param \dots Further arguments to be passed to generic methods

#' @export
survData <- function(x, ...){
  UseMethod("survData")
}
