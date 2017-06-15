#' Fitting models for survival analysis using Bayesian inference
#'
#' Survival analysis using Bayesian inference for \code{reproFitTT}, \code{survFitTT},
#' \code{survFitTKTD} and \code{gm_survFitTKTD} objects.
#'
#' @param x an object used to select a method
#' @param \dots Further arguments to be passed to generic methods

#' @export
survFit <- function(x, ...){
  UseMethod("survFit")
}
