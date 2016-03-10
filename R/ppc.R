#' Posterior predictive check plot
#' 
#' Plots posterior predictive check for \code{reproFitTT}, \code{survFitTT} and
#' \code{survFitTKTD} objects.
#' 
#' @param x an object used to select a method
#' @param \dots Further arguments to be passed to generic methods

#' @export
ppc <- function(x, ...){
  UseMethod("ppc")
}
